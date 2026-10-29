function S = ieeglab_spectral(EEG, varargin)
% ieeglab_spectral() - Frequency-domain analysis of iEEG: power spectrum of each
%                      contact, or time-frequency power with a baseline and the
%                      power of a band (high-frequency activity by default).
%                      Compute only; ieeglab_plot_spectral draws the result.
%
% Usage:
%   S = ieeglab_spectral(EEG, 'type', 'spectrum')               % Welch, 1-200 Hz
%   S = ieeglab_spectral(EEG, 'type', 'timefreq')               % 5:5:200 Hz, 0.2 s windows
%   S = ieeglab_spectral(EEG, 'type', 'timefreq', 'freqs', 4:2:40, 'cycles', 5, ...
%                        'baseline', [-500 -100], 'norm', 'relchange', 'band', [8 12])
%
% Options (name/value):
%   'type'        'spectrum' - power spectral density of each contact (Welch:
%                              Hann windows, 50% overlap, averaged). Continuous
%                              or epoched data. Units uV^2/Hz.
%                 'timefreq' - power of each contact over time and frequency,
%                              averaged over trials. Epoched data. Default.
%   'freqs'       spectrum: [low high] range, default [1 200]. Time-frequency:
%                 the frequencies, default 5:5:200 (FieldTrip's intracranial
%                 tutorial, Stolk et al., 2018). Frequencies at or above the
%                 Nyquist frequency are dropped.
%   'window'      seconds. Spectrum: segment length, default 1. Time-frequency:
%                 fixed Hann window, default 0.2 (as the FieldTrip tutorial).
%   'cycles'      time-frequency: a window of this many cycles at each
%                 frequency instead of a fixed one (wavelet-like). Default [].
%   'timestep'    time-frequency: ms between output times. Default 10.
%   'baseline'    time-frequency: [start end] ms. Only times whose whole window
%                 lies inside it are used, so no post-stimulus sample enters the
%                 baseline. Default [epoch start, 0].
%   'norm'        time-frequency: 'db' (10 log10 power/baseline, default),
%                 'relchange' (power/baseline - 1, as FieldTrip's tutorial),
%                 'ratio', or 'none'.
%   'band'        [low high] Hz averaged into one band. Default [70 150], high-
%                 frequency activity (HFA). Spectrum: band power in uV^2.
%   'band_window' time-frequency: [start end] ms averaged into one value per
%                 contact (for the brain map). Default [0, epoch end].
%   'condition'   only the epochs of this condition (event type at time 0, or
%                 stimulation site for CCEP data). Default all.
%   'stim_exclude' spectrum of continuous CCEP data: [start end] ms around each
%                 pulse that no segment may overlap, so the spectrum is that of
%                 the ongoing activity, not of the artifact and the responses.
%                 Default [-50 1000]; [] keeps every segment.
%   'artifact_window' time-frequency of CCEP data: [start end] ms around the
%                 pulse; times whose window overlaps it are NaN (grey), since a
%                 window that contains the artifact shows power at every
%                 frequency. Default [-5 15]; [] keeps them.
%   'keeptrials'  time-frequency: also return single-trial power. Default false.
%   'max_mb'      memory per block of epochs, MB. Default 512.
%   'verbose'     default true
%
% On CCEP data each contact's own stimulation epochs are left out of its
% average (its artifact would dominate), as N1 and CRP do; S.n_trials gives the
% epochs each contact averages.
%
% Time-frequency power is the squared magnitude of the data convolved with
% Hann-tapered complex sinusoids (FieldTrip's mtmconvol with a Hann taper),
% scaled so that a sinusoid of amplitude A has power A^2 (uV^2). The data are
% transformed once per block of epochs and reused for every frequency. Times
% whose window does not fit inside the epoch are NaN.
%
% Output S: type, freqs, labels, bad, n_trials, condition, power (contacts x
% freqs, or contacts x freqs x times), units, band, band_value (one value per
% contact), and for time-frequency also times (ms), tf (normalized power),
% norm, baseline, band_tc (normalized band power, contacts x times),
% band_window, window_s (per frequency) and, with keeptrials, trials
% (contacts x freqs x times x epochs, single, every epoch) and use (contacts x
% epochs, the epochs each contact averages).
%
% Cedric Cannard, iEEGLAB, 2026

p = inputParser;
p.addParameter('type', 'timefreq', @(x) any(strcmpi(x, {'spectrum', 'timefreq'})));
p.addParameter('freqs', []);
p.addParameter('window', []);
p.addParameter('cycles', []);
p.addParameter('timestep', 10);
p.addParameter('baseline', []);
p.addParameter('norm', 'db', @(x) any(strcmpi(x, {'db', 'relchange', 'ratio', 'none'})));
p.addParameter('band', [70 150]);
p.addParameter('band_window', []);
p.addParameter('condition', '');
p.addParameter('stim_exclude', [-50 1000]);
p.addParameter('artifact_window', [-5 15]);
p.addParameter('keeptrials', false);
p.addParameter('max_mb', 512);
p.addParameter('verbose', true);
p.parse(varargin{:});
o = p.Results;
o.type = lower(char(o.type)); o.norm = lower(char(o.norm));

fs = double(EEG.srate);
C = EEG.nbchan;
labels = {EEG.chanlocs.labels};
bad = false(C, 1);
if isfield(EEG.chanlocs, 'status')
    bad = cellfun(@(x) ~isempty(x) && strcmpi(char(string(x)), 'bad'), {EEG.chanlocs.status})';
end
epoched = EEG.trials > 1 || (isfield(EEG, 'epoch') && ~isempty(EEG.epoch));
if strcmp(o.type, 'timefreq') && ~epoched
    error('ieeglab_spectral:needsEpochs', ...
        'Time-frequency power needs epoched data. For continuous data use ''type'', ''spectrum''.');
end

% ---------------- epochs used, and which each contact averages ----------------
N = size(EEG.data, 3);
use = true(C, N);                                   % use(c, n): contact c averages epoch n
condition = string(o.condition);
isCCEP = strcmp(ieeglab_detect_mode(EEG), 'ccep');
if epoched
    [sites, stimIdx] = ieeglab_epoch_sites(EEG);
    if condition ~= ""
        want = condition;
        if isCCEP, want = ieeglab_canonical_site({char(condition)}, labels); end
        keep = sites == want;
        if ~any(keep)
            error('ieeglab_spectral:noCondition', 'No epoch of condition "%s".', condition);
        end
        use(:, ~keep) = false;
    end
    if isCCEP
        for n = 1:N, use(stimIdx{n}, n) = false; end
    end
end
nTrials = sum(use, 2);
if o.verbose && epoched && any(nTrials < N)
    fprintf('[spectral] %d epochs; each contact averages %d to %d of them (condition, own stimulation left out).\n', ...
        N, min(nTrials), max(nTrials));
end

S = struct('type', o.type, 'labels', {labels}, 'bad', bad, 'n_trials', nTrials, ...
    'condition', char(condition), 'srate', fs, 'ccep', isCCEP);
if strcmp(o.type, 'spectrum')
    S = local_spectrum(EEG, o, S, use, epoched);
else
    S = local_timefreq(EEG, o, S, use);
end
end

% =============================== spectrum ===============================

function S = local_spectrum(EEG, o, S, use, epoched)
fs = S.srate;
if isempty(o.freqs), o.freqs = [1 200]; end
if isempty(o.window), o.window = 1; end
[C, T, N] = size(EEG.data);
L = round(o.window * fs);
if L > T
    L = T;
    warning('ieeglab_spectral:shortData', 'Data shorter than the window: one segment of %.3f s.', T / fs);
end
step = max(1, round(L / 2));
starts = 1:step:(T - L + 1);
if ~epoched
    % continuous data: no segment across a discontinuity (EEGLAB 'boundary' event)
    bnd = [];
    if isfield(EEG, 'event') && ~isempty(EEG.event)
        isB = strcmpi(cellfun(@(x) char(string(x)), {EEG.event.type}, 'UniformOutput', false), 'boundary');
        bnd = [EEG.event(isB).latency];
    end
    for b = bnd(:)'
        starts = starts(~(starts < b & starts + L - 1 > b));
    end
    if isempty(starts), error('ieeglab_spectral:noSegments', 'No segment without a discontinuity.'); end
    if S.ccep && ~isempty(o.stim_exclude)
        % nor across a pulse, its artifact or its response
        [~, si] = ieeglab_canonical_site(cellfun(@(x) char(string(x)), {EEG.event.type}, 'UniformOutput', false), S.labels);
        lat = [EEG.event(~cellfun(@isempty, si)).latency];
        w = round(double(o.stim_exclude) / 1000 * fs);
        n0 = numel(starts);
        for b = lat(:)'
            starts = starts(starts + L - 1 < b + w(1) | starts > b + w(2));
        end
        if isempty(starts)
            error('ieeglab_spectral:noSegments', ['No %.2f s segment fits between the pulses outside %s ms around ' ...
                'each; shorten ''window'' or ''stim_exclude''.'], L / fs, mat2str(o.stim_exclude));
        end
        if o.verbose
            fprintf('[spectral] %d of %d segments are away from the %d pulses (%s ms around each left out).\n', ...
                numel(starts), n0, numel(lat), mat2str(o.stim_exclude));
        end
    end
end
w = 0.5 - 0.5 * cos(2 * pi * (0:L-1)' / L);                % periodic Hann
nf = floor(L / 2) + 1;
scale = 2 * ones(nf, 1) / (fs * sum(w.^2));                 % one-sided PSD
scale(1) = scale(1) / 2;
if mod(L, 2) == 0, scale(end) = scale(end) / 2; end         % Nyquist bin

% Psum: (n_freq x n_chan x n_epochs), summed over segments. Blocks of segments:
% one FFT call for every segment, contact and epoch of the block.
nS = numel(starts);
nb = max(1, floor(o.max_mb * 2^20 / (2 * 16 * L * C * N)));
Psum = zeros(nf, C, N);
for b0 = 1:nb:nS
    sb = starts(b0:min(nS, b0 + nb - 1));
    I = (0:L-1)' + sb;                                       % (L x n_seg) sample indices
    seg = reshape(double(EEG.data(:, I(:), :)), C, L, numel(sb), N);
    seg = reshape(permute(seg, [2 1 3 4]), L, []);           % (L x C*n_seg*N)
    seg = (seg - mean(seg, 1)) .* w;
    Fs = fft(seg);
    Psum = Psum + reshape(sum(reshape(abs(Fs(1:nf, :)).^2, nf, C, numel(sb), N), 3), nf, C, N);
end
Pep = Psum .* scale / nS;                                    % (n_freq x n_chan x n_epochs)
U = reshape(double(use), 1, C, N);                           % use: (n_chan x n_epochs)
P = sum(Pep .* U, 3) ./ reshape(S.n_trials, 1, C);           % (n_freq x n_chan)
P(:, S.n_trials == 0) = NaN;

f = (0:nf-1)' * fs / L;
keepF = f >= min(o.freqs) & f <= max(o.freqs);
S.freqs = f(keepF)';
S.power = P(keepF, :)';                                     % (n_chan x n_freq)
S.units = 'uV^2/Hz';
S.window_s = L / fs;
S.n_segments = numel(starts);
S.band = o.band;
S.band_value = nan(C, 1);
if ~isempty(o.band)
    inB = f >= o.band(1) & f <= o.band(2);
    if any(inB)
        S.band_value = sum(P(inB, :), 1)' * fs / L;         % band power, uV^2
    elseif o.verbose
        warning('ieeglab_spectral:bandOutside', 'Band %s Hz is outside 0-%g Hz.', mat2str(o.band), fs / 2);
    end
end
end

% ============================= time-frequency =============================

function S = local_timefreq(EEG, o, S, use)
fs = S.srate;
if isempty(o.freqs), o.freqs = 5:5:200; end
if isempty(o.window), o.window = 0.2; end
[C, T, N] = size(EEG.data);
times = double(EEG.times(:))';

freqs = double(o.freqs(:))';
over = freqs >= fs / 2 | freqs <= 0;
if any(over)
    if o.verbose
        warning('ieeglab_spectral:aboveNyquist', 'Frequencies dropped (at or above %g Hz, or not positive): %s', ...
            fs / 2, mat2str(freqs(over)));
    end
    freqs = freqs(~over);
end
if isempty(freqs), error('ieeglab_spectral:noFreqs', 'No frequency below the Nyquist frequency (%g Hz).', fs / 2); end
F = numel(freqs);

% window per frequency: odd number of samples, centred on its output time
if isempty(o.cycles), len = repmat(o.window, 1, F); else, len = o.cycles ./ freqs; end
hL = round(len * fs / 2);                                   % half-length, samples
Lw = 2 * hL + 1;

% output times, on samples
toi = times(1):o.timestep:times(end);
sIdx = unique(min(T, max(1, round((toi - times(1)) / 1000 * fs) + 1)));
toi = times(sIdx);
Nt = numel(sIdx);
valid = (sIdx - hL(:) >= 1) & (sIdx + hL(:) <= T);         % (n_freq x n_times)
if ~any(valid(:))
    error('ieeglab_spectral:epochTooShort', ...
        'The epochs (%.3f s) are shorter than every window; shorten ''window'' or ''cycles''.', T / fs);
end

% kernels: Hann-tapered complex sinusoids, FFT once
nfft = local_fastlen(T + max(Lw) - 1);
Hf = zeros(nfft, F);
gain = zeros(1, F);
for k = 1:F
    w = 0.5 - 0.5 * cos(2 * pi * (1:Lw(k))' / (Lw(k) + 1));   % symmetric Hann (MATLAB hanning)
    tk = (-hL(k):hL(k))' / fs;
    Hf(:, k) = fft(w .* exp(1i * 2 * pi * freqs(k) * tk), nfft);
    gain(k) = 2 / sum(w);                                    % amplitude A -> |coefficient| A
end

% blocks of epochs: nfft x (contacts x epochs) complex, three times (transform,
% product, inverse)
perEpoch = 3 * 16 * nfft * C;
nB = max(1, floor(o.max_mb * 2^20 / perEpoch));
Psum = zeros(C, F, Nt);
if o.keeptrials, Ptr = nan(C, F, Nt, N, 'single'); end
for b0 = 1:nB:N
    tr = b0:min(N, b0 + nB - 1);
    n = numel(tr);
    Z = reshape(permute(double(EEG.data(:, :, tr)), [2 1 3]), T, C * n);   % (T x C*n)
    Z = Z - mean(Z, 1);
    Zf = fft(Z, nfft);
    U = reshape(double(use(:, tr)), 1, C * n);
    for k = 1:F
        v = find(valid(k, :));
        if isempty(v), continue; end
        Y = ifft(Zf .* Hf(:, k));
        cf = Y(sIdx(v) + hL(k), :) * gain(k);                % (n_valid x C*n)
        pw = real(cf).^2 + imag(cf).^2;
        if o.keeptrials
            Ptr(:, k, v, tr) = single(permute(reshape(pw, numel(v), C, n), [2 4 1 3]));
        end
        pw = reshape(pw .* U, numel(v), C, n);
        Psum(:, k, v) = Psum(:, k, v) + reshape(sum(pw, 3)', C, 1, numel(v));
    end
end
% CCEP: a window that contains the stimulation artifact shows power at every
% frequency; those times are left out (NaN), as times outside the epoch are
halfMs = hL / fs * 1000;
if S.ccep && ~isempty(o.artifact_window)
    hitsArtifact = toi - halfMs(:) <= o.artifact_window(2) & toi + halfMs(:) >= o.artifact_window(1);
    valid = valid & ~hitsArtifact;
    S.artifact_window = o.artifact_window;
end
P = Psum ./ S.n_trials;                                     % (n_chan x n_freq x n_times)
P(S.n_trials == 0, :, :) = NaN;
P(repmat(reshape(~valid, 1, F, Nt), C, 1, 1)) = NaN;
if o.keeptrials, Ptr(repmat(reshape(~valid, 1, F, Nt), C, 1, 1, N)) = NaN; end

% ---------------- baseline normalization ----------------
if isempty(o.baseline), o.baseline = [times(1) 0]; end
inBase = toi >= o.baseline(1) + halfMs(:) & toi <= o.baseline(2) - halfMs(:) & valid;   % (n_freq x n_times)
noBase = ~any(inBase, 2);
if any(noBase) && ~strcmp(o.norm, 'none') && o.verbose
    warning('ieeglab_spectral:noBaseline', ...
        'No window fits inside the baseline %s ms at %s Hz: those frequencies are NaN after normalization.', ...
        mat2str(o.baseline), mat2str(freqs(noBase)));
end
B = nan(C, F);
for k = find(~noBase)'
    B(:, k) = mean(P(:, k, inBase(k, :)), 3, 'omitnan');
end
switch o.norm
    case 'db',        TF = 10 * log10(P ./ B);
    case 'relchange', TF = P ./ B - 1;
    case 'ratio',     TF = P ./ B;
    case 'none',      TF = P;
end

S.freqs = freqs;
S.times = toi;
S.power = P;
S.tf = TF;
S.units = 'uV^2';
S.norm = o.norm;
S.baseline = o.baseline;
S.window_s = Lw / fs;
if o.keeptrials
    S.trials = Ptr;          % every epoch, including those a contact does not average
    S.use = use;             % (contacts x epochs): the epochs each contact averages
end

% ---------------- band ----------------
S.band = o.band;
if isempty(o.band_window), o.band_window = [0 times(end)]; end
S.band_window = o.band_window;
S.band_tc = nan(C, Nt);
S.band_value = nan(C, 1);
if ~isempty(o.band)
    inB = freqs >= o.band(1) & freqs <= o.band(2);
    if any(inB)
        S.band_tc = reshape(mean(TF(:, inB, :), 2, 'omitnan'), C, Nt);
        inW = toi >= o.band_window(1) & toi <= o.band_window(2);
        S.band_value = mean(S.band_tc(:, inW), 2, 'omitnan');
    elseif o.verbose
        warning('ieeglab_spectral:bandOutside', 'No analysed frequency in the band %s Hz.', mat2str(o.band));
    end
end
end

function n = local_fastlen(n)
% Smallest length >= n with no prime factor above 7 (fast FFT sizes)
while any(factor(n) > 7), n = n + 1; end
end
