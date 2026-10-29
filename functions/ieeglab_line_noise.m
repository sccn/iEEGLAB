function [EEG, info] = ieeglab_line_noise(EEG, opt)
% ieeglab_line_noise() - Remove power-line noise by subtracting fitted sines,
%                        without filtering, so nothing rings around the pulses.
%
% Usage:
%   EEG         = ieeglab_line_noise(EEG, struct('freqs', 60))
%   [EEG, info] = ieeglab_line_noise(EEG, opt)
%
% Options:
%   .freqs      Hz, one or more (60, or [60 120 180] for harmonics). Default 60.
%   .window     s, length of the windows the sines are fitted in. Default 4
%               (windows overlap by half and are cross-faded).
%   .mask_ms    [start stop] ms around each stimulation left out of the fit.
%               Default [-10 300]. Used on CCEP data only (event types that
%               name a pair of contacts, see ieeglab_detect_mode).
%   .stim_mask_ms  [start stop] ms around each stimulation left out of the fit
%               of the contacts stimulated (the channels that use them).
%               Default [-10 2000]: on a stimulated contact the artifact's
%               recovery carries tens of uV at the line harmonics for 1 to 2 s
%               (tutorial sEEG, RA4 at 180 Hz: 64 uV 0.3-1 s after the pulse,
%               5 uV 1-2 s, 1.3 uV later); fitted with it, the sine put 13 to
%               21 uV of 180 Hz back into the rest of the recording. Where the
%               pulses are too close for it (less than 0.5 s left in a
%               window), mask_ms is used for them too.
%   .search_hz  the mains frequency is searched within +-search_hz of each
%               value (it can run a few hundredths of a Hz off). Default 0.1.
%   .verbose    logical. Default true.
%
% info: freqs (nominal), freq_found (Hz, per frequency, median over the
%   segments, weighted by length), amplitude (uV per channel, median over
%   windows, first frequency), masked (fraction of samples left out of the fit),
%   skipped (fraction of samples in no window with enough data away from the
%   pulses to fit; their line noise stays).
%
% Why not a notch filter, or ZapLine
% ---------------------------------
% A notch removes the narrow band around 60 Hz of everything, including each
% stimulation response and the brain's own 60 Hz activity, and a 2 Hz notch is a
% 3.3 s filter at 2048 Hz: what it removes is spread up to 1.6 s before and
% after the pulse. On the tutorials (artifact blanked first, each contact's own
% stimulation epochs left out, as every analysis does), it changes the trial
% average 10 to 50 ms after the pulse by 0.76 uV RMS (sEEG) and 0.94 uV (ECoG),
% the fit by 0.32 and 0.44 uV; in the stimulated contacts' own epochs, where the
% artifact remains, by 32 uV. On synthetic CCEP data with a known clean signal,
% the error left 5 to 50 ms after the pulse is 0.48 uV with the notch and
% 0.12 uV with the fit (stimulated contacts: 10.9 and 0.20 uV).
% ZapLine (de Cheveigne, 2020) and ZapLine-plus (Klug & Kloosterman, 2022)
% remove the spatial components that carry the line noise, estimated over the
% whole recording, pulses included. On CCEP data the pulses take part in those
% components: they changed the tutorials' trial average 10 to 50 ms after the
% pulse by 9 to 19 uV (ZapLine) and 9 to 400 uV (ZapLine-plus), left 18 to
% 21 uV of error after the pulse on synthetic data with line noise from two
% shared sources, and on the tutorial sEEG left most of the 60 Hz (36 and 14
% times the neighbouring frequencies, median over contacts; the fit 1.3).
%
% Here, in each window, a sine at each line frequency (its amplitude and phase)
% is fitted by least squares to the samples away from the pulses, then subtracted
% from all samples: FieldTrip's DFT filter (Oostenveld et al., 2011) in
% cross-faded windows, as CleanLine (Mullen, 2012) after Mitra and Pesaran
% (1999), so that slow changes of the line noise are followed (Leske & Dalal,
% 2019), with the fit restricted to the samples away from the pulses
% (least squares on unevenly sampled data; Lomb, 1976). The fit never sees a pulse or its response, so nothing is
% spread around them: on the full ds004696 sub-02 recording, 0.22 uV RMS is taken
% out of the trial average 900 to 500 ms before the pulse and 0.21 uV 10 to 50 ms
% after it (notch: 0.53 and 1.93 uV); what is left at 60 Hz equals the activity at
% the neighbouring frequencies.
%
% The frequency is found once per continuous stretch, from all of it (demodulated
% to 0.1 s blocks, so weak line noise still stands out), not window by window: in
% 4 s of weak line noise the best frequency is a noise peak, and fitting there
% removes brain activity and leaves line noise (tutorial ECoG: 1.26 uV left at
% 50 Hz, against 1.11 uV with the frequency found from the whole stretch and
% 1.06 uV at the neighbouring frequencies). Epochs and stretches shorter than
% 20 s use the frequency given.
%
% Assumes the line noise is steady over a few seconds (amplitude and phase), as
% mains noise usually is; a burst shorter than the window is only partly removed.
% Data split by 'boundary' events are fitted in segments, epoched data one epoch
% at a time (one window when an epoch is shorter than 1.5 windows). Samples that
% are NaN are left out of the fit and stay NaN.
%
% References
%   de Cheveigne, A. (2020). ZapLine: A simple and effective method to remove
%     power line artifacts. NeuroImage, 207, 116356.
%   Klug, M., & Kloosterman, N. A. (2022). Zapline-plus: A Zapline extension for
%     automatic and adaptive removal of frequency-specific noise artifacts in
%     M/EEG. Human Brain Mapping. https://doi.org/10.1002/hbm.25832
%   Leske, S., & Dalal, S. S. (2019). Reducing power line noise in EEG and MEG
%     data via spectrum interpolation. NeuroImage, 189, 763-776.
%   Lomb, N. R. (1976). Least-squares frequency analysis of unequally spaced
%     data. Astrophysics and Space Science, 39(2), 447-462.
%   Mitra, P. P., & Pesaran, B. (1999). Analysis of dynamic brain imaging data.
%     Biophysical Journal, 76(2), 691-708.
%   Mullen, T. (2012). CleanLine EEGLAB plugin. NITRC.
%   Oostenveld, R., Fries, P., Maris, E., & Schoffelen, J.-M. (2011). FieldTrip:
%     Open source software for advanced analysis of MEG, EEG, and invasive
%     electrophysiological data. Computational Intelligence and Neuroscience,
%     2011, 156869.
%
% Cedric Cannard, iEEGLAB, 2026

if nargin < 2 || isempty(opt), opt = struct(); end
def = struct('freqs', 60, 'window', 4, 'mask_ms', [-10 300], 'stim_mask_ms', [-10 2000], 'search_hz', 0.1, 'verbose', true);
fn = fieldnames(def);
for i = 1:numel(fn)
    if ~isfield(opt, fn{i}) || isempty(opt.(fn{i})), opt.(fn{i}) = def.(fn{i}); end
end
fs = double(EEG.srate);
f = double(opt.freqs(:))';
high = f >= fs/2 - 1 | f <= 0 | ~isfinite(f);
if any(high)
    warning('ieeglab_line_noise:aboveNyquist', ...
        'Line frequencies %s Hz skipped: not between 0 and Nyquist (%g Hz) minus 1 Hz.', mat2str(f(high)), fs/2);
    f = f(~high);
end
info = struct('freqs', f, 'freq_found', f, 'amplitude', NaN(EEG.nbchan, 1), 'masked', 0, 'skipped', 0);
if isempty(f) || isempty(EEG.data), return; end

[C, T, N] = size(EEG.data);
if C == 0, return; end
L = T * N;                                                % samples, epochs one after the other

% Segments fitted on their own: epochs, or the stretches between boundary events
if N > 1
    segs = [(0:N-1)' * T + 1, (1:N)' * T];
else
    b = [];
    if isfield(EEG, 'event') && ~isempty(EEG.event)
        isB = strcmpi(cellfun(@(x) char(string(x)), {EEG.event.type}, 'UniformOutput', false), 'boundary');
        b = sort(floor([EEG.event(isB).latency]));
        b = b(b >= 1 & b < L);
    end
    segs = [[1; b(:) + 1], [b(:); L]];
end

% Samples left out of the fit: around each stimulation of CCEP data, for every
% contact (keep), and longer for the contacts stimulated (stimIv: per channel,
% [first last] sample intervals)
keep = true(1, L);
stimIv = repmat({zeros(0, 2)}, C, 1);
if strcmp(ieeglab_detect_mode(EEG), 'ccep')
    ev = EEG.event(~strcmpi(cellfun(@(x) char(string(x)), {EEG.event.type}, 'UniformOutput', false), 'boundary'));
    w = round(double(opt.mask_ms) / 1000 * fs);
    wStim = round(double(opt.stim_mask_ms) / 1000 * fs);
    [~, si] = ieeglab_canonical_site(cellfun(@(x) char(string(x)), {ev.type}, 'UniformOutput', false), EEG.chanlocs);
    for e = 1:numel(ev)
        l = round(double(ev(e).latency));
        if ~isscalar(l) || ~isfinite(l), continue; end           % an event without a latency masks nothing
        keep(max(1, l + w(1)):min(L, l + w(2))) = false;
        for ch = si{e}(:)'
            stimIv{ch}(end+1, :) = [max(1, l + wStim(1)), min(L, l + wStim(2))];
        end
    end
end
info.masked = 1 - mean(keep);

win = max(round(opt.window * fs), 8);
hop = round(win / 2);
hann = 0.5 - 0.5 * cos(2 * pi * (0:win-1) / win);             % periodic: halves overlap-add to 1
minGood = max(round(0.5 * fs), ceil(10 * fs / min(f)));

% Pass 1: per segment, the mains frequency, then the sines of every window, fitted
% on the data as they are. Pass 2 subtracts them, cross-faded. Only the
% coefficients are kept in between, so the memory needed does not grow with the
% recording.
W = struct('idx', {}, 'w', {}, 'f', {}, 'c', {}, 't0', {});
found = zeros(size(segs, 1), numel(f));
for s = 1:size(segs, 1)
    a = segs(s, 1); z = segs(s, 2);
    fk = f;
    if opt.search_hz > 0 && z - a + 1 >= 20 * fs
        fk = local_frequency(EEG, a, z, keep, f, fs, opt.search_hz);
    end
    found(s, :) = fk;
    one = z - a + 1 <= 1.5 * win;                             % short: one window, unweighted
    if one, starts = a; else, starts = a - hop:hop:z; end
    for k = 1:numel(starts)
        if one
            idx = a:z; wk = ones(1, numel(idx));
        else
            idx = max(a, starts(k)):min(z, starts(k) + win - 1);
            wk = hann(idx - starts(k) + 1);
        end
        X = double(EEG.data(:, idx));                         % EEG.data(:, idx) runs over epochs too
        g = keep(idx) & all(isfinite(X), 1);
        if nnz(g) < minGood, continue; end
        t = (idx - a) / fs;                                   % time from the segment start
        A = local_design(t, fk);
        c = A(g, :) \ X(:, g)';                               % 2*nfreq x C
        % the stimulated contacts of a pulse in this window: fitted on their own,
        % without stim_mask_ms around that pulse (its artifact tail)
        for ch = find(cellfun(@(v) any(v(:, 1) <= idx(end) & v(:, 2) >= idx(1)), stimIv))'
            gc = g;
            for r = find(stimIv{ch}(:, 1) <= idx(end) & stimIv{ch}(:, 2) >= idx(1))'
                gc(idx >= stimIv{ch}(r, 1) & idx <= stimIv{ch}(r, 2)) = false;
            end
            if nnz(gc) >= minGood, c(:, ch) = A(gc, :) \ X(ch, gc)'; end   % else: pulses too close, mask_ms only
        end
        W(end+1) = struct('idx', [idx(1) idx(end)], 'w', wk, 'f', fk, 'c', c, 't0', a); %#ok<AGROW>
    end
end
if isempty(W)
    warning('ieeglab_line_noise:nothingFitted', ...
        'No window had enough data away from the pulses (%.1f s needed) to fit the line noise; data unchanged.', minGood / fs);
    return
end

ws = zeros(1, L, 'single');                                   % cross-fade weights, summed
for k = 1:numel(W), ws(W(k).idx(1):W(k).idx(2)) = ws(W(k).idx(1):W(k).idx(2)) + W(k).w; end
info.skipped = mean(ws == 0);                                 % no fitted window: line noise stays
cls = class(EEG.data);
for k = 1:numel(W)
    idx = W(k).idx(1):W(k).idx(2);
    r = W(k).w ./ double(ws(idx)); r(ws(idx) == 0) = 0;
    EEG.data(:, idx) = EEG.data(:, idx) - cast((local_design((idx - W(k).t0) / fs, W(k).f) * W(k).c)' .* r, cls);
end

len = segs(:, 2) - segs(:, 1) + 1;
for q = 1:numel(f)                                            % median over segments, weighted by length
    [v, o] = sort(found(:, q)); cw = cumsum(len(o)) / sum(len);
    info.freq_found(q) = v(find(cw >= 0.5, 1));
end
amp = cell2mat(arrayfun(@(x) hypot(x.c(1, :), x.c(2, :))', W, 'UniformOutput', false));
info.amplitude = median(amp, 2);
EEG.ieeglab.line_noise = info;

if opt.verbose
    [mx, im] = max(info.amplitude);
    lab = sprintf('channel %d', im);
    if isfield(EEG, 'chanlocs') && numel(EEG.chanlocs) >= im, lab = EEG.chanlocs(im).labels; end
    msk = '';
    if info.masked > 0, msk = sprintf(', %.0f%% of samples around the pulses left out of the fit', 100 * info.masked); end
    fprintf(['[line noise] Sines at %s Hz (found %s Hz) fitted in %g s windows and subtracted (no filter)%s. ' ...
             'Amplitude at %g Hz: median %.1f uV per channel, largest %.1f uV (%s).\n'], ...
        mat2str(f), mat2str(info.freq_found, 6), opt.window, msk, f(1), median(info.amplitude), mx, lab);
    if info.skipped > 0.01
        fprintf('[line noise] %.0f%% of the samples are in no window with enough data away from the pulses: their line noise stays.\n', ...
            100 * info.skipped);
    end
end
end

function fk = local_frequency(EEG, a, z, keep, f, fs, search)
% The mains frequency near each f, from the whole segment: the data are shifted
% by -f and summed in 0.1 s blocks (samples left out of the fit set to 0), so an
% offset d shows as a slow rotation exp(2i*pi*d*t) of the block sums; d is the
% one with the most power over all channels, on a 0.001 Hz grid. 30 s at a time.
blk = max(1, round(0.1 * fs));
nb = floor((z - a + 1) / blk);
d = -search:0.001:search;
fk = f;
for q = 1:numel(f)
    Z = complex(zeros(size(EEG.data, 1), nb));
    chunk = blk * max(1, floor(30 * fs / blk));
    for c0 = 0:chunk:nb * blk - 1
        idx = a + c0 + (0:min(chunk, nb * blk - c0) - 1);
        X = double(EEG.data(:, idx));
        X(:, ~keep(idx)) = 0; X(~isfinite(X)) = 0;
        Y = X .* exp(-2i * pi * f(q) * (idx - a) / fs);
        Z(:, c0 / blk + (1:numel(idx) / blk)) = squeeze(sum(reshape(Y, size(Y, 1), blk, []), 2));
    end
    tb = ((0:nb-1) * blk + (blk - 1) / 2) / fs;
    P = zeros(1, numel(d));
    for j = 1:numel(d)                                        % one offset at a time: memory stays small
        P(j) = sum(abs(Z * exp(-2i * pi * d(j) * tb')).^2);
    end
    [~, m] = max(P);
    fk(q) = f(q) + d(m);
end
end

function A = local_design(t, f)
A = zeros(numel(t), 2 * numel(f));
for q = 1:numel(f)
    A(:, 2*q-1) = cos(2 * pi * f(q) * t);
    A(:, 2*q)   = sin(2 * pi * f(q) * t);
end
end
