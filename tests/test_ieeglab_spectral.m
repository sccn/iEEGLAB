function tests = test_ieeglab_spectral
% Frequency-domain analysis (ieeglab_spectral, pop_ieeglab_spectral). Headless.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(mfilename('fullpath')));
tc.assumeNotEmpty(which('eeglab'), 'EEGLAB is not on the MATLAB path.');
if isempty(which('pop_loadset')), evalc('eeglab nogui'); end
addpath(root); addpath(fullfile(root, 'functions'));
tc.TestData.root = root;
end

function test_spectrum_units(tc)
% Welch PSD in uV^2/Hz: a sinusoid of amplitude A has band power A^2/2, white
% noise of SD s integrates to s^2 (Parseval).
fs = 1000; T = 20 * fs; t = (0:T-1) / fs;
rng(1);
X = [3 * sin(2 * pi * 40 * t) + 0.1 * randn(1, T); 2 * randn(1, T); randn(1, T)];
EEG = local_eeg(X, fs, 0);
S = ieeglab_spectral(EEG, 'type', 'spectrum', 'freqs', [0 500], 'band', [35 45], 'verbose', false);
tc.verifyEqual(S.band_value(1), 4.5, 'RelTol', 0.02, 'Band power of a 3 uV sinusoid must be 4.5 uV^2.');
df = S.freqs(2) - S.freqs(1);
tc.verifyEqual(sum(S.power(2, :)) * df, 4, 'RelTol', 0.05, 'White noise of SD 2 must integrate to 4 uV^2.');
[~, k] = max(S.power(1, :));
tc.verifyEqual(S.freqs(k), 40);
end

function test_spectrum_skips_discontinuities(tc)
% No Welch segment may span a 'boundary' event of continuous data.
fs = 500; X = randn(2, 10 * fs);
EEG = local_eeg(X, fs, 0);
S0 = ieeglab_spectral(EEG, 'type', 'spectrum', 'verbose', false);
EEG.event = struct('type', 'boundary', 'latency', 5 * fs + 0.5, 'duration', 0);
S1 = ieeglab_spectral(EEG, 'type', 'spectrum', 'verbose', false);
tc.verifyLessThan(S1.n_segments, S0.n_segments);
end

function test_timefreq_matches_direct_sum(tc)
% The FFT convolution must equal the direct sum of the Hann-tapered complex
% sinusoid over the window, for every output time, at two frequencies, with a
% fixed window and with cycles.
fs = 500; T = 600; N = 3; rng(2);
X = randn(2, T, N) * 10;
EEG = local_eeg(X, fs, -0.4);
for opts = {{'window', 0.2}, {'cycles', 4}}
    S = ieeglab_spectral(EEG, 'freqs', [20 60], 'norm', 'none', 'timestep', 30, 'verbose', false, opts{1}{:});
    for k = 1:2
        hL = (round(S.window_s(k) * fs) - 1) / 2;
        w = 0.5 - 0.5 * cos(2 * pi * (1:2*hL+1)' / (2 * hL + 2));
        ker = w .* exp(-1i * 2 * pi * S.freqs(k) * (-hL:hL)' / fs);
        for j = find(isfinite(S.power(1, k, :)))'
            s = find(abs(double(EEG.times) - S.times(j)) < 1e-9);
            pw = zeros(2, 1);
            for n = 1:N
                x = double(EEG.data(:, :, n)); x = x - mean(x, 2);   % the stored (single) data
                cf = x(:, s-hL:s+hL) * ker * 2 / sum(w);
                pw = pw + abs(cf).^2 / N;
            end
            tc.verifyEqual(squeeze(S.power(:, k, j)), pw, 'RelTol', 1e-9);
        end
    end
end
end

function test_timefreq_amplitude_and_baseline(tc)
% A 2 uV, 80 Hz sinusoid after time 0 only: power A^2 = 4 uV^2 at 80 Hz once
% the window is past 0, nothing before; in dB the post-stimulus 80 Hz power is
% far above the baseline and the band value is positive.
fs = 1000; T = 2000; N = 10; rng(3);
t = (0:T-1) / fs - 1;
X = 0.05 * randn(1, T, N) + reshape(repmat(2 * sin(2 * pi * 80 * t) .* (t >= 0), 1, N), 1, T, N);
EEG = local_eeg(X, fs, -1);
S = ieeglab_spectral(EEG, 'freqs', 60:10:100, 'norm', 'none', 'verbose', false);
k = find(S.freqs == 80);
tc.verifyEqual(S.power(1, k, find(S.times >= 500, 1)), 4, 'RelTol', 0.03);
tc.verifyLessThan(S.power(1, k, find(S.times >= -500, 1)), 0.01);
Sd = ieeglab_spectral(EEG, 'freqs', 60:10:100, 'baseline', [-800 0], 'band', [75 85], 'band_window', [200 800], 'verbose', false);
tc.verifyGreaterThan(Sd.band_value(1), 20, 'An 80 Hz burst must be > 20 dB above a near-silent baseline.');
tc.verifyTrue(all(isnan(Sd.tf(1, :, Sd.times > 900)), 'all'), 'Windows past the epoch end must be NaN.');
end

function test_baseline_without_whole_window_is_nan(tc)
% With 5 cycles at 5 Hz (1 s window) no window fits in a 0.5 s baseline.
fs = 250; X = randn(1, 500, 4);
EEG = local_eeg(X, fs, -0.5);
tc.verifyWarning(@() ieeglab_spectral(EEG, 'freqs', [5 40], 'cycles', 5), 'ieeglab_spectral:noBaseline');
[~, S] = evalc('ieeglab_spectral(EEG, ''freqs'', [5 40], ''cycles'', 5, ''verbose'', false)');
tc.verifyTrue(all(isnan(S.tf(:, 1, :)), 'all'));
end

function test_ccep_own_stimulation_left_out(tc)
% Tutorial sEEG (CCEP): each contact averages the epochs in which it was not
% stimulated; history line replays.
d = fullfile(tc.TestData.root, 'tutorial', 'seeg');
EEG = pop_loadset('filename', 'sub-02_ses-ieeg01_task-ccep_run-01_ieeg.set', 'filepath', d);
[~, EEG] = evalc('ieeglab_load(EEG, struct())');
EEG0 = EEG;                                                    % continuous
pre = struct('apply_highpass',true,'highpass',0.5,'apply_notch',false,'apply_lowpass',false, ...
    'apply_ds',false,'apply_epoch',true,'epoch_window',[-500 1000],'remove_bad_channels',true, ...
    'apply_baseline',true,'baseline_period',[-500 -50],'apply_car',false,'plot',false,'verbose',false);
[~, E] = evalc('ieeglab_preprocess(EEG, pre)');
[~, stimIdx] = ieeglab_epoch_sites(E);
nOwn = zeros(E.nbchan, 1);
for n = 1:E.trials, nOwn(stimIdx{n}) = nOwn(stimIdx{n}) + 1; end
[~, B, com] = evalc('pop_ieeglab_spectral(E, ''type'', ''timefreq'', ''freqs'', 10:10:50, ''plot'', false)');
S = B.ieeglab.spectral;
tc.verifyEqual(S.n_trials, E.trials - nOwn);
tc.verifyEqual(size(S.tf), [E.nbchan 5 numel(S.times)]);
EEG = E;
[~, B2] = evalc(com(7:end-1));
tc.verifyEqual(B2.ieeglab.spectral.tf, S.tf);
% windows that contain the pulse (0.2 s windows, artifact -5..15 ms) are NaN; others are not
tc.verifyTrue(all(isnan(S.tf(:, :, abs(S.times - 5) < 105)), 'all'));
tc.verifyTrue(any(isfinite(S.tf(:, :, S.times > 300)), 'all'));
% continuous CCEP spectrum: no segment overlaps -50..1000 ms around a pulse
[~, Sc] = evalc('ieeglab_spectral(EEG0, ''type'', ''spectrum'')');
[~, Sa] = evalc('ieeglab_spectral(EEG0, ''type'', ''spectrum'', ''stim_exclude'', [])');
tc.verifyLessThan(Sc.n_segments, Sa.n_segments);
tc.verifyError(@() ieeglab_spectral(EEG, 'type', 'timefreq', 'condition', 'XX1-XX2'), 'ieeglab_spectral:noCondition');
end

function test_spectrum_matches_pwelch(tc)
% Same estimate as MATLAB's pwelch (periodic Hann, 50% overlap, one-sided),
% continuous and epoched (mean over epochs). iEEGLAB also removes each
% segment's mean, which a Hann window spreads over the first two bins only.
tc.assumeNotEmpty(which('pwelch'), 'Signal Processing Toolbox not installed.');
fs = 512; L = fs; rng(5);
for N = [1 4]
    T = 3740 * (N == 1) + 900 * (N > 1);
    t = (0:T-1) / fs;
    X = 30 + randn(3, T, N) .* [1; 2; 3] + 4 * sin(2 * pi * 37.5 * t);
    EEG = local_eeg(X, fs, -0.3 * (N > 1));
    S = ieeglab_spectral(EEG, 'type', 'spectrum', 'freqs', [0 fs / 2], 'verbose', false);
    ref = zeros(L / 2 + 1, 3);
    for n = 1:N
        [p, f] = pwelch(double(EEG.data(:, :, n))', hann(L, 'periodic'), L / 2, L, fs);
        ref = ref + p / N;
    end
    tc.verifyEqual(S.freqs(:), f);
    k = f >= 2 * fs / L;
    tc.verifyEqual(S.power(:, k)', ref(k, :), 'RelTol', 1e-9, sprintf('%d epoch(s): differs from pwelch.', N));
end
end

function test_timefreq_matches_fieldtrip(tc)
% Same time-frequency power as FieldTrip's mtmconvol with a Hann taper and the
% same window, up to FieldTrip's constant scaling (its taper has unit norm and
% the coefficients are multiplied by sqrt(2/n)): same times, same NaN edges.
tc.assumeNotEmpty(which('ft_freqanalysis'), 'FieldTrip is not on the path.');
fs = 500; T = 600; N = 5; rng(6);
X = randn(3, T, N) * 5;
EEG = local_eeg(X, fs, -0.4);
S = ieeglab_spectral(EEG, 'freqs', [10 25 50 100], 'norm', 'none', 'timestep', 20, 'verbose', false);
data = struct('label', {{'LA1'; 'LA2'; 'LA3'}}, 'fsample', fs, ...
    'trial', {squeeze(num2cell(double(EEG.data), [1 2]))'}, 'time', {repmat({double(EEG.times) / 1000}, 1, N)});
cfg = struct('method', 'mtmconvol', 'taper', 'hanning', 'output', 'pow', 'foi', S.freqs, ...
    't_ftimwin', S.window_s, 'toi', S.times / 1000, 'keeptrials', 'no', 'feedback', 'no');
[~, F] = evalc('ft_freqanalysis(cfg, data)');
tc.verifyEqual(F.freq, S.freqs, 'AbsTol', 1e-9);
tc.verifyEqual(isnan(F.powspctrm), isnan(S.power), 'Windows outside the epoch must be NaN in both.');
n = round(S.window_s(1) * fs);
w = hanning(n);
K = (2 / sum(w))^2 / (2 / n / sum(w.^2));
r = S.power ./ F.powspctrm;
tc.verifyEqual(r(isfinite(r)), K * ones(nnz(isfinite(r)), 1), 'RelTol', 1e-8);
end

function test_timefreq_cycles_amplitude(tc)
% With a number of cycles, each frequency has its own window and a sinusoid
% of amplitude A still has power A^2.
fs = 1000; T = 3000; N = 4;
t = (0:T-1) / fs - 1.5;
f0 = [10 40 120]; A = [1 2 3];
X = reshape(repmat(A(:) .* sin(2 * pi * f0(:) * t), 1, 1, N), 3, T, N);
EEG = local_eeg(X, fs, -1.5);
S = ieeglab_spectral(EEG, 'freqs', f0, 'cycles', 6, 'norm', 'none', 'verbose', false);
tc.verifyEqual(S.window_s, (2 * round(6 ./ f0 * fs / 2) + 1) / fs);
j = find(S.times >= 0, 1);
for k = 1:3
    tc.verifyEqual(S.power(k, k, j), A(k)^2, 'RelTol', 0.01, sprintf('%g Hz', f0(k)));
end
end

function test_blocks_and_single_trials(tc)
% Splitting the epochs (or segments) into blocks changes nothing, and the
% single-trial power averages to the power over the epochs each contact uses.
fs = 500; rng(8);
EEG = local_eeg(randn(3, 400, 7) * 10, fs, -0.3);
S1 = ieeglab_spectral(EEG, 'freqs', [20 40 80], 'keeptrials', true, 'verbose', false);
S2 = ieeglab_spectral(EEG, 'freqs', [20 40 80], 'max_mb', 1e-4, 'verbose', false);
tc.verifyEqual(S2.power, S1.power, 'RelTol', 1e-12);
m = mean(double(S1.trials), 4);
tc.verifyEqual(m(isfinite(m)), S1.power(isfinite(S1.power)), 'RelTol', 1e-6);
tc.verifyEqual(isnan(m), isnan(S1.power));
P1 = ieeglab_spectral(EEG, 'type', 'spectrum', 'window', 0.2, 'verbose', false);
P2 = ieeglab_spectral(EEG, 'type', 'spectrum', 'window', 0.2, 'max_mb', 1e-4, 'verbose', false);
tc.verifyEqual(P2.power, P1.power, 'RelTol', 1e-12);
end

function test_condition_and_nyquist(tc)
% 'condition' gives the result of the dataset holding only those epochs;
% frequencies at or above Nyquist are dropped with a warning.
fs = 500; rng(9);
X = randn(2, 500, 6) * 10;
EEG = local_eeg(X, fs, -0.5);
for n = 2:2:6
    EEG.event(n).type = 'noise'; EEG.epoch(n).eventtype = {'noise'};
end
S = ieeglab_spectral(EEG, 'freqs', [20 60], 'condition', 'noise', 'verbose', false);
S0 = ieeglab_spectral(local_eeg(X(:, :, 2:2:6), fs, -0.5), 'freqs', [20 60], 'verbose', false);
tc.verifyEqual(S.power, S0.power, 'RelTol', 1e-12);
tc.verifyEqual(S.n_trials, [3; 3]);
tc.verifyWarning(@() ieeglab_spectral(EEG, 'freqs', [100 250 300]), 'ieeglab_spectral:aboveNyquist');
[~, S] = evalc('ieeglab_spectral(EEG, ''freqs'', [100 250 300], ''verbose'', false)');
tc.verifyEqual(S.freqs, 100);
end

function test_pulses_on_removed_contacts_still_excluded(tc)
% Continuous CCEP data whose stimulated contacts were removed (bad): the
% pulses still count, so no segment overlaps them, and the dataset is still CCEP.
fs = 500; T = 120 * fs; rng(12);
EEG = local_eeg(randn(4, T), fs, 0);
for c = 1:4, EEG.chanlocs(c).labels = sprintf('RA%d', c + 2); end    % RA3..RA6; RA1, RA2 removed
lat = (5:4:115) * fs;
EEG.event = struct('type', repmat({'RA1-RA2'}, 1, numel(lat)), 'latency', num2cell(lat), 'duration', 0);
tc.verifyEqual(ieeglab_detect_mode(EEG), 'ccep');
S = ieeglab_spectral(EEG, 'type', 'spectrum', 'verbose', false);
L = fs; starts = 1:L/2:(T - L + 1);
w = round([-50 1000] / 1000 * fs);
for b = lat, starts = starts(starts + L - 1 < b + w(1) | starts > b + w(2)); end
tc.verifyEqual(S.n_segments, numel(starts));
E = EEG; E.event(1).type = 'face-house';             % ordinary conditions stay what they are
[~, info] = ieeglab_detect_mode(E);
tc.verifyFalse(ismember("face-house", info.pair_types));
end

function test_nan_samples_leave_out_their_windows(tc)
% A blanked artifact left as NaN: windows that touch it are NaN, every other
% window equals the result on the same data with the NaN samples filled.
fs = 500; T = 600; N = 4; rng(13);
X = randn(2, T, N) * 10;
EEG = local_eeg(X, fs, -0.4);
k0 = find(EEG.times >= -2 & EEG.times <= 10);
EEG.data(:, k0, :) = NaN;
EEG.data(2, :, 3) = NaN;                              % one contact, one whole epoch
F = EEG;
for c = 1:2
    for n = 1:N
        x = F.data(c, :, n); x(isnan(x)) = mean(x, 'omitnan'); F.data(c, :, n) = x;
    end
end
F.data(2, :, 3) = 0;
S = ieeglab_spectral(EEG, 'freqs', [20 60], 'norm', 'none', 'verbose', false);
R = ieeglab_spectral(F, 'freqs', [20 60], 'norm', 'none', 'verbose', false);
h = (round(S.window_s(1) * fs) - 1) / 2 / fs * 1000;   % half window, ms
hit = abs(S.times - 4) <= 6 + h;                        % windows reaching -2..10 ms
tc.verifyTrue(all(isnan(S.power(:, :, hit)), 'all'));
ok = isfinite(S.power(1, :, :));
tc.verifyTrue(any(ok(:)), 'Windows away from the NaN samples must be finite.');
tc.verifyEqual(S.power(1, ok), R.power(1, ok), 'RelTol', 1e-9);
% contact 2 averages 3 epochs (its third is all NaN): its power equals the
% average over the other three
G = F; G.data = F.data(:, :, [1 2 4]); G.trials = 3; G.epoch = F.epoch([1 2 4]); G.event = F.event([1 2 4]);
R2 = ieeglab_spectral(G, 'freqs', [20 60], 'norm', 'none', 'verbose', false);
ok2 = isfinite(S.power(2, :, :));
tc.verifyEqual(S.power(2, ok2), R2.power(2, ok2), 'RelTol', 1e-9);
% spectrum: the segments with a NaN sample are left out, per contact and epoch
tc.assumeNotEmpty(which('pwelch'), 'Signal Processing Toolbox not installed.');
Sp = ieeglab_spectral(EEG, 'type', 'spectrum', 'window', 0.2, 'freqs', [0 250], 'verbose', false);
L = 100; starts = 1:50:(T - L + 1);
clean = starts + L - 1 < k0(1) | starts > k0(end);
ref = zeros(L / 2 + 1, 1);
for n = 1:N
    seg = double(EEG.data(1, :, n));
    P = zeros(L / 2 + 1, 1);
    for s = starts(clean)
        [p, f] = pwelch(seg(s:s + L - 1)', hann(L, 'periodic'), 0, L, fs);
        P = P + p / nnz(clean);
    end
    ref = ref + P / N;
end
kk = f >= 2 * fs / L;
tc.verifyEqual(Sp.power(1, kk)', ref(kk), 'RelTol', 1e-9);
tc.verifyEqual(Sp.n_trials, [4; 3]);
end

function test_band_power_integrates_white_noise(tc)
% Band power of white noise of variance 1: the one-sided PSD (2/fs per Hz)
% times the band width, with no extra half bin at each edge.
fs = 1000; rng(14);
EEG = local_eeg(randn(1, 400 * fs), fs, 0);
S = ieeglab_spectral(EEG, 'type', 'spectrum', 'band', [8 12], 'verbose', false);
tc.verifyEqual(S.band_value, 4 * 2 / fs, 'RelTol', 0.05);
end

function test_epoched_ccep_spectrum_avoids_the_pulse(tc)
% Epoched CCEP data: segments stay outside -50..1000 ms around the pulse; the
% default window shrinks to fit the part before it; an explicit window that
% does not fit is an error.
d = fullfile(tc.TestData.root, 'tutorial', 'seeg');
EEG = pop_loadset('filename', 'sub-02_ses-ieeg01_task-ccep_run-01_ieeg.set', 'filepath', d);
[~, EEG] = evalc('ieeglab_load(EEG, struct())');
pre = struct('apply_epoch', true, 'epoch_window', [-1000 1000], 'remove_bad_channels', true, ...
    'plot', false, 'verbose', false);
[~, E] = evalc('ieeglab_preprocess(EEG, pre)');
[~, S] = evalc('ieeglab_spectral(E, ''type'', ''spectrum'')');
i0 = find(E.times >= 0, 1);
L = round(S.window_s * E.srate);
tc.verifyEqual(L, i0 - round(0.05 * E.srate) - 1, 'The longest stretch before -50 ms.');
tc.verifyEqual(S.n_segments, 1);
tc.verifyError(@() ieeglab_spectral(E, 'type', 'spectrum', 'window', 1, 'verbose', false), 'ieeglab_spectral:noSegments');
[~, Sall] = evalc('ieeglab_spectral(E, ''type'', ''spectrum'', ''stim_exclude'', [])');
tc.verifyGreaterThan(Sall.n_segments, S.n_segments);
end

function test_history_replays_exact_frequencies(tc)
% The history line of the frequency analysis gives back the same frequencies,
% even when they are not round numbers.
fs = 500; rng(15);
EEG = local_eeg(randn(2, 600, 3), fs, -0.4);
f = (3:12) * 7 / 3;
[~, ~, com] = evalc('pop_ieeglab_spectral(EEG, ''type'', ''timefreq'', ''freqs'', f, ''plot'', false)');
S0 = ieeglab_spectral(EEG, 'freqs', f, 'verbose', false);
[~, B] = evalc(com(7:end-1));
tc.verifyEqual(B.ieeglab.spectral.freqs, S0.freqs, 'AbsTol', 1e-12);
end

% ------------------------------------------------------------------ helpers

function EEG = local_eeg(X, fs, xmin)
% Synthetic dataset from data X (channels x time x epochs); continuous when
% X has one epoch and xmin is 0.
[C, T, N] = size(X);
EEG = eeg_emptyset();
EEG.srate = fs; EEG.nbchan = C; EEG.pnts = T; EEG.trials = N;
EEG.data = single(X);
EEG.xmin = xmin; EEG.xmax = xmin + (T - 1) / fs;
EEG.times = (xmin + (0:T-1) / fs) * 1000;
for c = C:-1:1
    cl(c) = struct('labels', sprintf('LA%d', c), 'type', 'SEEG', 'X', c, 'Y', 0, 'Z', 0, 'status', 'good'); %#ok<AGROW>
end
EEG.chanlocs = cl;
EEG.event = struct('type', {}, 'latency', {}, 'duration', {});
if N > 1
    i0 = find(EEG.times >= 0, 1);
    for n = 1:N
        EEG.event(n) = struct('type', 'tone', 'latency', (n - 1) * T + i0, 'duration', 0);
        EEG.epoch(n).event = n; EEG.epoch(n).eventtype = {'tone'}; EEG.epoch(n).eventlatency = {0};
    end
end
end
