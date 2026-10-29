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
