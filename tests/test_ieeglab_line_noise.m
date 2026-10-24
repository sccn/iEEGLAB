function tests = test_ieeglab_line_noise
% test_ieeglab_line_noise - Line-noise removal by fitted sines (ieeglab_line_noise):
% it removes the line noise, takes nothing extra around the pulses (a notch
% does), and handles boundaries, epochs, NaN and frequencies above Nyquist.
% Synthetic data with a known clean signal, plus the tutorial sEEG.
%
% Cedric Cannard, iEEGLAB, 2026

tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(mfilename('fullpath')));
tc.TestData.root = root;
tc.assumeNotEmpty(which('eeglab'), 'EEGLAB is not on the MATLAB path.');
if isempty(which('pop_loadset')), evalc('eeglab nogui'); end
addpath(root); addpath(fullfile(root, 'functions'));
end

% ======================= helpers =======================

function [EEG, S, line, lat] = local_ccep(dur, fs)
% Six contacts, a pulse every 2.5 s with a residual artifact spike (2000 uV,
% 1 ms) and a large response; mains at 60.02 Hz, 20 uV, amplitude varying by 20%.
rs = RandStream('mt19937ar', 'Seed', 1);
C = 6; n = round(dur * fs); t = (0:n-1) / fs;
S = 5 * randn(rs, C, n);
lat = round((1:2.5:dur-1.5) * fs) + 1;
tau = (0:round(0.6 * fs)) / fs;
resp = -300 * exp(-((tau - 0.020) / 0.005).^2) + 100 * exp(-((tau - 0.150) / 0.05).^2);
resp(1:round(0.001 * fs)) = resp(1:round(0.001 * fs)) + 2000;
for l = lat
    S(:, l + (0:numel(tau)-1)) = S(:, l + (0:numel(tau)-1)) + (1 + (0:C-1)' / C) .* resp;
end
ph = 2 * pi * rand(rs, C, 1);
line = 20 * (1 + 0.2 * sin(2 * pi * t / 20)) .* sin(2 * pi * 60.02 * t + ph);
labels = compose('C%d', 1:C);
ev = struct('type', 'C1-C2', 'latency', num2cell(lat));
EEG = struct('data', single(S + line), 'srate', fs, 'nbchan', C, 'pnts', n, 'trials', 1, ...
    'chanlocs', struct('labels', labels), 'event', ev, 'xmin', 0);
end

function r = local_rms_at(E, lat, fs, win_ms)
% RMS of E over windows win_ms around each latency, single trials and trial average
k = round(win_ms(1) / 1000 * fs):round(win_ms(2) / 1000 * fs);
X = zeros(size(E, 1), numel(k), numel(lat));
for i = 1:numel(lat), X(:, :, i) = E(:, lat(i) + k); end
r = [sqrt(mean(X.^2, 'all')), sqrt(mean(mean(X, 3).^2, 'all'))];
end

% ======================= tests =======================

function test_removes_line_noise_and_nothing_around_the_pulses(tc)
fs = 2048;
[EEG, S, line, lat] = local_ccep(60, fs);
Y = ieeglab_line_noise(EEG, struct('freqs', 60, 'verbose', false));
err = double(Y.data) - S;                       % what is left of the line noise, and any damage
lineRMS = sqrt(mean(line.^2, 'all'));
far  = local_rms_at(err, lat, fs, [-900 -500]);
near = local_rms_at(err, lat, fs, [5 50]);
% Measured 23 Oct: 0.32 and 0.17 uV per trial, 0.03 uV in the average around the pulses (line
% RMS 14.3 uV). Without the pulses left out of the fit: 0.96, 0.89 and 0.84 uV.
tc.verifyLessThan(far(1), 0.05 * lineRMS, 'Line noise left 900-500 ms before the pulses.');
tc.verifyLessThan(near(1), 0.05 * lineRMS, 'Line noise left, or response changed, around the pulses.');
tc.verifyLessThan(near(2), 0.3, 'The trial average around the pulses changed.');
tc.verifyLessThan(near(2), 1.5 * far(2) + 0.2, 'More taken out around the pulses than elsewhere.');
tc.verifyEqual(Y.ieeglab.line_noise.freq_found, 60.02, 'AbsTol', 0.0015);
tc.verifyEqual(Y.ieeglab.line_noise.amplitude, 20 * ones(6, 1), 'RelTol', 0.25);
% The reason it is the default: a notch takes the response's own 60 Hz out and spreads it
if ~isempty(which('pop_eegfiltnew'))
    [~, N] = evalc('pop_eegfiltnew(EEG, ''locutoff'', 59, ''hicutoff'', 61, ''usefftfilt'', 1, ''revfilt'', 1)');
    nearN = local_rms_at(double(N.data) - S, lat, fs, [5 50]);
    tc.verifyGreaterThan(nearN(2), 5 * near(2), 'The notch should change the average around the pulses far more.');
end
end

function test_boundaries_are_fitted_apart(tc)
% A boundary where the mains phase jumps by half a cycle (data removed there)
fs = 1000; n = 40 * fs; t = (0:n-1) / fs;
rs = RandStream('mt19937ar', 'Seed', 2);
S = 3 * randn(rs, 4, n);
line = 15 * sin(2 * pi * 50 * t + pi * (t >= 20));
EEG = struct('data', S + line, 'srate', fs, 'nbchan', 4, 'pnts', n, 'trials', 1, ...
    'chanlocs', struct('labels', {'A', 'B', 'C', 'D'}), ...
    'event', struct('type', 'boundary', 'latency', 20 * fs + 0.5), 'xmin', 0);
Y = ieeglab_line_noise(EEG, struct('freqs', 50, 'verbose', false));
err = Y.data - S;
blk = reshape(sqrt(mean(reshape(err, 4, fs / 2, []).^2, [1 2])), 1, []);   % RMS per half second
tc.verifyLessThan(max(blk), 0.15 * 15 / sqrt(2), 'Line noise left near the boundary or elsewhere.');
tc.verifyEqual(Y.ieeglab.line_noise.masked, 0, 'Nothing is masked in data without stimulation.');
end

function test_epoched_data(tc)
fs = 2048;
[EEG, S, ~, lat] = local_ccep(80, fs);
k = -fs:fs-1;                                    % epochs of 2 s around each pulse
X = zeros(EEG.nbchan, numel(k), numel(lat), 'single'); R = zeros(size(X));
for i = 1:numel(lat), X(:, :, i) = EEG.data(:, lat(i) + k); R(:, :, i) = S(:, lat(i) + k); end
E = EEG; E.data = X; E.pnts = numel(k); E.trials = numel(lat); E.xmin = -1;
E.event = struct('type', 'C1-C2', 'latency', num2cell((0:numel(lat)-1) * numel(k) + fs + 1), ...
    'epoch', num2cell(1:numel(lat)));
Y = ieeglab_line_noise(E, struct('freqs', 60, 'verbose', false));
tc.verifySize(Y.data, size(X));
err = double(Y.data) - R;
far = sqrt(mean(err(:, 1:round(0.5 * fs), :).^2, 'all'));
near = sqrt(mean(mean(err(:, fs + (round(0.005 * fs):round(0.05 * fs)), :), 3).^2, 'all'));
tc.verifyLessThan(far, 0.15 * 20 / sqrt(2), 'Line noise left in the epochs.');
tc.verifyLessThan(near, 2, 'The trial average around the pulse changed.');
end

function test_nan_and_frequencies_above_nyquist(tc)
fs = 2048;
[EEG, S] = local_ccep(30, fs);
EEG.data(3, 10000:12000) = NaN;
Y = tc.verifyWarning(@() ieeglab_line_noise(EEG, struct('freqs', [60 1200], 'verbose', false)), ...
    'ieeglab_line_noise:aboveNyquist');
tc.verifyTrue(all(isnan(Y.data(3, 10000:12000))), 'NaN samples stay NaN.');
tc.verifyEqual(nnz(isnan(Y.data)), 2001, 'No other sample becomes NaN.');
ok = [1:9000, 13000:size(S, 2)];
tc.verifyLessThan(sqrt(mean((double(Y.data(:, ok)) - S(:, ok)).^2, 'all')), 2, 'Line noise left.');
EEG.srate = 100;                                % 60 Hz above Nyquist: nothing done
Z = tc.verifyWarning(@() ieeglab_line_noise(EEG, struct('freqs', 60, 'verbose', false)), ...
    'ieeglab_line_noise:aboveNyquist');
tc.verifyTrue(isequaln(Z.data, EEG.data), 'Data changed although no frequency is usable.');
end

function test_preprocess_default_on_the_tutorial(tc)
% ieeglab_preprocess removes the line noise with the sine fit unless told
% 'line_method', 'notch'; on the tutorial sEEG the fit takes about the same out
% of the trial average around the pulse as before it, the notch several times more.
d = fullfile(tc.TestData.root, 'tutorial', 'seeg');
EEG = pop_loadset('filename', 'sub-02_ses-ieeg01_task-ccep_run-01_ieeg.set', 'filepath', d);
[~, EEG] = evalc('ieeglab_load(EEG, struct(''event_field'', ''electrical_stimulation_site''))');
o = struct('apply_blank', true, 'apply_epoch', true, 'epoch_window', [-1000 1000], ...
    'remove_bad_channels', true, 'plot', false, 'verbose', false);
[~, E0] = evalc('ieeglab_preprocess(EEG, setfield(o, ''apply_notch'', false))'); %#ok<SFLD>
o.apply_notch = true; o.notch = 60;
[~, Ef] = evalc('ieeglab_preprocess(EEG, o)');
o.line_method = 'notch';
[~, En] = evalc('ieeglab_preprocess(EEG, o)');
tc.verifyTrue(isfield(Ef.ieeglab, 'line_noise'), 'The default is the sine fit.');
tc.verifyFalse(isfield(En.ieeglab, 'line_noise'), 'line_method ''notch'' uses the notch.');
t = E0.times;
avgRMS = @(E, w) sqrt(mean(mean(E0.data(:, t >= w(1) & t <= w(2), :) - E.data(:, t >= w(1) & t <= w(2), :), 3).^2, 'all'));
fit = [avgRMS(Ef, [-900 -500]), avgRMS(Ef, [10 50])];
notch = [avgRMS(En, [-900 -500]), avgRMS(En, [10 50])];
tc.verifyLessThan(fit(2), 2 * fit(1), 'The fit takes more out around the pulse than before it.');
tc.verifyGreaterThan(notch(2), 3 * fit(2), 'Expected the notch to change the average around the pulse more.');
end

function test_event_without_latency_masks_nothing(tc)
fs = 2048;
[EEG, S] = local_ccep(30, fs);
EEG.event(end+1) = struct('type', 'C1-C2', 'latency', NaN);
Y = ieeglab_line_noise(EEG, struct('freqs', 60, 'verbose', false));
tc.verifyLessThan(Y.ieeglab.line_noise.masked, 0.2, 'A NaN latency must not mask the recording.');
tc.verifyLessThan(sqrt(mean((double(Y.data) - S).^2, 'all')), 1, 'Line noise left.');
end

function test_mains_frequency_off_by_a_twentieth_of_a_hz(tc)
% Mains at 49.95 Hz, given 50: found from the whole recording, removed as well
fs = 1000; n = 120 * fs; t = (0:n-1) / fs;
rs = RandStream('mt19937ar', 'Seed', 3);
S = 3 * randn(rs, 4, n);
line = 4 * sin(2 * pi * 49.95 * t + (1:4)');
EEG = struct('data', S + line, 'srate', fs, 'nbchan', 4, 'pnts', n, 'trials', 1, ...
    'chanlocs', struct('labels', {'A', 'B', 'C', 'D'}), 'event', [], 'xmin', 0);
Y = ieeglab_line_noise(EEG, struct('freqs', 50, 'verbose', false));
tc.verifyEqual(Y.ieeglab.line_noise.freq_found, 49.95, 'AbsTol', 0.0015);
ry = sqrt(mean((Y.data - S).^2, 'all'));
tc.verifyLessThan(ry, 0.05 * 4 / sqrt(2), 'Line noise left.');
Z = ieeglab_line_noise(EEG, struct('freqs', 50, 'search_hz', 0, 'verbose', false));
tc.verifyGreaterThan(sqrt(mean((Z.data - S).^2, 'all')), 2 * ry, ...
    'At the nominal frequency more should be left (11% of the line noise when measured).');
end
