function tests = test_ieeglab_band_stats
% Band-power permutation statistics (ieeglab_band_stats) and the FieldTrip
% import (ieeglab_import_fieldtrip). Headless.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(mfilename('fullpath')));
tc.assumeNotEmpty(which('eeglab'), 'EEGLAB is not on the MATLAB path.');
if isempty(which('pop_loadset')), evalc('eeglab nogui'); end
addpath(root); addpath(fullfile(root, 'functions'));
end

function test_exact_matches_brute_force(tc)
% With all 2^K sign patterns, t and p must equal a direct computation that
% builds each epoch's statistic from S.trials, flips the signs and recomputes
% mean and SD for every pattern (the variance shortcut must not change
% anything), for both tails and one-sided; contact 2 averages 6 of the 9 epochs.
rng(4);
K = 9; C = 4;
S = local_S(C, K);
S.trials = S.trials .* reshape(1 + 0.3 * (1:C)' .* (S.times >= 0), C, 1, [], 1);   % effects of several sizes
S.use(2, [1 3 5]) = false;
D = local_d_ref(S, [0 400], [-500 0]);
for tail = {'both', 'positive', 'negative'}
    [T, R] = ieeglab_band_stats(S, 'window', [0 400], 'baseline', [-500 0], 'tail', tail{1}, ...
        'min_trials', 2, 'correct', 'none', 'verbose', false);
    tc.verifyTrue(all(R.exact));
    tc.verifyEqual(R.d(S.use), D(S.use), 'AbsTol', 1e-10, 'The epochs'' statistic, from S.trials directly.');
    for c = 1:C
        [tObs, pRef] = local_brute(D(c, S.use(c, :)), tail{1});
        tc.verifyEqual(T.t(c), tObs, 'RelTol', 1e-10);
        tc.verifyEqual(T.p(c), pRef, 'AbsTol', 1e-12, sprintf('contact %d, tail %s', c, tail{1}));
        tc.verifyEqual(T.n_trials(c), nnz(S.use(c, :)));
    end
end
end

function test_few_epochs_among_many_are_exact(tc)
% 20 epochs (random patterns) but a contact that averages only 5: its test
% uses all 32 patterns of its own epochs, and its smallest possible p (2/32)
% cannot reach 0.05, so it is not testable.
rng(9);
S = local_S(3, 20);
S.trials(:, :, S.times >= 0, :) = S.trials(:, :, S.times >= 0, :) * 3;
S.use(2, 6:end) = false;
[T, R] = ieeglab_band_stats(S, 'window', [0 400], 'baseline', [-500 0], 'n_perm', 2000, 'verbose', false);
tc.verifyEqual(R.exact', [false true false]);
tc.verifyEqual(R.n_patterns', [2001 32 2001]);
D = local_d_ref(S, [0 400], [-500 0]);
[~, pRef] = local_brute(D(2, 1:5), 'both');
tc.verifyEqual(T.p(2), pRef, 'AbsTol', 1e-12);
tc.verifyEqual(T.p_min', [1/2001 2/32 1/2001], 'AbsTol', 1e-15);
tc.verifyEqual(T.testable', [true false true]);
end

function test_statistic_is_log_band_power_change(tc)
% Noise-free power falling with frequency (40 to 0 dB), the same in every
% epoch: power x 10 at the lowest band frequency after time 0 is +10 dB there
% and 0 dB at the four others, so the statistic (mean over frequencies) is 2 dB;
% the log of the band's mean power would give 9.6 dB. Identical epochs: t is
% infinite and p is the smallest possible, 2 / 2^12, not 0.
K = 12;
S = local_S(2, K);
S.trials = repmat(single(reshape(10.^(4:-1:0), 1, [], 1, 1)), 2, 1, numel(S.times), K);
S.trials(1, 1, S.times >= 0, :) = S.trials(1, 1, S.times >= 0, :) * 10;
T = ieeglab_band_stats(S, 'window', [0 400], 'baseline', [-500 0], 'verbose', false);
tc.verifyEqual(T.band_db(1), 2, 'AbsTol', 1e-9);
tc.verifyEqual(T.band_db(2), 0, 'AbsTol', 1e-9);
tc.verifyEqual(T.t(1), Inf);
tc.verifyEqual(T.p(1), 2 / 2^K, 'AbsTol', 1e-15);
tc.verifyEqual(T.p(2), 1);
end

function test_missing_frequencies_do_not_bias(tc)
% Power falls with frequency and does not change. As with 'cycles', the lowest
% frequencies have no value near the end of the window and some values are
% blanked: the statistic must stay 0 (a mean over the frequencies left at each
% time would follow the 1/f slope).
K = 10;
S = local_S(1, K);
S.trials = repmat(single(reshape(1 ./ S.freqs.^2, 1, [], 1, 1)), 1, 1, numel(S.times), K);
S.trials(1, 1:2, S.times > 300, :) = NaN;
S.trials(1, 3, S.times < -300, 1:4) = NaN;
S.window_s = [0.3 0.25 0.1 0.1 0.1];
[T, R] = ieeglab_band_stats(S, 'window', [0 400], 'baseline', [-500 0], 'verbose', false);
tc.verifyEqual(T.band_db, 0, 'AbsTol', 1e-6);
tc.verifyEqual(R.baseline_mask(1, :), S.times >= -500 + 149.5 & S.times <= -149.5, ...
    'Per-frequency baseline times, as ieeglab_spectral (half-window (L - 1) / 2 samples).');
end

function test_detects_injected_hfa_only_where_it_is(tc)
% End to end from epochs: a 100 Hz burst 100-400 ms after time 0 on two
% contacts of ten, 40 epochs. Those two come out significant (increase), the
% others not.
fs = 1000; T = 1500; N = 40; C = 10; rng(5);
t = (0:T-1) / fs - 0.5;
X = randn(C, T, N);
burst = 0.8 * sin(2 * pi * 100 * t) .* (t >= 0.1 & t <= 0.4);
X(1:2, :, :) = X(1:2, :, :) + reshape(burst, 1, T);
EEG = local_eeg(X, fs, -0.5);
S = ieeglab_spectral(EEG, 'freqs', 70:10:150, 'window', 0.1, 'keeptrials', true, 'verbose', false);
Tb = ieeglab_band_stats(S, 'window', [100 400], 'baseline', [-500 0], 'n_perm', 2000, 'verbose', false);
tc.verifyEqual(find(Tb.significant)', [1 2]);
tc.verifyTrue(all(Tb.band_db(1:2) > 1));

% Single-trial power of the band's frequencies only: the same test
S2 = ieeglab_spectral(EEG, 'freqs', 40:10:150, 'window', 0.1, 'keeptrials', 'band', 'verbose', false);
tc.verifyEqual(S2.trials_freqs, 70:10:150);
tc.verifyEqual(size(S2.trials, 2), 9);
T2 = ieeglab_band_stats(S2, 'window', [100 400], 'baseline', [-500 0], 'n_perm', 2000, 'verbose', false);
tc.verifyEqual(T2.p, Tb.p);
tc.verifyEqual(T2.band_db, Tb.band_db, 'AbsTol', 1e-12);
tc.verifyError(@() ieeglab_band_stats(S2, 'band', [40 150], 'verbose', false), 'ieeglab_band_stats:bandNotKept');

% From the menu function: the table is kept, the single-trial power is not
[~, E2, com] = evalc(['pop_ieeglab_spectral(EEG, ''type'', ''timefreq'', ''freqs'', 70:10:150, ''window'', 0.1, ' ...
    '''band_window'', [100 400], ''stats'', {''n_perm'', 2000}, ''plot'', false)']);
R = E2.ieeglab.spectral;
tc.verifyFalse(isfield(R, 'trials'));
tc.verifyEqual(R.stats.p, Tb.p);
tc.verifyEqual(R.stats_opt.window, [100 400]);
tc.verifySubstring(com, '''stats'', {''n_perm'', 2000}');
end

function test_null_is_calibrated(tc)
% No effect anywhere: across 400 contacts, uncorrected p < 0.05 for about 5%
% of them (limit 3 SE above), and FDR finds almost nothing.
rng(6);
S = local_S(400, 20);
T = ieeglab_band_stats(S, 'window', [0 400], 'baseline', [-500 0], 'n_perm', 2000, 'verbose', false);
rate = mean(T.p < 0.05);
tc.verifyLessThan(rate, 0.05 + 3 * sqrt(0.05 * 0.95 / 400));
tc.verifyGreaterThan(rate, 0.01);
tc.verifyLessThanOrEqual(nnz(T.significant), 2);

% Corrections against their definitions (p of the 400 contacts, ties included)
p = T.p; n = numel(p);
Tb = ieeglab_band_stats(S, 'window', [0 400], 'baseline', [-500 0], 'n_perm', 2000, 'correct', 'bonferroni', 'verbose', false);
tc.verifyEqual(Tb.p_adj, min(1, p * n), 'AbsTol', 1e-15);
[ps, ord] = sort(p); bh = zeros(n, 1);
for i = 1:n, bh(ord(i)) = min(1, min(ps(i:n) .* n ./ (i:n)')); end
tc.verifyEqual(T.p_adj, bh, 'AbsTol', 1e-12);
Ty = ieeglab_band_stats(S, 'window', [0 400], 'baseline', [-500 0], 'n_perm', 2000, 'correct', 'by', 'verbose', false);
by = zeros(n, 1);
for i = 1:n, by(ord(i)) = min(1, min(ps(i:n) .* n ./ (i:n)' * sum(1 ./ (1:n)))); end
tc.verifyEqual(Ty.p_adj, by, 'AbsTol', 1e-12);
end

function test_input_errors(tc)
% Wrong 'keeptrials' or 'stats' values are refused with their own message; a
% test band of its own keeps the trials it needs.
rng(10);
EEG = local_eeg(randn(3, 1000, 12), 1000, -0.5);
tc.verifyError(@() ieeglab_spectral(EEG, 'freqs', 70:10:150, 'keeptrials', [], 'verbose', false), 'ieeglab_spectral:keeptrials');
tc.verifyError(@() ieeglab_spectral(EEG, 'freqs', 70:10:150, 'keeptrials', 'all', 'verbose', false), 'ieeglab_spectral:keeptrials');
tc.verifyError(@() pop_ieeglab_spectral(EEG, 'freqs', 70:10:150, 'stats', 'on', 'plot', false), 'pop_ieeglab_spectral:stats');
tc.verifyError(@() pop_ieeglab_spectral(EEG, 'freqs', 70:10:150, 'band', [], 'stats', true, 'plot', false), ...
    'pop_ieeglab_spectral:statsNeedBand');
[~, E] = evalc(['pop_ieeglab_spectral(EEG, ''freqs'', 30:10:150, ''window'', 0.1, ''stats'', ' ...
    '{''band'', [30 60], ''window'', [0 300], ''baseline'', [-400 0], ''n_perm'', 500}, ''plot'', false)']);
tc.verifyEqual(E.ieeglab.spectral.stats_opt.band, [30 60]);
tc.verifyFalse(isfield(E.ieeglab.spectral, 'trials'));
end

function test_keeptrials_band_with_cycles(tc)
% Windows that grow at low frequencies ('cycles'): single-trial power of the
% band only gives the same test as all of it.
rng(11);
EEG = local_eeg(randn(4, 1500, 15), 1000, -0.5);
S1 = ieeglab_spectral(EEG, 'freqs', 4:2:30, 'cycles', 3, 'band', [4 12], 'keeptrials', true, 'verbose', false);
S2 = ieeglab_spectral(EEG, 'freqs', 4:2:30, 'cycles', 3, 'band', [4 12], 'keeptrials', 'band', 'verbose', false);
o = {'window', [0 600], 'baseline', [-500 0], 'per_time', true, 'n_perm', 500, 'verbose', false};
[T1, R1] = ieeglab_band_stats(S1, o{:});
[T2, R2] = ieeglab_band_stats(S2, o{:});
tc.verifyTrue(isequaln(T1, T2));
tc.verifyTrue(isequaln(R1.p_tc, R2.p_tc));
end

function test_unused_epochs_and_bad_contacts(tc)
% Epochs a contact does not average (S.use, its own stimulation on CCEP data)
% are left out of its test; bad contacts and contacts with too few epochs are
% not tested and do not enter the correction.
S = local_S(3, 12);
S.use(1, 1:4) = false;
S.use(3, 1:9) = false;
S.bad(2) = true;
T = ieeglab_band_stats(S, 'window', [0 400], 'baseline', [-500 0], 'verbose', false);
tc.verifyEqual(T.n_trials(1), 8);
tc.verifyTrue(isnan(T.p(2)) && isnan(T.p(3)));
tc.verifyEqual(T.p_adj(1), T.p(1), 'One contact tested: nothing to correct.');
end

function test_per_time(tc)
% Time-resolved test: significant only where the effect is.
rng(7);
S = local_S(2, 30);
late = S.times >= 200 & S.times <= 400;
S.trials(1, :, late, :) = S.trials(1, :, late, :) * 4;
[~, R] = ieeglab_band_stats(S, 'window', [0 400], 'baseline', [-500 0], 'per_time', true, 'n_perm', 2000, 'verbose', false);
sig = R.sig_tc(1, :);
tc.verifyTrue(all(sig(R.times_tc >= 220 & R.times_tc <= 380)));
tc.verifyFalse(any(sig(R.times_tc <= 150)));
tc.verifyFalse(any(R.sig_tc(2, :)));
end

function test_fieldtrip_import(tc)
% FieldTrip raw data and elec in cm with the labels in another order: data
% kept, coordinates matched by label and converted to mm, one event at time 0
% per epoch carrying its trialinfo row, types from patterns.
fs = 200; rng(8);
time = -0.5:1/fs:1;
data = struct('label', {{'G1'; 'G2'; 'D1'}}, 'fsample', fs, 'trialinfo', [500; 650; 720]);
for k = 1:3, data.trial{k} = randn(3, numel(time)); data.time{k} = time; end
elec = struct('label', {{'D1'; 'G2'; 'G1'}}, 'elecpos', [1 2 3; 4 5 6; 7 8 9], 'unit', 'cm', 'coordsys', 'acpc');
EEG = ieeglab_import_fieldtrip(data, 'elec', elec, 'types', {'G*', 'ECOG'; '*', 'SEEG'}, 'verbose', false);
tc.verifyEqual(size(EEG.data), [3 numel(time) 3]);
tc.verifyEqual(double(EEG.data(:, :, 2)), data.trial{2}, 'AbsTol', 1e-5);
tc.verifyEqual([EEG.chanlocs(1).X EEG.chanlocs(1).Y EEG.chanlocs(1).Z], [70 80 90], 'AbsTol', 1e-9);   % G1 (eeg_checkset rounds through spherical)
tc.verifyEqual([EEG.chanlocs(3).X EEG.chanlocs(3).Y EEG.chanlocs(3).Z], [10 20 30], 'AbsTol', 1e-9);   % D1
tc.verifyEqual({EEG.chanlocs.type}, {'ECOG', 'ECOG', 'SEEG'});
tc.verifyEqual(EEG.ieeglab.coordsys, 'acpc');
tc.verifyEqual(numel(EEG.event), 3);
tc.verifyEqual(EEG.event(3).trialinfo, 720);
tc.verifyEqual(EEG.times(EEG.event(2).latency - EEG.pnts), 0, 'AbsTol', 1e-9);
tc.verifyEqual(ieeglab_detect_mode(EEG), 'erp');
[~, info] = ieeglab_brain_surfaces(EEG, 'check_fit', false);
tc.verifyEqual(info.space, 'acpc');
end

% ===================== helpers =====================

function S = local_S(C, K)
% A minimal ieeglab_spectral time-frequency output with single-trial power:
% 5 frequencies in 70-150 Hz, times -500..800 ms every 20 ms, lognormal power.
times = -500:20:800;
S = struct('type', 'timefreq', 'labels', {compose('C%d', 1:C)}, 'bad', false(C, 1), 'srate', 1000, ...
    'freqs', [70 90 110 130 150], 'times', times, 'window_s', repmat(0.101, 1, 5), ...
    'band', [70 150], 'band_window', [0 800], 'baseline', [-500 0], 'use', true(C, K));
S.trials = single(exp(0.5 * randn(C, 5, numel(times), K)));
end

function D = local_d_ref(S, win, base)
% Each epoch's statistic, written out: per contact, epoch and frequency, the
% mean of 10 log10(power) over the window minus its mean over the baseline
% times whose whole window fits (half-window (L - 1) / 2 samples), then the
% mean over frequencies.
[C, F, ~, N] = size(S.trials);
D = nan(C, N);
for c = 1:C
    for n = 1:N
        x = nan(1, F);
        for f = 1:F
            h = (S.window_s(f) * S.srate - 1) / 2 / S.srate * 1000;
            b = S.times >= base(1) + h & S.times <= base(2) - h;
            w = S.times >= win(1) & S.times <= win(2);
            L = 10 * log10(double(squeeze(S.trials(c, f, :, n))))';
            x(f) = mean(L(w), 'omitnan') - mean(L(b), 'omitnan');
        end
        D(c, n) = mean(x, 'omitnan');
    end
end
end

function [tObs, p] = local_brute(d, tail)
% One-sample t and its p over every sign pattern of d, recomputed directly.
K = numel(d);
B = dec2bin(0:2^K-1, K) == '1';
tp = zeros(2^K, 1);
for j = 1:2^K
    x = d .* (1 - 2 * B(j, :));
    tp(j) = mean(x) / (std(x) / sqrt(K));
end
tObs = mean(d) / (std(d) / sqrt(K));
switch tail
    case 'both',     p = mean(abs(tp) >= abs(tObs) - 1e-12);
    case 'positive', p = mean(tp >= tObs - 1e-12);
    case 'negative', p = mean(-tp >= -tObs - 1e-12);
end
end

function EEG = local_eeg(X, fs, xmin)
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
i0 = find(EEG.times >= 0, 1);
for n = 1:N
    EEG.event(n) = struct('type', 'tone', 'latency', (n - 1) * T + i0, 'duration', 0);
    EEG.epoch(n).event = n; EEG.epoch(n).eventtype = {'tone'}; EEG.epoch(n).eventlatency = {0};
end
end
