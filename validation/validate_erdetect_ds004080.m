function out = validate_erdetect_ds004080(dsDir, erdDir)
% validate_erdetect_ds004080() - iEEGLAB's N1 detection against ER-detect
%       (van den Boom et al., 2025) on ds004080 (van Blooijs et al., 2023),
%       sub-ccepAgeUMCU02 run-041456: ECoG, 2048 Hz, single-pulse stimulation.
%
% Usage:
%   out = validate_erdetect_ds004080('C:\data\ds004080', 'C:\data\erdetect_ds004080')
%
% dsDir   the ds004080 BIDS folder (at least sub-ccepAgeUMCU02, with the
%         electrodes, channels and events files)
% erdDir  ER-detect's output for this participant, made with its defaults
%         (version 2.6.2; on Windows set PYTHONUTF8=1 first):
%           python -m erdetect <dsDir> <erdDir> --participant_label ccepAgeUMCU02
%         It holds sub-ccepAgeUMCU02_ses-1_task-SPESclin_run-041456/erdetect_data.mat.
%
% ER-detect's defaults: no filtering or re-referencing, epochs -1 to 2 s,
% median baseline -0.5 to -0.02 s, both directions of a stimulation pair
% together, and a response where the trial average (taken as it is) has a
% negative peak 9 to 90 ms after the pulse larger than 3.4 x max(SD of -1 to
% -0.1 s, 50 uV). Its peaks come from a peak finder with 20 uV selectivity
% over 0 to 500 ms, the largest peak in the window must be negative, and peaks
% above 3000 uV are rejected; iEEGLAB takes the largest local extremum in the
% window (ieeglab_detect_n1 help).
% Three comparisons, each over every stimulation site x recording contact of
% this run (8 sites x 84 contacts):
%   A  same input: ER-detect's own trial averages through ieeglab_detect_n1
%      with ER-detect's rule ('sd', 'peak_from' 'zero'): the detection alone
%   B  iEEGLAB end to end with ER-detect's choices: ieeglab_import,
%      ieeglab_load, ieeglab_preprocess (no filter, epochs -1 to 2 s, median
%      baseline -500 to -20 ms), the same rule
%   C  iEEGLAB's defaults (the step-by-step tutorial): artifact blanking, line
%      noise removed at 50/100/150 Hz, epochs -1 to 1 s, outlier trials
%      removed, median baseline -500 to -50 ms, common average, N1 by sign-flip
%      permutation test with FDR across contacts of each site, window 10 to
%      100 ms (so detections after 90 ms cannot be ER-detect's; counted)
% Reported: agreement (Cohen's kappa) of the detections, and where both detect,
% latency and amplitude. ER-detect 2.6.2 reports latencies one sample early
% (detection.py, peak index), so latencies are also compared after that shift.
% An ER-detect response at a pair iEEGLAB did not test counts as missed.
%
% References:
%   van Blooijs, D., van den Boom, M. A., van der Aar, J. F., Huiskamp, G. M.,
%     Castegnaro, G., Demuru, M., Zweiphenning, W. J. E. M., van Eijsden, P.,
%     Miller, K. J., Leijten, F. S. S., & Hermes, D. (2023). Developmental
%     trajectory of transmission speed in the human brain. Nature Neuroscience,
%     26(4), 537-541.
%   van den Boom, M. A., Gregg, N. M., Ojeda Valencia, G., Lundstrom, B. N., Miller,
%     K. J., van Blooijs, D., Huiskamp, G. J. M., Leijten, F. S. S., Worrell,
%     G. A., & Hermes, D. (2025). ER-detect: A pipeline for robust detection of
%     early evoked responses in BIDS-iEEG electrical stimulation data. Journal of
%     Neuroscience Methods, 418, 110389.
%
% Cedric Cannard, iEEGLAB, 2026

t0 = tic;
run = 'sub-ccepAgeUMCU02_ses-1_task-SPESclin_run-041456';
erdFile = fullfile(erdDir, run, 'erdetect_data.mat');
if ~isfile(erdFile)
    error('validate_erdetect_ds004080:noErdetect', ['No %s. Run ER-detect first:\n' ...
        '  python -m erdetect %s %s --participant_label ccepAgeUMCU02'], erdFile, dsDir, erdDir);
end
E = load(erdFile);
ref.channels = cellstr(string(E.channel_labels(:)));
ref.sites = cellstr(string(E.stimpair_labels(:)));
ref.latency = double(E.neg_peak_latency_ms);           % channels x sites, NaN = no response
ref.amplitude = double(E.neg_peak_amplitudes);
fs = double(E.sampling_rate);
rule = struct('method', 'sd', 'threshold', 3.4, 'min_baseline_sd', 50, 'polarity', 'negative', ...
    'peak_from', 'zero', 'n1_window', [9 90], 'baseline', [-1000 -100], 'min_trials', 2, 'verbose', false);
out.erdetect_responses = nnz(isfinite(ref.latency));

% ---------------- A: same input ----------------
nC = numel(ref.channels); nS = numel(ref.sites);
X = zeros(nC, size(E.ccep_average, 3), 2 * nS);
for s = 1:nS, X(:, :, 2*s-1:2*s) = repmat(reshape(E.ccep_average(:, s, :), nC, []), 1, 1, 2); end
EEG = eeg_emptyset;
EEG.data = X; [EEG.nbchan, EEG.pnts, EEG.trials] = size(X);
EEG.srate = fs; EEG.times = double(E.epoch_time_s(:))' * 1000;
EEG.xmin = EEG.times(1) / 1000; EEG.xmax = EEG.times(end) / 1000;
EEG.chanlocs = struct('labels', ref.channels', 'type', repmat({'ECOG'}, 1, nC));
i0 = find(EEG.times >= 0, 1);
EEG.event = struct('type', repelem(ref.sites', 2), 'latency', num2cell(((1:EEG.trials) - 1) * EEG.pnts + i0), ...
    'epoch', num2cell(1:EEG.trials));
EEG = eeg_checkset(EEG, 'eventconsistency');
[~, T] = ieeglab_detect_n1(EEG, rule);
out.A = local_compare('A  same input (ER-detect''s trial averages), ER-detect''s rule', T, ref, fs);

% ---------------- B and C: iEEGLAB from the raw recording ----------------
rec = fullfile(dsDir, 'sub-ccepAgeUMCU02', 'ses-1', 'ieeg', [run '_ieeg.vhdr']);
EEG0 = ieeglab_import(rec);
EEG0 = ieeglab_load(EEG0, struct('event_field', 'electrical_stimulation_site'));

EEG = ieeglab_preprocess(EEG0, struct('apply_epoch', true, 'epoch_window', [-1000 2000], ...
    'remove_bad_channels', true, 'apply_baseline', true, 'baseline_method', 'median', ...
    'baseline_period', [-500 -20], 'plot', false));
[~, T] = ieeglab_detect_n1(EEG, rule);
out.B = local_compare('B  iEEGLAB end to end, ER-detect''s choices and rule', T, ref, fs);

EEG = ieeglab_preprocess(EEG0, struct('apply_blank', true, 'apply_notch', true, 'notch', [50 100 150], ...
    'apply_epoch', true, 'epoch_window', [-1000 1000], 'reject_trials', true, 'remove_rare_cond', true, ...
    'min_trials', 5, 'remove_bad_channels', true, 'apply_baseline', true, 'baseline_period', [-500 -50], 'plot', false));
EEG = pop_ieeglab_reref(EEG, 'method', 'car');
EEG = ieeglab_stats_subject(EEG, struct('run_n1', true, 'run_crp', false, 'run_matrix', false, 'plot', false));
out.C = local_compare('C  iEEGLAB defaults, permutation test', EEG.ieeglab.n1.table, ref, fs);
out.seconds = toc(t0);
fprintf('Done in %.0f s.\n', out.seconds);
end

function r = local_compare(name, T, ref, fs)
% Site x contact pairs present in both; detections, then latency and amplitude
[inC, ic] = ismember(T.channel, ref.channels);
[inS, is] = ismember(T.site, ref.sites);
k = inC & inS;
r.name = name;
r.not_in_erdetect = height(T) - nnz(k);
T = T(k, :); idx = sub2ind(size(ref.latency), ic(k), is(k));
dI = logical(T.significant); dE = isfinite(ref.latency(idx));
% ER-detect responses at pairs iEEGLAB did not test count as missed
untested = find(isfinite(ref.latency(:)));
untested = setdiff(untested, idx);
r.erdetect_untested_by_ieeglab = numel(untested);
dI = [dI; false(numel(untested), 1)]; dE = [dE; true(numel(untested), 1)];
n = numel(dI);
a = nnz(dI & dE); b = nnz(dI & ~dE); c = nnz(~dI & dE); d = nnz(~dI & ~dE);
po = (a + d) / n; pe = ((a + b) * (a + c) + (c + d) * (b + d)) / n^2;
r.n_pairs = n; r.both = a; r.ieeglab_only = b; r.erdetect_only = c; r.neither = d;
r.kappa = (po - pe) / (1 - pe);
fprintf('\n%s\n   %d pairs: both %d, iEEGLAB only %d, ER-detect only %d, neither %d; kappa %.3f\n', ...
    name, n, a, b, c, d, r.kappa);
if r.not_in_erdetect > 0, fprintf('   (%d pairs not in ER-detect''s output left out)\n', r.not_in_erdetect); end
if r.erdetect_untested_by_ieeglab > 0
    fprintf('   (%d ER-detect responses at pairs iEEGLAB did not test, counted as missed)\n', r.erdetect_untested_by_ieeglab);
end
onlyI = dI(1:height(T)) & ~dE(1:height(T));
r.ieeglab_only_after_90ms = nnz(onlyI & T.n1_latency_ms > 90);
if r.ieeglab_only_after_90ms > 0
    fprintf('   %d of the %d iEEGLAB-only detections are after 90 ms (outside ER-detect''s window)\n', ...
        r.ieeglab_only_after_90ms, nnz(onlyI));
end
dI = dI(1:height(T)); dE = dE(1:height(T));                % per-pair details: tested pairs only
both = dI & dE;
if a > 0
    dl = T.n1_latency_ms(both) - ref.latency(idx(both));
    ds = dl - 1000 / fs;                                   % after ER-detect's one-sample shift
    da = T.n1_amplitude_uv(both) - ref.amplitude(idx(both));
    r.latency_same_after_shift = nnz(abs(ds) < 1e-6);
    r.latency_within_5ms = nnz(abs(dl) <= 5);
    r.latency_max_abs_diff_ms = max(abs(dl));
    cc = corrcoef(T.n1_amplitude_uv(both), ref.amplitude(idx(both))); r.amplitude_r = cc(1, 2);
    r.amplitude_median_abs_diff_uv = median(abs(da));
    r.amplitude_max_abs_diff_uv = max(abs(da));
    fprintf(['   both detect: latency equal after the one-sample shift %d/%d, within 5 ms %d/%d (max |diff| %.2f ms); ' ...
        'amplitude r %.4f, median |diff| %.2f uV, max %.2f uV\n'], r.latency_same_after_shift, a, r.latency_within_5ms, a, ...
        r.latency_max_abs_diff_ms, r.amplitude_r, r.amplitude_median_abs_diff_uv, r.amplitude_max_abs_diff_uv);
end
dis = find(dI ~= dE);
r.disagree = T(dis, intersect({'site', 'channel', 'n1_amplitude_uv', 'n1_latency_ms', 'z', 'p_adj'}, ...
    T.Properties.VariableNames, 'stable'));
r.disagree.erdetect_latency_ms = ref.latency(idx(dis));
r.disagree.erdetect_amplitude_uv = ref.amplitude(idx(dis));
if ~isempty(dis) && numel(dis) <= 12, disp(r.disagree); end
end
