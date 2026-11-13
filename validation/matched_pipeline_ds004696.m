function S = matched_pipeline_ds004696(dsDir, outDir, only)
% matched_pipeline_ds004696() - The published CRP decisions of ds004696 sub-02
%       (Ojeda Valencia et al., 2023) reproduced with iEEGLAB at the paper's
%       settings, then each setting changed, one at a time, to iEEGLAB's.
%
% Usage:
%   S = matched_pipeline_ds004696('C:\data\ds004696', 'C:\data\ref_compare')
%   S = matched_pipeline_ds004696(dsDir, outDir, {'paper'})   % some variants only
%
% dsDir    the ds004696 BIDS folder with sub-02 and derivatives/stats/sub-02
% outDir   the folder of compare_references_ds004696: its preprocessed epochs
%          and its 'hapwave' and 'carla' tables are reused (and made first when
%          missing). Writes one table per variant (matched_<variant>.tsv: site,
%          channel, trials, tau_R, p) and matched_summary.json. A variant whose
%          table exists is not run again, so the script resumes after an
%          interruption.
% only     variants to run (default: all); the summary covers every table present
%
% The paper's settings ('paper'): no filtering, stimulation blanking or
% line-noise removal; contacts marked bad in channels.tsv kept out of the
% reference and the tests, every other contact kept; no trial rejection;
% epochs -1 to 1 s with a per-trial median baseline strictly between -500 and
% -50 ms (iEEGLAB takes the samples nearest to the window's ends, so the
% window is given as [-499.5 -50.3] at 2048 Hz); the paper's re-referencing,
% before the baseline (in each block of 64 channels of the recording, the 20%
% of good contacts with the lowest variance over the trials from 15 to 500 ms,
% 'varsubset'; each trial's offset counts in that variance); CRP from 15 ms to
% 1 s, p at tau_R (t-test), and the
% Benjamini-Yekutieli correction over the paper's 509 pairs (21 sites). iEEGLAB
% leaves out site RF4-RF5, which stimulates a bad contact; the published p of
% its 25 pairs complete the paper's family, so the correction runs over the
% same 509 p values.
%
% 'paper_epochs' is the same with the published epochs' timing: they place the
% stimulation sample of events.tsv one sample after time zero, iEEGLAB at zero
% (the BIDS onset); every other variant uses iEEGLAB's.
%
% Each other variant changes one setting to the one of Part A of the
% step-by-step page:
%   baseline     baseline from -500 to -50 ms, both ends included (two samples more)
%   contacts     contacts without coordinates and bad contacts removed
%   line_noise   line noise at 60, 120 and 180 Hz removed (sine fit)
%   blanking     stimulation artifact blanked
%   rejection    outlier trials removed
%   order        re-referencing after the baseline, as the re-referencing step
%                of the step-by-step page (which then removes the baseline again)
%   carla        CARLA instead of the paper's re-referencing
%   permutation  permutation test instead of the t-test at tau_R
%   bh           Benjamini-Hochberg instead of Benjamini-Yekutieli
%   all_pairs    the correction over all pairs of the recording, not the paper's
% (with about 10 trials per site the sign-flip permutation test cannot give p
% below 2/2^K, 0.002 for 10 trials, which limits how many pairs pass
% Benjamini-Yekutieli over 509 pairs), and five rows combine them:
%   permutation_bh   the permutation test with Benjamini-Hochberg
%   statistics       iEEGLAB's statistics (permutation, Benjamini-Hochberg over
%                    all pairs of the recording) on the paper's preprocessing
%                    and re-referencing
%   preprocessing    Part A's preprocessing (the first five changes), then the
%                    paper's re-referencing after the baseline
%   part_a_paper_ref Part A with the paper's re-referencing (permutation test,
%                    Benjamini-Hochberg over all pairs of the recording)
%   part_a           Part A as it is (CARLA), iEEGLAB's defaults
% Every variant is compared with the paper's decisions on the pairs both
% analyses test: same decision, Cohen's kappa, the 2 x 2 table, and tau_R.
%
% References:
%   Huang, H., et al. (2024). CARLA: Adjusted common average referencing for
%     cortico-cortical evoked potential data. Journal of Neuroscience Methods,
%     407, 110153.
%   Ojeda Valencia, G., et al. (2023). Signatures of electrical stimulation
%     driven network interactions in the human limbic system. Journal of
%     Neuroscience, 43(39), 6697-6711.
%
% Cedric Cannard, iEEGLAB, 2026

t00 = tic;
stem = 'sub-02_ses-ieeg01_task-ccep_run-01';
rec = fullfile(dsDir, 'sub-02', 'ses-ieeg01', 'ieeg', [stem '_ieeg.mefd']);
P = load(fullfile(dsDir, 'derivatives', 'stats', 'sub-02', [stem '_crp.mat']));
pub = local_published(P);
fprintf('Published CRP: %d pairs at %d sites, %d respond (Benjamini-Yekutieli over the %d pairs)\n', ...
    numel(pub.key), numel(P.average_ccep_names), nnz(pub.significant), numel(pub.key));
if ~all(cellfun(@isfile, fullfile(outDir, {'preprocessed.mat', 'pairs_hapwave.tsv', 'pairs_carla.tsv'})))
    compare_references_ds004696(dsDir, outDir, {'hapwave', 'carla'});
end

% ---------------- variants: name, change, preprocessing, reference, test ----------------
paperPre = struct('apply_blank', false, 'apply_notch', false, 'apply_highpass', false, ...
    'apply_epoch', true, 'epoch_window', [-1000 1000], 'reject_trials', false, 'remove_rare_cond', false, ...
    'remove_no_coords', false, 'remove_bad_channels', false, 'apply_baseline', true, ...
    'baseline_method', 'median', 'baseline_period', [-499.5 -50.3], 'plot', false);   % nearest samples: the open interval
% reference: re-referencing inside the preprocessing, before the baseline, as the
% paper (its contacts are ranked by their variance over all trials, so each
% trial's offset counts); '_after': on the baseline-corrected epochs, as the
% re-referencing step of the step-by-step page (which then removes the baseline again)
V = {
    'paper_epochs', 'the paper''s epochs (stimulation sample at +1/srate)', struct('event_shift', -1),                 'varsubset', 'ttest'
    'paper',       'none (the paper''s settings)',                          struct(),                                   'varsubset', 'ttest'
    'baseline',    'baseline [-500 -50] ms (iEEGLAB: both ends included)',  struct('baseline_period', [-500 -50]),      'varsubset', 'ttest'
    'contacts',    'contacts without coordinates and bad contacts removed', struct('remove_no_coords', true, ...
                                                                                   'remove_bad_channels', true),   'varsubset', 'ttest'
    'blanking',    'stimulation artifact blanked',                          struct('apply_blank', true),                'varsubset', 'ttest'
    'rejection',   'outlier trials removed',                                struct('reject_trials', true),              'varsubset', 'ttest'
    'line_noise',  'line noise at 60, 120 and 180 Hz removed',              struct('apply_notch', true, ...
                                                                                   'notch', [60 120 180]),         'varsubset', 'ttest'
    'order',       're-referencing after the baseline',                     struct(),                                   'varsubset_after', 'ttest'
    'carla',       'CARLA instead of the paper''s re-referencing',          struct(),                                   'carla',     'ttest'
    'permutation', 'permutation test instead of the t-test at tau_R',       struct(),                                   'varsubset', 'permutation'
    };
if nargin < 3 || isempty(only), only = [V(:, 1)', {'preprocessing'}]; end
only = cellstr(only);

EEGc = []; EEGp = []; oP = [];
for v = 1:size(V, 1)
    if ~ismember(V{v, 1}, only), continue; end
    fT = fullfile(outDir, ['matched_' V{v, 1} '.tsv']);
    if isfile(fT), fprintf('%s: done before\n', V{v, 1}); continue; end
    t0 = tic;
    o = paperPre; d = V{v, 3};
    for f = fieldnames(d)', o.(f{1}) = d.(f{1}); end
    ref = V{v, 4}; after = endsWith(ref, '_after');
    if ~after, o.apply_car = true; o.car_method = ref; o.car_fraction = 0.2; end
    if ~isempty(EEGp) && isequal(o, oP)
        EEG = EEGp;
    else
        if isempty(EEGc)
            [~, EEGc] = evalc('ieeglab_import(rec)'); %#ok<ASGLU> used in the next evalc
            [~, EEGc] = evalc('ieeglab_load(EEGc, struct(''event_field'', ''electrical_stimulation_site''))');
            fprintf('Imported: %d contacts, %d events (%.0f s)\n', EEGc.nbchan, numel(EEGc.event), toc(t0));
        end
        oP = o; EEGs = EEGc;
        if isfield(o, 'event_shift')     % move every event by this many samples (the data are not copied)
            lat = num2cell([EEGs.event.latency] + o.event_shift); [EEGs.event.latency] = lat{:};
            o = rmfield(o, 'event_shift'); %#ok<NASGU> used in the evalc below
        end
        [~, EEG] = evalc('ieeglab_preprocess(EEGs, o)'); clear EEGs
        EEGp = EEG;
    end
    if after, EEG = local_reref(EEG, erase(ref, '_after')); end
    tPre = toc(t0);
    t1 = tic; T = local_stats(EEG, V{v, 5});
    writetable(T, fT, 'FileType', 'text', 'Delimiter', '\t');
    fprintf('%s: %d contacts x %d epochs, preprocessing and re-referencing %.0f s, CRP %.0f s, %d pairs\n', ...
        V{v, 1}, EEG.nbchan, EEG.trials, tPre, toc(t1), height(T));
    clear EEG
end
clear EEGc EEGp

fT = fullfile(outDir, 'matched_preprocessing.tsv');
if ismember('preprocessing', only) && ~isfile(fT)
    t0 = tic;
    L = load(fullfile(outDir, 'preprocessed.mat'), 'EEG'); EEG = L.EEG; clear L
    if ~isfield(EEG.etc.ieeglab_import, 'recording_channels')    % saved before the import recorded it
        I = ieeglab_import(rec, 'info', true);
        EEG.etc.ieeglab_import.recording_channels = I.labels;
    end
    EEG = local_reref(EEG, 'varsubset');
    T = local_stats(EEG, 'ttest');
    writetable(T, fT, 'FileType', 'text', 'Delimiter', '\t');
    fprintf('preprocessing: %d contacts x %d epochs, %.0f s, %d pairs\n', EEG.nbchan, EEG.trials, toc(t0), height(T));
    clear EEG
end

% ---------------- summary: every table present ----------------
R = [V(:, 1:2), strcat('matched_', V(:, 1)), repmat({'paper', 'by'}, size(V, 1), 1)
    {   % row, change, table, correction family, method
    'bh',               'Benjamini-Hochberg instead of Benjamini-Yekutieli',              'matched_paper',         'paper', 'bh'
    'all_pairs',        'correction over all pairs of the recording',                     'matched_paper',         'all',   'by'
    'permutation_bh',   'permutation test with Benjamini-Hochberg',                       'matched_permutation',   'paper', 'bh'
    'statistics',       'iEEGLAB''s statistics (permutation, BH over all pairs)',         'matched_permutation',   'all',   'bh'
    'preprocessing',    'Part A''s preprocessing, re-referenced after the baseline',      'matched_preprocessing', 'paper', 'by'
    'part_a_paper_ref', 'Part A with the paper''s re-referencing',                        'pairs_hapwave',         'all',   'bh'
    'part_a',           'Part A as it is (CARLA, iEEGLAB''s defaults)',                   'pairs_carla',           'all',   'bh'
    }];
S = struct([]);
for r = 1:size(R, 1)
    f = fullfile(outDir, [R{r, 3} '.tsv']);
    if ~isfile(f), continue; end
    T = readtable(f, 'FileType', 'text', 'Delimiter', '\t');
    T = T(isfinite(T.p), :);
    c = local_compare(T, pub, R{r, 4}, R{r, 5});
    s = struct('variant', R{r, 1}, 'change', R{r, 2}, 'correction', ...
        sprintf('%s over %s', upper(R{r, 5}), local_if(strcmp(R{r, 4}, 'paper'), 'the paper''s pairs', 'all pairs')));
    for f = fieldnames(c)', s.(f{1}) = c.(f{1}); end
    S = [S, s]; %#ok<AGROW>
end
fid = fopen(fullfile(outDir, 'matched_summary.json'), 'w');
fprintf(fid, '%s', jsonencode(S, 'PrettyPrint', true)); fclose(fid);
disp(struct2table(rmfield(S, {'change', 'correction'})));
fprintf('Done in %.0f min.\n', toc(t00) / 60);
end

function EEG = local_reref(EEG, method)
switch method
    case 'varsubset', [~, EEG] = evalc('pop_ieeglab_reref(EEG, ''method'', ''varsubset'', ''fraction'', 0.2)');
    case 'carla',     [~, EEG] = evalc('pop_ieeglab_reref(EEG, ''method'', ''carla'')');
end
end

function T = local_stats(EEG, test) %#ok<INUSL> EEG is used in the evalc below
% CRP of every site-contact pair; the correction is applied later, per family
o = struct('run_n1', false, 'run_crp', true, 'run_matrix', false, 'plot', false, 'verbose', false, ...
    'crp_test', test, 'correct', 'by', 'correct_scope', 'subject', 'crp_window', [15 1000]); %#ok<NASGU>
[~, E] = evalc('ieeglab_stats_subject(EEG, o)');
T = E.ieeglab.stats.table(:, {'site', 'channel', 'n_trials', 'tR_ms', 'p'});
end

function s = local_compare(T, pub, family, method)
% Decisions of T (p adjusted over the family) against the paper's, on the pairs both test
key = string(T.site) + "|" + string(T.channel);
[inT, j] = ismember(key, pub.key);
if strcmp(family, 'paper')       % the paper's 509 p values, iEEGLAB's where it tests the pair
    [inP, i] = ismember(pub.key, key);
    pf = pub.p; pf(inP) = T.p(i(inP));
    q = local_adjust(pf, method);
    padj = nan(height(T), 1); padj(inT) = q(j(inT));
    nFam = numel(pf);
else
    padj = local_adjust(T.p, method); nFam = height(T);
end
a = padj(inT) < 0.05; b = pub.significant(j(inT)) == 1;
n = nnz(inT); po = mean(a == b); pe = mean(a) * mean(b) + mean(~a) * mean(~b);
dt = abs(T.tR_ms(inT) - pub.tR_ms(j(inT)));
lp = abs(log10(max(T.p(inT), realmin)) - log10(max(pub.p(j(inT)), realmin)));
s = struct('pairs_tested', height(T), 'family_size', nFam, 'pairs_compared', n, ...
    'same_decision_pct', 100 * po, 'kappa', (po - pe) / (1 - pe), 'respond_both', nnz(a & b), ...
    'respond_only_ieeglab', nnz(a & ~b), 'respond_only_published', nnz(~a & b), 'respond_neither', nnz(~a & ~b), ...
    'tR_same_sample_pct', 100 * mean(dt < 0.2), 'tR_median_abs_diff_ms', median(dt), ...
    'p_max_abs_log10_diff', max(lp));
end

function q = local_adjust(p, method)
% Benjamini-Hochberg ('bh') or Benjamini-Yekutieli ('by') adjusted p values
m = numel(p); [ps, o] = sort(p(:));
c = 1; if strcmp(method, 'by'), c = sum(1 ./ (1:m)); end
qs = min(1, flipud(cummin(flipud(ps * m * c ./ (1:m)'))));
q = zeros(m, 1); q(o) = qs;
end

function pub = local_published(P)
% Every pair of the published CRP: key "site|contact", p at tau_R, the paper's
% decision (Benjamini-Yekutieli over all its pairs, below 0.05), tau_R in ms
key = strings(0, 1); p = []; tR = [];
for s = 1:numel(P.average_ccep_names)
    cs = ieeglab_canonical_site(P.average_ccep_names(s), P.channel_names);
    for ch = 1:size(P.crp_out, 1)
        Q = P.crp_out(ch, s);
        if isempty(Q.crp_parms) || ~isfield(Q.crp_projs, 'p_value_tR'), continue; end
        key(end + 1, 1) = cs(1) + "|" + string(P.channel_names{ch}); %#ok<AGROW>
        p(end + 1, 1) = Q.crp_projs.p_value_tR; %#ok<AGROW>
        tR(end + 1, 1) = Q.crp_parms.tR * 1000; %#ok<AGROW>
    end
end
pub = struct('key', key, 'p', p, 'significant', double(local_adjust(p, 'by') < 0.05), 'tR_ms', tR);
end

function s = local_if(c, a, b)
if c, s = a; else, s = b; end
end
