function tests = test_ieeglab_regressions
% Regression tests for the issues fixed in v1.1.0.
% Headless and compute-only: nothing here opens a figure or a dialog.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(mfilename('fullpath')));
tc.assumeNotEmpty(which('eeglab'), 'EEGLAB is not on the MATLAB path.');
if isempty(which('pop_loadset')), evalc('eeglab nogui'); end
addpath(root); addpath(fullfile(root, 'functions'));
tc.TestData.root = root;
d = fullfile(root, 'tutorial', 'seeg');
EEG = pop_loadset('filename','sub-02_ses-ieeg01_task-ccep_run-01_ieeg.set','filepath',d);
[~, EEG] = evalc('ieeglab_load(EEG, struct())');
tc.TestData.C = EEG;                                   % continuous, annotated
tc.TestData.opt = struct('apply_highpass',true,'highpass',0.5,'apply_notch',false, ...
    'apply_lowpass',false,'apply_epoch',true,'epoch_window',[-500 1000], ...
    'apply_car',false,'apply_baseline',true,'baseline_period',[-500 -50],'verbose',false);
[~, tc.TestData.E] = evalc('ieeglab_preprocess(EEG, tc.TestData.opt)');
end

% ---------------------------------------------------------------- statistics

function test_n1_permutation_is_calibrated_on_noise(tc)
% On pure noise almost nothing may be significant (the old rule took the max
% over the window without accounting for it). Enough pairs that a rate of 0.10
% cannot pass: the limit is 3 standard errors above 0.05.
E = tc.TestData.E;
p = []; sig = [];
for seed = 1:30
    rng(seed);
    E.data = single(50 * randn(size(E.data)));
    [~, A] = evalc('ieeglab_detect_n1(E, struct(''n_perm'',500,''verbose'',false,''seed'',seed))');
    p = [p; A.ieeglab.n1.table.p]; %#ok<AGROW>
    sig = [sig; A.ieeglab.n1.table.significant]; %#ok<AGROW>
end
tc.assertGreaterThan(numel(p), 900);
tc.verifyLessThan(mean(sig), 0.05, sprintf('%d of %d noise pairs called significant.', sum(sig), numel(p)));
tc.verifyLessThan(mean(p < 0.05), 0.05 + 3 * sqrt(0.05 * 0.95 / numel(p)), ...
    sprintf('Uncorrected p < 0.05 on %.1f%% of noise pairs.', 100 * mean(p < 0.05)));
end

function test_crp_p_accounts_for_duration_selection(tc)
E = tc.TestData.E;
p = []; sig = [];
for seed = 1:30
    rng(seed);
    E.data = single(50 * randn(size(E.data)));
    o = struct('run_n1',false,'run_crp',true,'run_matrix',false,'plot',false,'verbose',false, ...
        'n_perm',200,'crp_window',[15 200],'seed',seed);
    [~, A] = evalc('ieeglab_stats_subject(E, o)');
    p = [p; A.ieeglab.stats.table.p]; %#ok<AGROW>
    sig = [sig; A.ieeglab.stats.table.significant]; %#ok<AGROW>
end
tc.assertGreaterThan(numel(p), 900);
tc.verifyLessThan(mean(p < 0.05), 0.05 + 3 * sqrt(0.05 * 0.95 / numel(p)), ...
    sprintf('Selection-corrected CRP p < 0.05 on %.1f%% of noise pairs.', 100 * mean(p < 0.05)));
tc.verifyLessThan(mean(sig), 0.05);
end

function test_crp_detects_real_responses(tc)
E = tc.TestData.E;
o = struct('run_n1',false,'run_crp',true,'run_matrix',false,'plot',false,'verbose',false,'n_perm',200);
[~, E] = evalc('ieeglab_stats_subject(E, o)');
T = E.ieeglab.stats.table;
tc.verifyGreaterThan(sum(T.significant), 10, 'Real CCEPs must still be detected.');
tc.verifyTrue(all(T.testable), 'BH and 200 permutations leave these pairs testable.');
end

function test_too_few_trials_count_as_not_tested(tc)
% With 5 trials a sign-flip p cannot go below 2/2^5 = 0.0625: such a pair has
% no result, and "no response" would read as a finding.
E = tc.TestData.E;
sites = ieeglab_epoch_sites(E);
keep = [];
for s = unique(sites)
    k = find(sites == s, 5);
    keep = [keep k]; %#ok<AGROW>
end
[~, E] = evalc('pop_select(E, ''trial'', keep)');
o = struct('run_n1',false,'run_crp',true,'run_matrix',true,'matrix_source','crp','plot',false,'verbose',false);
A = tc.verifyWarning(@() ieeglab_stats_subject(E, o), 'ieeglab_stats_subject:untestable');
T = A.ieeglab.stats.table; M = A.ieeglab.ccep_matrix;
tc.verifyEqual(unique(T.p_min), 0.0625);
tc.verifyFalse(any(T.testable));
tc.verifyFalse(any(T.significant));
tc.verifyEqual(M.n_tested, 0);
tc.verifyEqual(M.n_untestable, height(T));
tc.verifyTrue(all(isnan(M.response(:))), 'Not tested, not "no response".');
% N1 is one-sided: 1/2^5 = 0.031 can pass; both signs cannot
n = struct('polarity','abs','verbose',false);
[~, A] = evalc('ieeglab_detect_n1(E, n)');
tc.verifyEqual(unique(A.ieeglab.n1.table.p_min), 0.0625);
tc.verifyFalse(any(A.ieeglab.n1.table.testable));
n.polarity = 'negative';
[~, A] = evalc('ieeglab_detect_n1(E, n)');
tc.verifyEqual(unique(A.ieeglab.n1.table.p_min), 0.03125);
tc.verifyTrue(all(A.ieeglab.n1.table.testable));
end

function test_bonferroni_needs_enough_permutations(tc)
% 1/(n_perm + 1) x number of tests must stay below alpha. More permutations do
% not help a site whose sign patterns are all used: 2/2^10 x 36 tests > 0.05.
E = tc.TestData.E;
o = struct('run_n1',false,'run_crp',true,'run_matrix',false,'plot',false,'verbose',false, ...
    'n_perm',200,'correct','bonferroni','crp_window',[15 200]);
A = tc.verifyWarning(@() ieeglab_stats_subject(E, o), 'ieeglab_stats_subject:untestable');
T = A.ieeglab.stats.table;
tc.verifyGreaterThan(1 / 201 * height(T), 0.05, 'Test data: more tests than 200 permutations can correct.');
tc.verifyFalse(any(T.testable));
o.n_perm = 2000;
[~, A] = evalc('ieeglab_stats_subject(E, o)');
T = A.ieeglab.stats.table;
tc.assertTrue(any(2.^T.n_trials <= 2001), 'Test data: a site with all sign patterns used.');
tc.verifyEqual(T.testable, 2.^T.n_trials > 2001);
end

function test_no_statistics_toolbox_needed(tc)
% Every function of the toolbox the analyses could reach fails when called
d = tempname; mkdir(d);
tc.addTeardown(@rmdir, d, 's');
for f = {'ttest','tcdf','prctile','quantile','corr','datasample'}
    fid = fopen(fullfile(d, [f{1} '.m']), 'w');
    fprintf(fid, 'function varargout = %s(varargin)\nerror(''test:toolbox'', ''needs the toolbox'');\nend\n', f{1});
    fclose(fid);
end
addpath(d);
tc.addTeardown(@rmpath, d);
E = tc.TestData.E;
o = struct('run_n1',false,'run_crp',true,'run_matrix',false,'plot',false,'verbose',false,'n_perm',200);
for t = {'permutation','ttest'}
    o.crp_test = t{1};
    A = tc.verifyWarningFree(@() ieeglab_stats_subject(E, o));
    tc.verifyTrue(all(isfinite(A.ieeglab.stats.table.p)));
end
for m = {'carla','varsubset'}
    [~, ~, out] = evalc('ieeglab_car(E, struct(''car_method'', m{1}, ''car_nboot'', 10, ''verbose'', false))');
    tc.verifyNotEmpty(out);
end
end

function test_crp_t_test_p_matches_ttest(tc)
tc.assumeNotEmpty(which('ttest'), 'Statistics Toolbox not installed.');
rng(5);
t = (1:200) / 512;
V = randn(200, 12) + 0.3 * sin(2 * pi * 4 * t(:));
[~, c] = run_CRP(V, t, struct('verbose', false));
[~, p] = ttest(c.S_all(c.stat_indices, c.tR_index), 0, 'tail', 'right');
tc.verifyEqual(c.p_value_tR, p, 'RelTol', 1e-9);
[~, p] = ttest(-c.S_all(c.stat_indices, c.tR_index), 0, 'tail', 'right');   % negative t
tc.verifyEqual(1 - c.p_value_tR, p, 'RelTol', 1e-6);
end

function test_quantile_matches_matlab(tc)
rng(6);
for n = [1 2 3 101]
    x = randn(1, n);
    for q = [0 0.05 0.2 0.5 0.95 1]
        tc.verifyEqual(ieeglab_quantile(x, q), quantile(x, q), 'AbsTol', 1e-12);
        tc.verifyEqual(ieeglab_quantile(x, q), prctile(x, 100 * q), 'AbsTol', 1e-12);
    end
end
tc.verifyEqual(ieeglab_quantile([1 NaN 3], 0.5), 2, 'NaN values are ignored.');
end

% ---------------------------------------------------------------- preprocessing

function test_marked_bad_contact_trials_always_dropped(tc)
% A contact marked bad at load (as a channels.tsv status would) must lose its
% stimulation trials with default options, not only when some unrelated
% bad-channel option is on.
C = tc.TestData.C;
k = find(strcmpi({C.chanlocs.labels}, 'RA5'), 1);
tc.assumeNotEmpty(k);
for c = 1:numel(C.chanlocs), C.chanlocs(c).status = 'good'; end
C.chanlocs(k).status = 'bad'; C.chanlocs(k).status_description = 'noisy (test)';
tc.assumeTrue(any(ieeglab_events_using(C.event, {'RA5'})), 'No RA5 stimulation in the tutorial.');
[~, E] = evalc('ieeglab_preprocess(C, tc.TestData.opt)');
tc.verifyTrue(any(strcmpi({E.chanlocs.labels}, 'RA5')), 'A marked channel is kept (marked), not removed.');
[~, stimIdx] = ieeglab_epoch_sites(E);
badIdx = find(strcmpi({E.chanlocs.labels}, 'RA5'));
tc.verifyFalse(any(cellfun(@(s) any(ismember(s, badIdx)), stimIdx)), ...
    'A trial stimulating a marked-bad contact survived default preprocessing.');
o = tc.TestData.opt; o.exclude_soz = true;       % no SOZ labels: must change nothing
[~, E2] = evalc('ieeglab_preprocess(C, o)');
tc.verifyEqual(E2.trials, E.trials, 'exclude_soz with no SOZ contact changed the trials.');
end

function test_char_event_filter_equals_cell(tc)
C = tc.TestData.C; o = tc.TestData.opt;
o.event_filters = struct('type', 'RA3-RA4');
[~, A] = evalc('ieeglab_preprocess(C, o)');
o.event_filters = struct('type', {{'RA3-RA4'}});
[~, B] = evalc('ieeglab_preprocess(C, o)');
tc.verifyEqual(A.trials, B.trials);
tc.verifyGreaterThan(A.trials, 10);
tc.verifyTrue(all(ieeglab_epoch_sites(A) == "RA3-RA4"));
end

function test_filter_removing_everything_errors(tc)
C = tc.TestData.C; o = tc.TestData.opt;
o.event_filters = struct('type', {{'NOPE'}});
tc.verifyError(@() ieeglab_preprocess(C, o), 'ieeglab_preprocess:allEventsFiltered');
end

function test_rare_conditions_counted_per_site(tc)
% A site recorded with both polarities must be counted as one site. Record
% half of the RA4-RA5 pulses as RA5-RA4.
C = tc.TestData.C;
k = find(strcmp({C.event.type}, 'RA4-RA5'));
[C.event(k(1:2:end)).type] = deal('RA5-RA4');
types = string({C.event.type});
nSite = nnz(types == "RA4-RA5" | types == "RA5-RA4");
nMax  = max(nnz(types == "RA4-RA5"), nnz(types == "RA5-RA4"));
o = tc.TestData.opt; o.remove_rare_cond = true; o.min_trials = nMax + 1;
[~, E] = evalc('ieeglab_preprocess(C, o)');
tc.verifyEqual(nnz(ieeglab_epoch_sites(E) == "RA4-RA5"), nSite, ...
    'A site with enough trials in total was removed because each polarity alone was rare.');
end

function test_boundary_events_survive_selection(tc)
C = tc.TestData.C; o = tc.TestData.opt;
lat = C.event(find(strcmp({C.event.type}, 'RA3-RA4'), 1, 'last')).latency;
B = C;
B.event(end+1).type = 'boundary';
B.event(end).latency = lat + round(0.1 * C.srate);
B.event(end).duration = 0;
B = eeg_checkset(B, 'eventconsistency');
o.remove_rare_cond = true; o.min_trials = 5;
o.event_filters = struct('type', {{'RA3-RA4'}});
[~, E0] = evalc('ieeglab_preprocess(C, o)');
[~, E1] = evalc('ieeglab_preprocess(B, o)');
tc.verifyEqual(E1.trials, E0.trials - 1, 'The epoch spanning the boundary must be rejected.');
end

function test_numeric_bad_channels_resolved_before_removal(tc)
C = tc.TestData.C;
lab = {C.chanlocs.labels};
i2 = 2;
C.chanlocs(i2).X = []; C.chanlocs(i2).Y = []; C.chanlocs(i2).Z = [];
o = tc.TestData.opt;
o.remove_no_coords = true; o.remove_bad_channels = true;
o.bad_channels = 3;                                % channel 3 OF THE INPUT
[~, E] = evalc('ieeglab_preprocess(C, o)');
out = {E.chanlocs.labels};
tc.verifyFalse(any(strcmpi(out, lab{3})), sprintf('%s (index 3 of the input) should be removed.', lab{3}));
tc.verifyTrue(any(strcmpi(out, lab{4})), sprintf('%s must not be removed in place of %s.', lab{4}, lab{3}));
tc.verifyFalse(any(strcmpi(out, lab{2})), 'The channel without coordinates should be removed.');
end

function test_exclude_soz_alone_removes_zone(tc)
C = tc.TestData.C;
k = find(strcmpi({C.chanlocs.labels}, 'RA4'), 1);
tc.assumeNotEmpty(k);
C.chanlocs(k).clinical_zone = 'SOZ';
o = tc.TestData.opt; o.exclude_soz = true;
[~, E] = evalc('ieeglab_preprocess(C, o)');
tc.verifyFalse(any(strcmpi({E.chanlocs.labels}, 'RA4')), 'exclude_soz must remove the SOZ contact.');
marked = {C.chanlocs(strcmpi({C.chanlocs.status}, 'bad')).labels};
tc.verifyTrue(all(ismember(marked, {E.chanlocs.labels})), ...
    'Clinician-marked channels stay (marked) unless remove_bad_channels is on.');
end

function test_stored_legacy_aliases_do_not_override(tc)
E = tc.TestData.E;
E.ieeglab.opt.apply_acar = false;
E.ieeglab.opt.acar_timewin = [15 500];
[~, E2] = evalc('ieeglab_preprocess(E, struct(''apply_car'',true,''car_method'',''car'',''verbose'',false))');
tc.verifyEqual(E2.ref, 'CAR', 'An explicit apply_car=true must re-reference.');
end

function test_stored_steps_are_not_reapplied(tc)
% Re-running only the baseline must not high-pass, re-epoch or re-reference.
E = tc.TestData.E;
o = struct('apply_baseline', true, 'baseline_period', [-400 -50], 'verbose', false);
[~, E2] = evalc('ieeglab_preprocess(E, o)');
tc.verifyEqual(E2.trials, E.trials);
D = double(E2.data) - double(E.data);
% Single-precision data: allow rounding at the data's own scale (a second
% high-pass or re-reference would change the waveform by microvolts).
tol = 10 * double(eps(single(max(abs(E.data(:))))));
tc.verifyLessThan(max(std(D, 0, 2), [], 'all'), tol, 'Only a per-trial constant may change.');
tc.verifyFalse(isfield(E2.ieeglab.opt, 'plot'), 'The plot switch must not be stored.');
end

function test_preprocess_history_replays(tc)
C = tc.TestData.C;
[~, E, com] = evalc('ieeglab_preprocess(C, tc.TestData.opt)');
tc.verifySubstring(com, 'ieeglab_preprocess(EEG, struct(');
EEG = C; %#ok<NASGU>
evalc(com);
tc.verifyEqual(EEG.trials, E.trials);
tc.verifyEqual(double(EEG.data), double(E.data), 'AbsTol', 1e-6);
end

function test_blank_nan_survives_filtering(tc)
C = tc.TestData.C; o = tc.TestData.opt;
o.apply_blank = true; o.blank_method = 'nan'; o.blank_window = [-1.5 5];
[~, E] = evalc('ieeglab_preprocess(C, o)');
f = mean(isnan(E.data(:)));
tc.verifyGreaterThan(f, 0);
tc.verifyLessThan(f, 0.05, sprintf('%.1f%% of samples are NaN; filtering spread the blanks.', 100*f));
end

function test_blank_on_epoched_data(tc)
E = tc.TestData.E;
[~, ~, info] = evalc('ieeglab_blank_stim(E, struct(''verbose'',false))');
tc.verifyGreaterThanOrEqual(info.n_blanked, E.trials, 'Every epoch must be blanked, not just the first.');
end

% ---------------------------------------------------------------- re-referencing

function test_varsubset_matches_reference(tc)
% apply_ieeg_car (tests/reference) is the HAPwave implementation of the legacy
% lowest-variance subset CAR that car_method 'varsubset' reproduces.
root = tc.TestData.root;
addpath(fullfile(root, 'tests', 'reference'));
c = onCleanup(@() rmpath(fullfile(root, 'tests', 'reference')));
E = tc.TestData.E;
% Non-CCEP labelling, so there is one group and no contact exclusion
for k = 1:numel(E.event), E.event(k).type = 'stim'; end
for k = 1:numel(E.epoch), E.epoch(k).eventtype = repmat({'stim'}, size(E.epoch(k).eventtype)); end
for k = 1:numel(E.chanlocs), E.chanlocs(k).status = 'good'; end
frac = 0.3; win = [15 500];
[~, A] = evalc('ieeglab_car(E, struct(''car_method'',''varsubset'',''car_fraction'',frac,''car_timewin'',win,''verbose'',false))');
[~, ref] = evalc('apply_ieeg_car(double(E.data), double(E.times)/1000, 1:E.nbchan, frac, win/1000, false)');
tc.verifyEqual(double(A.data), ref, 'AbsTol', 1e-6);
tc.verifyEqual(A.ieeglab.car.n_groups, 1, 'Non-CCEP data must form a single reference group.');
end

function test_carla_ignores_all_nan_sample(tc)
E = tc.TestData.E;
sites = ieeglab_epoch_sites(E);
tr = find(sites == "RA3-RA4");
[~, stimIdx] = ieeglab_epoch_sites(E);
V = double(E.data(:, :, tr));
V(unique(vertcat(stimIdx{tr})), :, :) = NaN;
t = double(E.times) / 1000;
o = struct('winResp', [0.01 0.3], 'nboot', 20, 'verbose', false);
rng(3); [~, ~, s0] = ieeglab_carla(t, V, E.srate, o);
k = find(E.times >= 10, 1);
V(:, k, :) = NaN;
rng(3); [~, ~, s1] = ieeglab_carla(t, V, E.srate, o);
tc.verifyFalse(isfield(s1,'fallback') && s1.fallback);
tc.verifyGreaterThan(numel(s1.chsUsed), 2, 'CARLA collapsed to the 2-channel floor.');
tc.verifyGreaterThanOrEqual(numel(intersect(s0.chsUsed, s1.chsUsed)), 0.7 * numel(s0.chsUsed));
end

% ---------------------------------------------------------------- matrix, plots, export

function test_matrix_degrees_and_staleness(tc)
E = tc.TestData.E;
[~, E] = evalc('ieeglab_ccep_matrix(E, struct(''verbose'',false,''n_perm'',200))');
M = E.ieeglab.ccep_matrix;
never = all(isnan(M.response), 1)';
tc.verifyTrue(all(isnan(M.in_degree(never))), 'Never-tested contacts must have in_degree NaN, not 0.');
tc.verifyTrue(ieeglab_matrix_current(E));
[~, E2] = evalc('ieeglab_detect_n1(E, struct(''n_perm'',200,''n1_window'',[40 100],''verbose'',false))');
tc.verifyFalse(isfield(E2.ieeglab, 'ccep_matrix') && ieeglab_matrix_current(E2), ...
    'Recomputing N1 must invalidate the matrix built from the old N1.');
end

function test_topoplot_masks_nonsignificant_pairs(tc)
E = tc.TestData.E;
[~, E] = evalc('ieeglab_ccep_matrix(E, struct(''verbose'',false,''n_perm'',200))');
M = E.ieeglab.ccep_matrix;
[~, r] = max(sum(M.response == 0, 2));
[v, info] = ieeglab_topoplot(E, 'n1_latency', 'site', M.sites{r}, 'draw', false); %#ok<ASGLU>
[tf, loc] = ismember(upper({E.chanlocs.labels}), upper(M.channels));
resp = nan(numel(v),1); resp(tf) = M.response(r, loc(tf));
tc.verifyTrue(all(isnan(v(resp == 0))), 'Non-significant pairs must not be drawn as responses.');
tc.verifyTrue(all(isfinite(v(resp == 1))));
v2 = ieeglab_topoplot(E, 'n1_latency', 'site', M.sites{r}, 'draw', false, 'sig_only', false);
tc.verifyGreaterThan(nnz(isfinite(v2)), nnz(isfinite(v)));
end

function test_topoplot_latency_excludes_own_stimulation(tc)
E = tc.TestData.E;
sites = ieeglab_epoch_sites(E);
s = sites(1);
[vAll, iAll] = ieeglab_topoplot(E, 20, 'draw', false); %#ok<ASGLU>
[vSite, iSite] = ieeglab_topoplot(E, 20, 'site', char(s), 'draw', false);
tc.verifyNotEqual(vAll, vSite, '''site'' must restrict the average to that site''s trials.');
[~, stimIdx] = ieeglab_epoch_sites(E);
st = stimIdx{find(sites == s, 1)};
tc.verifyTrue(all(isnan(vSite(st))), 'The stimulated contacts have no trial left once their own are excluded.');
tc.verifyTrue(contains(strjoin(cellstr(iSite.title), ' '), char(s)));
tc.verifyError(@() ieeglab_topoplot(E, [900 5000], 'draw', false), 'ieeglab_topoplot:latencyOutside');
end

function test_export_without_status_fields_and_overwrite(tc)
E = tc.TestData.E;
E.chanlocs = rmfield(E.chanlocs, intersect(fieldnames(E.chanlocs), {'status','status_description','clinical_zone'}));
d = tempname; c = onCleanup(@() rmdir(d, 's'));
files = ieeglab_export(E, d, struct('formats', {{'tsv'}}, 'verbose', false));
tc.verifyNumElements(files, 1);
T = readtable(files{1}, 'FileType','text', 'Delimiter','\t', 'TextType','string');
tc.verifyTrue(all(T.status == "n/a"));
tc.verifyError(@() ieeglab_export(E, d, struct('overwrite', false, 'verbose', false)), 'ieeglab_export:exists');
tc.verifyEmpty(dir(fullfile(d, '*.json')), 'A refused export must write nothing.');
tc.verifyError(@() ieeglab_export(E, d, struct('formats', {{'csv'}}, 'verbose', false)), 'ieeglab_export:badFormat');
end

function test_export_history_replays(tc)
E = tc.TestData.E;
d = fullfile(tempname, 'it''s here'); c = onCleanup(@() rmdir(fileparts(d), 's'));
[files, com] = ieeglab_export(E, d, struct('formats', {{'tsv','json'}}, 'prefix', 'custom', 'verbose', false));
cellfun(@delete, files);
EEG = E; %#ok<NASGU>
evalc(com);
tc.verifyEqual(numel(dir(fullfile(d, 'custom_*'))), numel(files));
end

function test_literal_round_trip(tc)
o = struct('a', 'o''brien', 'b', {{'x','y'}}, 'c', [1 2; 3 4], 'd', true, ...
    'e', struct('type', {{'RA1-RA2'}}), 'f', "str", 'g', ["s1" "s2"]);
o2 = eval(ieeglab_literal(o));
tc.verifyEqual(o2.a, o.a); tc.verifyEqual(o2.b, o.b); tc.verifyEqual(o2.c, o.c);
tc.verifyEqual(o2.d, o.d); tc.verifyEqual(o2.e.type, {'RA1-RA2'});
tc.verifyEqual(o2.f, 'str'); tc.verifyEqual(o2.g, {'s1','s2'});
end
