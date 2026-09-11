function tests = test_ieeglab_audit_regressions
% Regression tests for the findings of the adversarial audit (v1.1.0).
% Headless and compute-only: nothing here opens a figure or a dialog.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(mfilename('fullpath')));
tc.assumeNotEmpty(which('eeglab'), 'EEGLAB is not on the MATLAB path.');
if isempty(which('pop_loadset')), evalc('eeglab nogui'); end
addpath(root); addpath(fullfile(root, 'functions'));
tc.TestData.root = root;
d = fullfile(root, 'tutorial', 'dataset_seeg');
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
% over the window without accounting for it).
E = tc.TestData.E;
rng(1);
E.data = single(50 * randn(size(E.data)));
[~, E] = evalc('ieeglab_detect_n1(E, struct(''n_perm'',500,''verbose'',false))');
T = E.ieeglab.n1.table;
tc.assertGreaterThan(height(T), 50);
tc.verifyLessThan(mean(T.significant), 0.05, ...
    sprintf('%d of %d noise pairs called significant.', sum(T.significant), height(T)));
tc.verifyLessThan(mean(T.p < 0.05), 0.10, 'Uncorrected p-values on noise must be roughly uniform.');
end

function test_crp_p_accounts_for_duration_selection(tc)
E = tc.TestData.E;
rng(2);
E.data = single(50 * randn(size(E.data)));
o = struct('run_n1',false,'run_crp',true,'run_matrix',false,'plot',false,'verbose',false, ...
    'n_perm',200,'crp_window',[15 200]);
[~, E] = evalc('ieeglab_stats_subject(E, o)');
T = E.ieeglab.stats.table;
tc.assertGreaterThan(height(T), 50);
tc.verifyLessThan(mean(T.p < 0.05), 0.12, ...
    sprintf('Selection-corrected CRP p < 0.05 on %.0f%% of noise pairs.', 100*mean(T.p < 0.05)));
tc.verifyLessThan(mean(T.significant), 0.05);
end

function test_crp_detects_real_responses(tc)
E = tc.TestData.E;
o = struct('run_n1',false,'run_crp',true,'run_matrix',false,'plot',false,'verbose',false,'n_perm',200);
[~, E] = evalc('ieeglab_stats_subject(E, o)');
tc.verifyGreaterThan(sum(E.ieeglab.stats.table.significant), 10, 'Real CCEPs must still be detected.');
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
o.event_filters = struct('type', 'ROP1-ROP2');
[~, A] = evalc('ieeglab_preprocess(C, o)');
o.event_filters = struct('type', {{'ROP1-ROP2'}});
[~, B] = evalc('ieeglab_preprocess(C, o)');
tc.verifyEqual(A.trials, B.trials);
tc.verifyGreaterThan(A.trials, 10);
tc.verifyTrue(all(ieeglab_epoch_sites(A) == "ROP1-ROP2"));
end

function test_filter_removing_everything_errors(tc)
C = tc.TestData.C; o = tc.TestData.opt;
o.event_filters = struct('type', {{'NOPE'}});
tc.verifyError(@() ieeglab_preprocess(C, o), 'ieeglab_preprocess:allEventsFiltered');
end

function test_rare_conditions_counted_per_site(tc)
% ROP2-ROP4 is recorded with both polarities; it must be counted as one site.
C = tc.TestData.C;
types = string({C.event.type});
tc.assumeTrue(any(types == "ROP2-ROP4") && any(types == "ROP4-ROP2"), 'Tutorial lacks the mixed-polarity site.');
nSite = nnz(types == "ROP2-ROP4" | types == "ROP4-ROP2");
nMax  = max(nnz(types == "ROP2-ROP4"), nnz(types == "ROP4-ROP2"));
o = tc.TestData.opt; o.remove_rare_cond = true; o.min_trials = nMax + 1;
tc.assumeLessThan(nMax + 1, nSite + 1);
[~, E] = evalc('ieeglab_preprocess(C, o)');
tc.verifyGreaterThan(nnz(ieeglab_epoch_sites(E) == "ROP2-ROP4"), nMax, ...
    'A site with enough trials in total was removed because each polarity alone was rare.');
end

function test_boundary_events_survive_selection(tc)
C = tc.TestData.C; o = tc.TestData.opt;
lat = C.event(find(strcmp({C.event.type}, 'ROP1-ROP2'), 1, 'last')).latency;
B = C;
B.event(end+1).type = 'boundary';
B.event(end).latency = lat + round(0.1 * C.srate);
B.event(end).duration = 0;
B = eeg_checkset(B, 'eventconsistency');
o.remove_rare_cond = true; o.min_trials = 5;
o.event_filters = struct('type', {{'ROP1-ROP2'}});
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

function test_varsubset_matches_archive(tc)
root = tc.TestData.root;
addpath(fullfile(root, 'archive'));
c = onCleanup(@() rmpath(fullfile(root, 'archive')));
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
tr = find(sites == "ROP1-ROP2");
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
