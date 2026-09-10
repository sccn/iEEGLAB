function tests = test_ieeglab_features
% test_ieeglab_features - Tests for bad channels, the headless loader, the
% connectivity matrix, export and electrode-value plotting.
%
% All headless. Plot functions are exercised through their 'draw', false
% path, which returns exactly what would be drawn: on some machines MATLAB's
% -batch mode hangs on ANY graphics object, so rendering itself is not tested
% here.
%
% Cedric Cannard, iEEGLAB, 2026

tests = functiontests(localfunctions);
end

% ======================= fixtures =======================

function setupOnce(tc)
root = fileparts(fileparts(mfilename('fullpath')));
tc.TestData.root = root;
tc.assumeNotEmpty(which('eeglab'), 'EEGLAB is not on the MATLAB path.');
if isempty(which('pop_loadset')), evalc('eeglab nogui'); end
addpath(root); addpath(fullfile(root, 'functions'));
set(0, 'DefaultFigureVisible', 'off');
tc.TestData.seegDir = fullfile(root, 'tutorial', 'dataset_seeg');
tc.TestData.set = 'sub-02_ses-ieeg01_task-ccep_run-01_ieeg.set';
tc.TestData.tmp = tempname; mkdir(tc.TestData.tmp);

% One loaded + preprocessed + analysed dataset shared by the read-only tests
EEG = local_load_headless(tc, struct());
opt = struct('apply_highpass',true,'highpass',0.5,'apply_notch',false,'apply_lowpass',false, ...
    'apply_ds',false,'apply_epoch',true,'epoch_window',[-500 1000], ...
    'apply_car',true,'car_method','carla','car_nboot',15, ...
    'apply_baseline',true,'baseline_period',[-500 -50],'plot',false,'verbose',false);
[~, E] = evalc('ieeglab_preprocess(EEG, opt)');
[~, E] = evalc(['ieeglab_stats_subject(E, struct(''run_n1'',true,''run_crp'',true,' ...
    '''run_matrix'',true,''min_trials'',8,''crp_window'',[15 400],''n1_window'',[15 100],' ...
    '''verbose'',false))']);
tc.TestData.E = E;
end

function teardownOnce(tc)
if isfolder(tc.TestData.tmp), rmdir(tc.TestData.tmp, 's'); end
end

function EEG = local_load_headless(tc, opt)
EEG = pop_loadset('filename', tc.TestData.set, 'filepath', tc.TestData.seegDir);
[~, EEG] = evalc('ieeglab_load(EEG, opt)');
end

function f = local_channels_tsv(tc, labels, bad, badWhy, ecg)
% Write a BIDS channels.tsv into the temp folder.
f = fullfile(tc.TestData.tmp, sprintf('channels_%d.tsv', randi(1e9)));
fid = fopen(f, 'w');
fprintf(fid, 'name\ttype\tunits\tstatus\tstatus_description\n');
for i = 1:numel(labels)
    ty = 'SEEG'; st = 'good'; why = 'n/a';
    k = find(strcmp(bad, labels{i}), 1);
    if ~isempty(k), st = 'bad'; why = badWhy{k}; end
    if any(strcmp(ecg, labels{i})), ty = 'ECG'; end
    fprintf(fid, '%s\t%s\tuV\t%s\t%s\n', labels{i}, ty, st, why);
end
fprintf(fid, 'NOTINDATA\tSEEG\tuV\tbad\tn/a\n');     % must be ignored gracefully
fclose(fid);
end

% ======================= shared parsing =======================

function test_site_tokens(tc)
tc.verifyEqual(ieeglab_site_tokens('RA1-RA2'), ["RA1" "RA2"], ...
    'char input used to be split into single characters (R A 1 R A 2).');
tc.verifyEqual(ieeglab_site_tokens("RA1-RA2"), ["RA1" "RA2"]);
tc.verifyEqual(ieeglab_site_tokens({'RA1-RA2'}), ["RA1" "RA2"]);
tc.verifyEqual(ieeglab_site_tokens("LA'3 / LA'4"), ["LA'3" "LA'4"], 'Apostrophes in sEEG labels must survive.');
tc.verifyEmpty(ieeglab_site_tokens('n/a'), 'BIDS n/a used to split into the tokens "n" and "a".');
tc.verifyEmpty(ieeglab_site_tokens(''));
end

function test_events_using_is_exact_not_substring(tc)
% Regression: contains(types, 'RA1') also matched 'RA10-RA9'.
types = {'RA1-RA2', 'RA10-RA9', 'RA2-RA3', 'face'};
m = ieeglab_events_using(types, {'RA1'});
tc.verifyEqual(m, [true false false false], ...
    'Removing contact RA1 must not drop RA10-RA9 trials.');
tc.verifyFalse(any(ieeglab_events_using({'RA1'}, {'RA1'})), ...
    'A single-token event is a condition name, not a stimulation pair.');
end

function test_epoch_sites_are_order_independent(tc)
E = tc.TestData.E;
[sites, stimIdx] = ieeglab_epoch_sites(E);
tc.verifyNumElements(sites, E.trials);
tc.verifyFalse(any(sites == "ROP4-ROP2"), 'ROP4-ROP2 must be canonicalised to ROP2-ROP4.');
tc.verifyTrue(any(sites == "ROP2-ROP4"));
tc.verifyTrue(all(cellfun(@numel, stimIdx) == 2), 'Every CCEP epoch should resolve both stimulated contacts.');
end

% ======================= headless loader =======================

function test_load_headless_finds_bids_siblings(tc)
EEG = local_load_headless(tc, struct());
tc.verifyEqual(numel(EEG.event), 169, 'All 169 stimulation events should be placed.');
tc.verifyTrue(all(arrayfun(@(c) ~isempty(c.X) && isfinite(c.X), EEG.chanlocs)), 'Every channel should get coordinates.');
tc.verifyEqual(EEG.ieeglab.opt.event_field, 'electrical_stimulation_site');
tc.verifyEqual(ieeglab_detect_mode(EEG), 'ccep');
ev = readtable(fullfile(tc.TestData.seegDir, 'sub-02_ses-ieeg01_task-ccep_run-01_events.tsv'), ...
    'FileType','text', 'Delimiter','	');
tc.verifyEqual(EEG.event(1).latency, min(ev.onset) * EEG.srate + 1, 'AbsTol', 1e-9, ...
    'Latencies must be 1-based samples: onset 0 s is sample 1.');
end

function test_load_history_is_rerunnable(tc)
EEG = pop_loadset('filename', tc.TestData.set, 'filepath', tc.TestData.seegDir);
[~, E1, com] = evalc('ieeglab_load(EEG, struct())');
tc.verifySubstring(com, 'ieeglab_load(EEG, struct(');
evalc(com);             % the history line is an assignment to EEG; run it as written
E2 = EEG;
tc.verifyEqual(numel(E2.event), numel(E1.event));
tc.verifyEqual({E2.chanlocs.labels}, {E1.chanlocs.labels});
end

function test_load_drops_events_marked_bad(tc)
src = fullfile(tc.TestData.seegDir, 'sub-02_ses-ieeg01_task-ccep_run-01_events.tsv');
T = readtable(src, 'FileType','text', 'Delimiter','\t');
T.status = cellstr(string(T.status));
T.status(1:3) = {'bad'};
f = fullfile(tc.TestData.tmp, 'events_with_bad.tsv');
writetable(T, f, 'FileType','text', 'Delimiter','\t');
EEG = local_load_headless(tc, struct('events_tsv', f));
tc.verifyEqual(numel(EEG.event), 166, 'Events with status=bad must be dropped at load.');
EEG = local_load_headless(tc, struct('events_tsv', f, 'drop_bad_events', false));
tc.verifyEqual(numel(EEG.event), 169);
end

function test_load_channel_selection_is_exact(tc)
% Removing RA1 must drop RA1-RA2 trials but keep RA10-RA9 trials.
EEG0 = local_load_headless(tc, struct());
keep = setdiff({EEG0.chanlocs.labels}, {'RA1'}, 'stable');
EEG = local_load_headless(tc, struct('chan_list', {keep}));
types = {EEG.event.type};
tc.verifyEqual(EEG.nbchan, EEG0.nbchan - 1);
tc.verifyFalse(any(strcmp(types, 'RA1-RA2')), 'Trials stimulating the removed contact must go.');
tc.verifyEqual(nnz(strcmp(types, 'RA10-RA9')), nnz(strcmp({EEG0.event.type}, 'RA10-RA9')), ...
    'RA10-RA9 trials must survive removing RA1 (substring-match regression).');
end

% ======================= bad channels =======================

function test_bad_channels_from_channels_tsv_mark(tc)
EEG = local_load_headless(tc, struct());
labels = {EEG.chanlocs.labels};
f = local_channels_tsv(tc, labels, {'RA3','ROP5'}, {'high impedance','n/a'}, {'ROP6'});
[~, E, T] = evalc('ieeglab_bad_channels(EEG, struct(''channels_tsv'', f, ''action'', ''mark''))');
tc.verifyEqual(E.nbchan, EEG.nbchan, 'mark must not remove channels.');
bad = cellstr(T.label(T.status == "bad"));
tc.verifyEqual(sort(bad(:))', sort({'RA3','ROP5','ROP6'}));
tc.verifyEqual(char(T.reason(T.label == "RA3")), 'high impedance', 'status_description must be kept as the reason.');
tc.verifySubstring(char(T.reason(T.label == "ROP6")), 'ECG', 'Non-iEEG types must be flagged.');
tc.verifyEqual(E.chanlocs(strcmp(labels,'RA3')).status, 'bad');
end

function test_bad_channels_remove_drops_their_stim_trials_exactly(tc)
EEG = local_load_headless(tc, struct());
labels = {EEG.chanlocs.labels};
f = local_channels_tsv(tc, labels, {'RA1'}, {'noisy'}, {});
nRA10 = nnz(strcmp({EEG.event.type}, 'RA10-RA9'));
[~, E] = evalc('ieeglab_bad_channels(EEG, struct(''channels_tsv'', f, ''action'', ''remove''))');
tc.verifyFalse(any(strcmp({E.chanlocs.labels}, 'RA1')));
tc.verifyFalse(any(strcmp({E.event.type}, 'RA1-RA2')), 'Stimulation of a bad contact must be dropped.');
tc.verifyEqual(nnz(strcmp({E.event.type}, 'RA10-RA9')), nRA10, 'RA10-RA9 must survive (exact matching).');
tc.verifyEqual(E.ieeglab.removed_channels, {'RA1'});
end

function test_marked_bad_channel_is_kept_out_of_the_reference(tc)
E = tc.TestData.E;
k = find(strcmp({E.chanlocs.labels}, 'RA5'));
E.chanlocs(k).status = 'bad';
[~, out] = ieeglab_car(E, struct('car_method','carla','car_nboot',10,'verbose',false));
for g = 1:numel(out)
    tc.verifyFalse(ismember(k, out(g).car_channels), ...
        sprintf('Site %s: a channel marked bad entered the reference.', out(g).group));
end
end

function test_mark_then_remove_honours_previous_marks(tc)
EEG = local_load_headless(tc, struct());
labels = {EEG.chanlocs.labels};
f = local_channels_tsv(tc, labels, {'ROP3'}, {'flat'}, {});
[~, E] = evalc('ieeglab_bad_channels(EEG, struct(''channels_tsv'', f, ''action'', ''mark''))');
% second pass reads no file at all, as preprocessing does
[~, E2] = evalc('ieeglab_bad_channels(E, struct(''channels_tsv'', '''', ''action'', ''remove''))');
tc.verifyFalse(any(strcmp({E2.chanlocs.labels}, 'ROP3')), 'A channel marked at load must be removable later.');
end

function test_seizure_zone_is_recorded_and_optionally_excluded(tc)
% None of the 16 tutorial channels carry a seizure_zone label, so label two of
% them in a copy of the real electrodes table.
EEG = local_load_headless(tc, struct());
Te = readtable(fullfile(tc.TestData.seegDir, 'sub-02_ses-ieeg01_electrodes.tsv'), ...
    'FileType','text', 'Delimiter','	', 'TextType','string');
Te.seizure_zone(Te.name == "RA4")  = "SOZ";
Te.seizure_zone(Te.name == "ROP2") = "IrritativeZone";
[~, E] = evalc('ieeglab_bad_channels(EEG, struct(''channels_tsv'','''', ''elec_tsv'', Te, ''action'',''mark''))');
lab = {E.chanlocs.labels};
tc.verifyEqual(E.chanlocs(strcmp(lab,'RA4')).clinical_zone, 'SOZ', 'The zone label must be recorded.');
tc.verifyEqual(E.chanlocs(strcmp(lab,'RA4')).status, 'good', 'Recording a zone must not by itself mark the contact bad.');
[~, E] = evalc('ieeglab_bad_channels(EEG, struct(''channels_tsv'','''', ''elec_tsv'', Te, ''exclude_soz'', true, ''action'',''mark''))');
tc.verifyEqual(E.chanlocs(strcmp(lab,'RA4')).status, 'bad', 'exclude_soz must mark SOZ contacts bad.');
tc.verifyEqual(E.chanlocs(strcmp(lab,'ROP2')).status, 'good', 'exclude_soz alone must not touch irritative-zone contacts.');
[~, E] = evalc('ieeglab_bad_channels(EEG, struct(''channels_tsv'','''', ''elec_tsv'', Te, ''exclude_irritative'', true, ''action'',''mark''))');
tc.verifyEqual(E.chanlocs(strcmp(lab,'ROP2')).status, 'bad');
end

function test_preprocess_bad_channel_step(tc)
EEG = local_load_headless(tc, struct());
opt = struct('remove_bad_channels',true,'bad_channels',{{'RA3'}},'apply_highpass',false, ...
    'apply_notch',false,'apply_lowpass',false,'apply_ds',false,'apply_epoch',true, ...
    'epoch_window',[-500 1000],'apply_car',false,'apply_baseline',false,'plot',false,'verbose',false);
[~, E] = evalc('ieeglab_preprocess(EEG, opt)');
tc.verifyFalse(any(strcmp({E.chanlocs.labels}, 'RA3')));
sites = ieeglab_epoch_sites(E);
tc.verifyFalse(any(sites == "RA2-RA3" | sites == "RA3-RA4"), 'Trials stimulating RA3 must be gone.');
tc.verifyTrue(any(sites == "RA10-RA9"));
end

% ======================= connectivity matrix =======================

function test_ccep_matrix_structure(tc)
M = tc.TestData.E.ieeglab.ccep_matrix;
E = tc.TestData.E;
tc.verifyEqual(numel(M.channels), E.nbchan);
tc.verifySize(M.response, [numel(M.sites) numel(M.channels)]);
tc.verifyTrue(all(ismember(M.response(~isnan(M.response)), [0 1])), 'response must be 0, 1 or NaN.');
tc.verifyEqual(M.out_degree, sum(M.response == 1, 2));
tc.verifyEqual(M.in_degree, sum(M.response == 1, 1)');
tc.verifyGreaterThanOrEqual(M.density, 0); tc.verifyLessThanOrEqual(M.density, 1);
% Stimulated contacts are "not measured" for their own site, never 0 or 1
for s = 1:numel(M.sites)
    tok = upper(ieeglab_site_tokens(M.sites{s}));
    cols = ismember(upper(M.channels), tok);
    tc.verifyTrue(all(isnan(M.response(s, cols))), sprintf('%s: stimulated contacts must be NaN.', M.sites{s}));
end
% Sites follow montage order, not alphabetical: RA10-RA9 (contacts 9,10) comes
% after RA8-RA9 (8,9), and every RA site precedes every ROP site.
tc.verifyLessThan(find(strcmp(M.sites,'RA8-RA9')), find(strcmp(M.sites,'RA10-RA9')), ...
    'RA10-RA9 must follow RA8-RA9 (montage order, not string order).');
firstROP = find(startsWith(M.sites, 'ROP'), 1);
lastRA = find(startsWith(M.sites, 'RA'), 1, 'last');
tc.verifyLessThan(lastRA, firstROP, 'Sites should be ordered by montage position.');
end

function test_ccep_matrix_agrees_with_n1_table(tc)
E = tc.TestData.E; M = E.ieeglab.ccep_matrix; T = E.ieeglab.n1.table;
tc.verifyEqual(M.n_significant, nnz(T.significant));
tc.verifyEqual(M.n_tested, height(T));
r = find(T.significant, 1);
si = strcmp(M.sites, T.site{r}); ci = strcmp(M.channels, T.channel{r});
tc.verifyEqual(M.amplitude_uv(si, ci), T.n1_amplitude_uv(r), 'AbsTol', 1e-12);
end

function test_ccep_matrix_from_crp(tc)
[~, E, M] = evalc('ieeglab_ccep_matrix(tc.TestData.E, struct(''source'',''crp'',''verbose'',false))');
tc.verifyEqual(M.source, 'crp');
tc.verifyTrue(isfield(M, 'tR_ms') && isfield(M, 'explained_var'));
tc.verifyEqual(M.n_significant, nnz(E.ieeglab.stats.table.significant));
end

function test_plot_matrix_compute_only(tc)
M = tc.TestData.E.ieeglab.ccep_matrix;
[img, info] = ieeglab_plot_ccep_matrix(M, 'amplitude_uv', 'draw', false);
tc.verifySize(img, size(M.response));
tc.verifyTrue(all(isnan(img(M.response ~= 1))), 'Only significant cells are coloured by default.');
tc.verifyTrue(all(img(M.response == 1) >= 0), 'Amplitude is shown as a magnitude.');
tc.verifyEmpty(findall(0,'Type','figure'));
tc.verifyError(@() ieeglab_plot_ccep_matrix(M, 'nonsense', 'draw', false), 'ieeglab_plot_ccep_matrix:badMetric');
tc.verifySubstring(info.title, sprintf('%d of %d', M.n_significant, M.n_tested));
end

% ======================= export =======================

function test_export_writes_readable_files(tc)
d = fullfile(tc.TestData.tmp, 'export');
[~, files] = evalc('ieeglab_export(tc.TestData.E, d)');
names = cellfun(@local_name, files, 'UniformOutput', false);
for want = {'desc-channels','desc-n1','desc-crp','desc-ccepresponse','desc-ccepamplitude','desc-cceplatency','desc-pipeline','_ieeglab.mat'}
    tc.verifyTrue(any(contains(names, want{1})), ['Missing export: ' want{1}]);
end
% TSVs round-trip, with BIDS n/a for missing values
f = files{contains(names, 'desc-n1')};
T = readtable(f, 'FileType','text', 'Delimiter','\t', 'TreatAsMissing','n/a');
tc.verifyEqual(height(T), height(tc.TestData.E.ieeglab.n1.table));
fm = files{contains(names, 'desc-ccepresponse')};
txt = fileread(fm);
tc.verifyTrue(startsWith(txt, sprintf('stimulation_site\t')), 'Matrix header must start with stimulation_site.');
tc.verifySubstring(txt, 'n/a', 'Unmeasured cells must be written as n/a.');
% JSON parses and records provenance
J = jsondecode(fileread(files{contains(names, 'desc-pipeline')}));
tc.verifyEqual(J.generated_by, 'iEEGLAB');
tc.verifyEqual(J.reference, 'CARLA');
tc.verifyTrue(isfield(J, 'ccep_matrix'));
% MAT loads
S = load(files{contains(names, '_ieeglab.mat')});
tc.verifyTrue(all(isfield(S, {'channels','n1','crp','ccep_matrix','pipeline'})));
end

function test_export_respects_formats_and_overwrite(tc)
d = fullfile(tc.TestData.tmp, 'export_tsv_only');
[~, files] = evalc('ieeglab_export(tc.TestData.E, d, struct(''formats'',{{''tsv''}}))');
tc.verifyFalse(any(endsWith(files, {'.json','.mat'})));
tc.verifyError(@() ieeglab_export(tc.TestData.E, d, struct('formats',{{'tsv'}},'overwrite',false,'verbose',false)), ...
    'ieeglab_export:exists');
end

function n = local_name(f)
[~, a, b] = fileparts(f); n = [a b];
end

% ======================= electrode values =======================

function test_topoplot_latency_values(tc)
E = tc.TestData.E;
[vals, info] = ieeglab_topoplot(E, 30, 'draw', false);
[~, k] = min(abs(E.times - 30));
tc.verifyEqual(vals, mean(double(E.data(:, k, :)), 3), 'AbsTol', 1e-9);
tc.verifyTrue(all(info.has_xyz));
[vals2] = ieeglab_topoplot(E, [20 60], 'draw', false);
idx = E.times >= 20 & E.times <= 60;
tc.verifyEqual(vals2, mean(mean(double(E.data(:, idx, :)), 3), 2), 'AbsTol', 1e-9);
tc.verifyError(@() ieeglab_topoplot(E, 99999, 'draw', false), 'ieeglab_topoplot:latencyOutside');
tc.verifyEmpty(findall(0,'Type','figure'));
end

function test_topoplot_connectivity_values(tc)
E = tc.TestData.E; M = E.ieeglab.ccep_matrix;
v = ieeglab_topoplot(E, 'in_degree', 'draw', false);
tc.verifyEqual(v, M.in_degree);
s = M.sites{1};
a = ieeglab_topoplot(E, 'n1_amplitude', 'site', s, 'draw', false);
tc.verifyEqual(a(:)', M.amplitude_uv(1,:), 'AbsTol', 1e-12);
tok = ieeglab_site_tokens(s);
rev = strjoin(fliplr(tok), '-');
a2 = ieeglab_topoplot(E, 'n1_amplitude', 'site', rev, 'draw', false);
tc.verifyEqual(a2, a, 'A reversed site name must resolve to the same site.');
tc.verifyError(@() ieeglab_topoplot(E, 'n1_amplitude', 'draw', false), 'ieeglab_topoplot:noSite');
end

function test_is_ieeg(tc)
tc.verifyTrue(ieeglab_is_ieeg(tc.TestData.E));
S = struct('chanlocs', struct('labels', {'Fz','Cz'}, 'type', {'EEG','EEG'}));
tc.verifyFalse(ieeglab_is_ieeg(S));
S.chanlocs(1).type = 'SEEG'; S.chanlocs(2).type = 'SEEG';
tc.verifyTrue(ieeglab_is_ieeg(S));
end

% ======================= orchestration =======================

function test_stats_stage_toggles(tc)
E0 = tc.TestData.E;
E0.ieeglab = rmfield(E0.ieeglab, intersect(fieldnames(E0.ieeglab), {'n1','stats','ccep_matrix'}));
[~, E] = evalc(['ieeglab_stats_subject(E0, struct(''run_n1'',true,''run_crp'',false,' ...
    '''run_matrix'',false,''min_trials'',8,''verbose'',false))']);
tc.verifyTrue(isfield(E.ieeglab, 'n1'));
tc.verifyFalse(isfield(E.ieeglab, 'stats'), 'CRP was switched off.');
tc.verifyFalse(isfield(E.ieeglab, 'ccep_matrix'), 'The matrix was switched off.');
d = fullfile(tc.TestData.tmp, 'orchestrated');
[~, E] = evalc(['ieeglab_stats_subject(E0, struct(''run_n1'',true,''run_crp'',false,' ...
    '''run_matrix'',true,''min_trials'',8,''export_dir'',d,''verbose'',false))']);
tc.verifyTrue(isfield(E.ieeglab, 'ccep_matrix'));
tc.verifyNotEmpty(dir(fullfile(d, '*_desc-ccepresponse_ieeglab.tsv')), 'export_dir must trigger an export.');
tc.verifyEmpty(findall(0,'Type','figure'));
end
