function tests = test_ieeglab_smoke
% test_ieeglab_smoke - Regression tests for iEEGLAB, runnable headless.
%
% Run with:
%   results = runtests('tests/test_ieeglab_smoke.m'); disp(results)
%
% Or from a CI runner:
%   matlab -batch "addpath('<eeglab>'); eeglab nogui; addpath(genpath(pwd)); ...
%                  r = runtests('tests/test_ieeglab_smoke.m'); ...
%                  assertSuccess(r);"
%
% Every test here is headless: none of them opens a dialog or a figure. That is
% only possible because ieeglab_load, ieeglab_preprocess, ieeglab_car and
% ieeglab_stats_subject all take an options struct.
%
% NOTE on CI: MATLAB's -batch mode on a machine with no display can hang on
% graphics calls (imagesc + colorbar was observed to hang indefinitely on
% Windows). Run CI with -softwareopengl or under a virtual display, and keep
% 'plot' false in every option struct.
%
% Cedric Cannard, iEEGLAB, 2026

tests = functiontests(localfunctions);
end

% ======================= fixtures =======================

function setupOnce(tc)
tc.TestData.root = fileparts(fileparts(mfilename('fullpath')));
tc.assumeNotEmpty(which('eeglab'), 'EEGLAB is not on the MATLAB path.');
if isempty(which('pop_loadset'))
    evalc('eeglab nogui');
end
addpath(tc.TestData.root);
addpath(fullfile(tc.TestData.root, 'functions'));

tc.TestData.seegDir = fullfile(tc.TestData.root, 'tutorial', 'seeg');
tc.TestData.ecogDir = fullfile(tc.TestData.root, 'tutorial', 'ecog');
set(0, 'DefaultFigureVisible', 'off');
end

function teardown(~)
close all force
end

function EEG = loadSeeg(tc, doEvents)
% sEEG CCEP dataset with coordinates and, optionally, events attached.
d = tc.TestData.seegDir;
EEG = pop_loadset('filename','sub-02_ses-ieeg01_task-ccep_run-01_ieeg.set','filepath',d);
elecs = readtable(fullfile(d,'sub-02_ses-ieeg01_electrodes.tsv'), ...
    'FileType','text','Delimiter','\t');
[~, EEG] = evalc('get_elec_coor(EEG, elecs)');
EEG.event = [];
if nargin > 1 && doEvents
    ev = readtable(fullfile(d,'sub-02_ses-ieeg01_task-ccep_run-01_events.tsv'), ...
        'FileType','text','Delimiter','\t');
    for i = 1:height(ev)
        EEG.event(i).type    = ev.electrical_stimulation_site{i};
        EEG.event(i).latency = ev.onset(i)*EEG.srate + 1;
    end
    EEG.ieeglab.opt.events = ev;
end
EEG = eeg_checkset(EEG, 'eventconsistency');
EEG = eeg_checkset(EEG);
end

function opt = defaultOpt()
opt = struct('apply_highpass',true, 'highpass',0.5, ...
             'apply_notch',false, 'apply_lowpass',false, 'apply_ds',false, ...
             'remove_no_coords',true, 'remove_rare_cond',false, ...
             'apply_epoch',true, 'epoch_window',[-500 1000], ...
             'apply_car',true, 'car_method','carla', 'car_nboot',20, ...
             'apply_baseline',true, 'baseline_period',[-500 -50], ...
             'plot',false, 'verbose',false);
end

% ======================= tests =======================

function test_no_syntax_errors(tc)
% Every shipped .m file must parse and have no undefined variables.
root = tc.TestData.root;
files = [dir(fullfile(root,'*.m')); dir(fullfile(root,'functions','*.m'))];
bad = {};
for i = 1:numel(files)
    f = fullfile(files(i).folder, files(i).name);
    m = checkcode(f, '-id');
    for j = 1:numel(m)
        if any(strcmp(m(j).id, {'SYNER','UNDEF','NODEF','PFOR'}))
            bad{end+1} = sprintf('%s:%d [%s] %s', files(i).name, m(j).line, m(j).id, m(j).message); %#ok<AGROW>
        end
    end
end
tc.verifyEmpty(bad, sprintf('Code Analyzer errors:\n%s', strjoin(bad, newline)));
end

function test_every_menu_callback_resolves(tc)
% Regression: ieeglab_stats_subject was referenced by the menu but never
% written, so the fourth menu item raised an error for every user.
for f = {'ieeglab_load','ieeglab_vis_elec','ieeglab_preprocess','ieeglab_stats_subject'}
    tc.verifyNotEmpty(which(f{1}), ...
        sprintf('Menu callback %s does not resolve on the path.', f{1}));
end
end

function test_install_check_passes(tc)
[ok, report] = ieeglab_check_install('quiet');
missing = {report(~[report.found] & strcmp({report.kind},'required')).name};
tc.verifyTrue(ok, sprintf('Missing required dependencies: %s', strjoin(missing, ', ')));
end

function test_install_check_needs_mef_mex(tc)
% A help-only read_mef_session_metadata.m (what an Apple Silicon Mac sees: the
% plugin has no MEX file for it) must count as missing, not as installed.
stub = tempname; mkdir(stub);
fid = fopen(fullfile(stub, 'read_mef_session_metadata.m'), 'w');
fprintf(fid, 'function metadata = read_mef_session_metadata(varargin)\n%% help only\n');
fclose(fid);
addpath(stub, '-begin'); c = onCleanup(@() cellfun(@feval, {@() rmpath(stub), @() rmdir(stub, 's')}));
[~, report] = ieeglab_check_install('quiet');
r = report(strcmp({report.name}, 'read_mef_session_metadata'));
tc.verifyFalse(r.found);
tc.verifySubstring(r.hint, 'matmef#building-from-source');
msg = '';
try, ieeglab_load_mefd(tempdir); catch e, msg = e.message; end
tc.verifySubstring(msg, 'no compiled (MEX) reader');
end

function test_no_duplicate_plot_ccep(tc)
% my_quick_ccep_plot.m declared "function plot_ccep", shadowing the real one.
tc.verifyEmpty(which('my_quick_ccep_plot'), ...
    'my_quick_ccep_plot should no longer exist; it duplicated plot_ccep.');
end

function test_plot_ccep_argument_guards(tc)
% The nargin thresholds used to be off by one against the argument positions,
% so plot_ccep(d,t,n,'single') silently drew a heatmap instead of erroring.
%
% This test deliberately exercises only the ARGUMENT VALIDATION, which returns
% before any drawing. Do not add a case here that actually renders: imagesc +
% colorbar hangs indefinitely under matlab -batch without a display, which is
% exactly the hazard noted in this file's header.
D = randn(8,120,10); tv = linspace(-0.5,1,120);
nm = arrayfun(@(k) sprintf('CH%d',k), 1:8, 'UniformOutput', false);
tc.verifyError(@() plot_ccep(D,tv,nm,'single'), 'plot_ccep:noChanIdx', ...
    'A 4-argument ''single'' call must error, not silently draw a heatmap.');

% And the guards themselves must reference the right argument positions.
src = fileread(which('plot_ccep'));
tc.verifySubstring(src, 'nargin < 4 || isempty(view_type)', ...
    'view_type is argument 4 and must be guarded by nargin < 4.');
end

function test_preprocess_runs_headless(tc)
% The keystone regression: this used to hang forever because the GUI guard was
% 'if nargin < 2' inside a one-input function.
EEG = loadSeeg(tc, true);
[E, com] = ieeglab_preprocess(EEG, defaultOpt());
tc.verifyGreaterThan(E.trials, 1, 'Data should be epoched.');
tc.verifyEqual(E.ref, 'CARLA');
tc.verifyNotEmpty(com, 'A history string is required for eeglab_new to store the dataset.');
tc.verifyEmpty(findall(0,'Type','figure'), 'No figures may be opened when plot=false.');
end

function test_stim_contact_never_in_own_reference(tc)
% THE scientific regression. Before the fix, 22 of 30 checked trials had a
% stimulated contact inside the common average, so the stimulation artifact was
% subtracted into every channel.
EEG = loadSeeg(tc, true);
E = ieeglab_preprocess(EEG, defaultOpt());
[~, out] = ieeglab_car(E, struct('car_method','carla','car_nboot',20,'car_neighbors',0,'verbose',false));

tc.verifyNotEmpty(out, 'No stimulation-site groups were formed.');
for g = 1:numel(out)
    tc.verifyEmpty(intersect(out(g).excluded_channels, out(g).car_channels), ...
        sprintf('Site %s: a stimulated contact is inside its own CAR reference set.', out(g).group));
    tc.verifyNumElements(out(g).excluded_channels, 2, ...
        sprintf('Site %s should exclude exactly the stimulated pair with neighbors = 0.', out(g).group));
end
% Default, as the CARLA authors on real data: also 2 contacts on each side, same lead
[~, out] = ieeglab_car(E, struct('car_method','carla','car_nboot',20,'verbose',false));
labels = {E.chanlocs.labels};
for g = 1:numel(out)
    tok = regexp(out(g).group, '^RA(\d+)-RA(\d+)$', 'tokens', 'once');
    n = str2double(tok);
    want = labels(ismember(labels, compose('RA%d', max(1, min(n) - 2):max(n) + 2)));
    tc.verifyEqual(sort(labels(out(g).excluded_channels)), sort(want), ...
        sprintf('Site %s: default CARLA must leave out the pair and 2 contacts on each side.', out(g).group));
    tc.verifyEmpty(intersect(out(g).excluded_channels, out(g).car_channels));
end
end

function test_freesurfer_surface_read(tc)
% A FreeSurfer triangle file (lh.pial) is read like a GIfTI.
d = tempname; mkdir(d); c = onCleanup(@() rmdir(d, 's'));
[V, F] = local_tetra();
local_write_fs(fullfile(d, 'lh.pial'), V, F);
g = ieeglab_read_gifti(fullfile(d, 'lh.pial'));
tc.verifyEqual(g.vertices, V);
tc.verifyEqual(g.faces, F + 1);
end

function test_brain_matches_coordinate_space(tc)
% Which brain is drawn, by the space in coordsystem.json: the participant's own
% pial surfaces, fsaverage for fsaverage / MNI coordinates, none otherwise.
d = tempname; mkdir(d); c = onCleanup(@() rmdir(d, 's'));
[V, F] = local_tetra();
fsdir = fullfile(d, 'fs', 'fsaverage', 'surf'); mkdir(fsdir);
local_write_fs(fullfile(fsdir, 'lh.pial'), V, F); local_write_fs(fullfile(fsdir, 'rh.pial'), V, F);
ieeg = fullfile(d, 'ds', 'sub-01', 'ses-1', 'ieeg'); mkdir(ieeg);
EEG = struct('filepath', ieeg, 'subject', 'sub-01', ...
    'chanlocs', struct('labels', {'A1', 'A2'}, 'X', {2, 1}, 'Y', {2, 3}, 'Z', {2, 2}));
env = {'SUBJECTS_DIR', 'FREESURFER_HOME', 'USERPROFILE', 'HOME'};
old = cellfun(@getenv, env, 'UniformOutput', false);
r = onCleanup(@() cellfun(@setenv, env, old));
setenv('SUBJECTS_DIR', fullfile(d, 'fs')); setenv('FREESURFER_HOME', ''); setenv('USERPROFILE', d); setenv('HOME', d);
fsFiles = fullfile(fsdir, {'lh.pial', 'rh.pial'});

local_coordsystem(ieeg, 'fsaverage');
[f, info] = ieeglab_brain_surfaces(EEG);
tc.verifyEqual(f, fsFiles); tc.verifyEqual(info.kind, 'fsaverage'); tc.verifyTrue(info.fit_ok);

local_coordsystem(ieeg, 'MNI152NLin2009cAsym');
[f, info] = ieeglab_brain_surfaces(EEG);
tc.verifyEqual(f, fsFiles); tc.verifyEqual(info.kind, 'fsaverage_approx'); tc.verifyTrue(info.template_ok);

local_coordsystem(ieeg, 'ACPC');
[f, info] = ieeglab_brain_surfaces(EEG);
tc.verifyEmpty(f, 'Participant coordinates must not be drawn on a template.');
tc.verifyEqual(info.kind, 'none'); tc.verifyFalse(info.template_ok);
tc.verifySubstring(info.message, 'will not line up');

own = fullfile(d, 'ds', 'derivatives', 'freesurfer', 'sub-01', 'surf'); mkdir(own);
local_write_fs(fullfile(own, 'lh.pial'), V, F); local_write_fs(fullfile(own, 'rh.pial'), V, F);
[f, info] = ieeglab_brain_surfaces(EEG);
tc.verifyEqual(f, fullfile(own, {'lh.pial', 'rh.pial'})); tc.verifyEqual(info.kind, 'participant');
rmdir(fullfile(d, 'ds', 'derivatives'), 's');

delete(fullfile(ieeg, '*_coordsystem.json'));             % no space given: fsaverage only if the contacts fit it
[f, info] = ieeglab_brain_surfaces(EEG);
tc.verifyEqual(f, fsFiles); tc.verifyEqual(info.kind, 'fsaverage');
far = EEG; [far.chanlocs.X] = deal(200);
tc.verifyEmpty(ieeglab_brain_surfaces(far));

local_coordsystem(ieeg, 'fsaverage');                      % contacts outside: said so
[~, info] = ieeglab_brain_surfaces(far);
tc.verifyFalse(info.fit_ok); tc.verifySubstring(info.message, 'outside');

setenv('SUBJECTS_DIR', '');                                % fsaverage not installed: how to get it
[f, info] = ieeglab_brain_surfaces(EEG);
tc.verifyEmpty(f); tc.verifyTrue(info.template_ok); tc.verifySubstring(info.message, 'fetch_fsaverage');
end

function test_crp_window_end_tolerance(tc)
% EEGLAB ends a [-1000 1000] ms epoch one sample before 1000 ms (and rounds at
% some rates): the default CRP window [15 1000] must not be reported as cut.
EEG = loadSeeg(tc, true);
opt = defaultOpt(); opt.apply_car = false; opt.epoch_window = [-1000 1000];
opt.apply_ds = true; opt.downsample = 600;
E = ieeglab_preprocess(EEG, opt);
txt = evalc('E2 = ieeglab_stats_subject(E, struct(''run_matrix'', false, ''plot'', false, ''n_perm'', 50, ''verbose'', false));');
tc.verifyFalse(contains(txt, 'ends after the epoch'), 'A window ending at the last sample is not cut.');
tc.verifyEqual(E2.ieeglab.stats.opt.crp_window, [15 1000]);
opt.epoch_window = [-1000 800];                            % a really shorter epoch: said so
E = ieeglab_preprocess(EEG, opt);
txt = evalc('ieeglab_stats_subject(E, struct(''run_matrix'', false, ''plot'', false, ''n_perm'', 50, ''verbose'', false));');
tc.verifyTrue(contains(txt, 'ends after the epoch'));
end

function test_crp_published_test_and_scope(tc)
% crp_test 'ttest' uses run_CRP's t-test at tau_R; 'by' is Benjamini-Yekutieli;
% correct_scope 'subject' corrects across all pairs, 'site' within each site.
EEG = loadSeeg(tc, true);
opt = defaultOpt(); opt.apply_car = false;
E = ieeglab_preprocess(EEG, opt);
base = struct('run_matrix', false, 'plot', false, 'verbose', false, 'n_perm', 200);
Et = ieeglab_stats_subject(E, setfield(setfield(base, 'crp_test', 'ttest'), 'correct', 'by')); %#ok<SFLD>
T = Et.ieeglab.stats.table; S = Et.ieeglab.stats.crp;
for k = 1:numel(S)
    i = find(strcmp(T.site, S(k).site) & strcmp(T.channel, S(k).channel));
    tc.verifyEqual(T.p(i), S(k).projs.p_value_tR, 'AbsTol', 1e-12);
end
n = height(T); [ps, o] = sort(T.p); q = ps .* n ./ (1:n)' * sum(1 ./ (1:n));
q = min(1, flipud(cummin(flipud(q)))); want = zeros(n, 1); want(o) = q;
tc.verifyEqual(T.p_adj, want, 'AbsTol', 1e-12, 'BY across all pairs');
Es = ieeglab_stats_subject(E, setfield(setfield(setfield(base, 'crp_test', 'ttest'), 'correct', 'by'), 'correct_scope', 'site')); %#ok<SFLD>
Ts = Es.ieeglab.stats.table;
for s = unique(Ts.site)'
    m = strcmp(Ts.site, s{1}); p = Ts.p(m); n = numel(p); [ps, o] = sort(p);
    q = ps .* n ./ (1:n)' * sum(1 ./ (1:n)); q = min(1, flipud(cummin(flipud(q)))); want = zeros(n, 1); want(o) = q;
    tc.verifyEqual(Ts.p_adj(m), want, 'AbsTol', 1e-12, 'BY within each site');
end
tc.verifyError(@() ieeglab_stats_subject(E, setfield(base, 'crp_test', 'nope')), 'ieeglab_stats_subject:badTest'); %#ok<SFLD>
end

function test_plot_crp_refuses_stale_results(tc)
% The CRP figure must not mix the current data with results fitted on other data.
EEG = loadSeeg(tc, true);
opt = defaultOpt(); opt.apply_car = false;
E = ieeglab_preprocess(EEG, opt);
E = ieeglab_stats_subject(E, struct('run_matrix', false, 'plot', false, 'verbose', false, 'n_perm', 50));
T = E.ieeglab.stats.table;
[~, info] = ieeglab_plot_crp(E, T.site{1}, lower(T.channel{1}), 'draw', false);   % contact names: any case
tc.verifyEqual(info.channel, T.channel{1});
sites = ieeglab_epoch_sites(E);
first = find(sites == string(T.site{1}), 1);
[~, E2] = evalc('pop_select(E, ''notrial'', first)');
tc.verifyError(@() ieeglab_plot_crp(E2, T.site{1}, T.channel{1}, 'draw', false), 'ieeglab_plot_crp:stale');
end

function test_freesurfer_scanner_shift(tc)
% FreeSurfer vertices are moved by the cras of the file's volume information
d = tempname; mkdir(d); c = onCleanup(@() rmdir(d, 's'));
[V, F] = local_tetra();
f = fullfile(d, 'lh.pial');
local_write_fs(f, V, F);
fid = fopen(f, 'a', 'ieee-be'); fwrite(fid, [0 0 0 2 0 0 0 20], 'uint8');
fprintf(fid, 'valid = 1  # volume info valid\ncras   = 1.5 -2 3\n'); fclose(fid);
g = ieeglab_read_gifti(f);
tc.verifyEqual(g.vertices, V + [1.5 -2 3], 'AbsTol', 1e-6);
tc.verifyEqual(g.c_ras, [1.5 -2 3]);
end

function test_remembered_surfaces_stay_in_their_dataset(tc)
% sub-02 of another dataset is another person: its surfaces are not offered
had = ispref('iEEGLAB', 'surf_files'); oldPref = [];
if had, oldPref = getpref('iEEGLAB', 'surf_files'); end
r = onCleanup(@() local_restore_pref(had, oldPref));
d = tempname; mkdir(d); c = onCleanup(@() rmdir(d, 's'));
fA = fullfile(d, 'A', 'pial.L.surf.gii'); mkdir(fileparts(fA)); fclose(fopen(fA, 'w'));
EA = struct('filepath', fullfile(d, 'A', 'sub-02', 'ieeg'), 'filename', 'sub-02_task-x_ieeg.set', 'subject', '');
EB = struct('filepath', fullfile(d, 'B', 'sub-02', 'ieeg'), 'filename', 'sub-02_task-x_ieeg.set', 'subject', '');
ieeglab_remembered_surfaces(EA, {fA});
tc.verifyEqual(ieeglab_remembered_surfaces(EA), {fA});
tc.verifyEmpty(ieeglab_remembered_surfaces(EB));
end

function test_coordsystem_of_the_loaded_electrodes(tc)
% With several coordinate systems next to the data, the one paired with the
% electrodes file that was loaded (same name, space- entity) decides.
d = tempname; mkdir(d); c = onCleanup(@() rmdir(d, 's'));
ieeg = fullfile(d, 'sub-01', 'ses-1', 'ieeg'); mkdir(ieeg);
local_json(fullfile(ieeg, 'sub-01_ses-1_coordsystem.json'), 'ACPC');
local_json(fullfile(ieeg, 'sub-01_ses-1_space-fsaverage_coordsystem.json'), 'fsaverage');
EEG = struct('filepath', ieeg, 'subject', 'sub-01', 'chanlocs', struct('labels', {'A1'}, 'X', {1}, 'Y', {1}, 'Z', {1}));
EEG.ieeglab.opt.elec_tsv = fullfile(ieeg, 'sub-01_ses-1_space-fsaverage_electrodes.tsv');
[~, info] = ieeglab_brain_surfaces(EEG, 'check_fit', false);
tc.verifyEqual(info.space, 'fsaverage');
EEG.ieeglab.opt.elec_tsv = fullfile(ieeg, 'sub-01_ses-1_electrodes.tsv');
[~, info] = ieeglab_brain_surfaces(EEG, 'check_fit', false);
tc.verifyEqual(info.space, 'ACPC');
end

function local_json(f, space)
fid = fopen(f, 'w'); fprintf(fid, '{"iEEGCoordinateSystem": "%s"}', space); fclose(fid);
end

function local_restore_pref(had, old)
if had, setpref('iEEGLAB', 'surf_files', old);
elseif ispref('iEEGLAB', 'surf_files'), rmpref('iEEGLAB', 'surf_files'); end
end

function test_results_repeat(tc)
% Fixed seeds: the permutation p-values and CARLA's reference are the same on
% every run, and MATLAB's global random state is left alone.
EEG = loadSeeg(tc, true);
opt = defaultOpt(); opt.apply_car = false;
E = ieeglab_preprocess(EEG, opt);
rng(1); g0 = rng;
[~, R1] = evalc('pop_ieeglab_reref(E, ''method'', ''carla'')');
[~, R2] = evalc('pop_ieeglab_reref(E, ''method'', ''carla'')');
tc.verifyEqual(R1.data, R2.data, 'CARLA must give the same reference on every run.');
o = struct('run_matrix', false, 'plot', false, 'verbose', false, 'n_perm', 200);
S1 = ieeglab_stats_subject(R1, o); S2 = ieeglab_stats_subject(R1, o);
tc.verifyEqual(S1.ieeglab.stats.table.p, S2.ieeglab.stats.table.p);
S3 = ieeglab_stats_subject(R1, setfield(o, 'seed', 7)); %#ok<SFLD>
tc.verifyNotEqual(S3.ieeglab.stats.table.p, S1.ieeglab.stats.table.p, 'Another seed, other permutations.');
tc.verifyEqual(rng, g0, 'The global random state must not change.');
end

function test_unflat_baseline_reported(tc)
% RA14 of the tutorial carries a slow deflection in every epoch (average 51 uV
% before the pulse, median contact 8 uV) without being noisy trial by trial.
% It is reported, and kept; a 0.5 Hz high-pass does not remove it (43 uV).
EEG = loadSeeg(tc, true);
opt = defaultOpt(); opt.apply_highpass = false; opt.apply_car = false;
E = ieeglab_preprocess(EEG, opt);
tc.verifyEqual(E.ieeglab.baseline_not_flat.label, {'RA14'});
tc.verifyTrue(any(strcmp({E.chanlocs.labels}, 'RA14')), 'A reported contact must be kept.');
E = ieeglab_preprocess(EEG, setfield(defaultOpt(), 'apply_car', false)); %#ok<SFLD> with the 0.5 Hz high-pass
tc.verifyEqual(E.ieeglab.baseline_not_flat.label, {'RA14'});
end

function test_carla_reference_size_is_adaptive(tc)
% CARLA chooses the set size from the data. A constant size across sites would
% mean we had silently fallen back to the fixed-fraction method.
EEG = loadSeeg(tc, true);
E = ieeglab_preprocess(EEG, defaultOpt());
[~, out] = ieeglab_car(E, struct('car_method','carla','car_nboot',20,'verbose',false));
sizes = [out.n_sel];
tc.verifyGreaterThan(numel(unique(sizes)), 1, ...
    sprintf('CARLA set sizes are all %d; expected them to vary per stimulation site.', sizes(1)));
end

function test_carla_matches_reference_on_synthetic_data(tc)
% Construct a montage where the answer is known: 10 quiet channels plus 4 with
% a large evoked response. CARLA should exclude the responsive ones.
rng(42);
nQuiet = 10; nResp = 4; T = 200; K = 30;
tt = linspace(-0.1, 0.5, T);
V = 5 * randn(nQuiet + nResp, T, K);                 % shared-scale noise
resp = 200 * exp(-((tt-0.08)/0.03).^2);              % a clear evoked bump
for ch = nQuiet+1 : nQuiet+nResp
    V(ch,:,:) = V(ch,:,:) + repmat(resp, 1, 1, K);
end
[~, ~, stats] = ieeglab_carla(tt, V, 1000, struct('nboot',40,'verbose',false,'notchFirst',false));

tc.verifyLessThanOrEqual(max(stats.chsUsed), nQuiet, ...
    sprintf('CARLA put a responsive channel in the reference: chsUsed = %s', mat2str(stats.chsUsed')));
tc.verifyGreaterThanOrEqual(numel(stats.chsUsed), 2, 'CARLA must use at least 2 channels.');
end

function test_carla_excludes_nan_channels(tc)
% NaN across all trials is the contract for "do not use in the reference".
rng(7);
V = randn(12, 150, 20);
V(3,:,:) = NaN;
V(9,:,:) = NaN;
tt = linspace(-0.1, 0.4, 150);
[Vout, ~, stats] = ieeglab_carla(tt, V, 1000, struct('nboot',20,'verbose',false,'notchFirst',false));
tc.verifyEmpty(intersect(stats.chsUsed, [3 9]), 'NaN channels must not enter the reference.');
tc.verifySize(Vout, size(V));
end

function test_continuous_mode_does_not_error(tc)
% A dataset with no events is legitimate (epilepsy / clinical monitoring) and
% used to hit "No events left after event filtering!".
EEG = loadSeeg(tc, false);
opt = defaultOpt();
[~, E] = evalc('ieeglab_preprocess(EEG, opt)');
tc.verifyEqual(E.trials, 1, 'Continuous data must stay continuous.');
tc.verifyNotEqual(E.pnts, 0);
end

function EEG = loadEcogAsErp(tc)
% ECoG dataset with its stimulation events relabelled as two ordinary
% conditions, which gives an event-related (non-CCEP) design.
d = tc.TestData.ecogDir;
base = 'sub-ccepAgeUMCU02_ses-1';
EEG = pop_loadset('filename',[base '_task-SPESclin_run-041456_ieeg.set'],'filepath',d);
elecs = readtable(fullfile(d,[base '_electrodes.tsv']),'FileType','text','Delimiter','\t');
ev    = readtable(fullfile(d,[base '_task-SPESclin_run-041456_events.tsv']),'FileType','text','Delimiter','\t');
ev    = ev(strcmp(ev.trial_type, 'electrical_stimulation'), :);
[~, EEG] = evalc('get_elec_coor(EEG, elecs)');
EEG.event = [];
cond = {'face', 'house'};
for i = 1:height(ev)
    EEG.event(i).type = cond{mod(i, 2) + 1};
    EEG.event(i).latency = ev.onset(i)*EEG.srate + 1;
end
EEG = eeg_checkset(EEG,'eventconsistency');
end

function test_ecog_non_ccep_pipeline(tc)
% Event types that are ordinary condition labels, not stimulation pairs.
EEG = loadEcogAsErp(tc);
opt = defaultOpt();
opt.epoch_window = [-200 600];
opt.car_nboot = 15;
[~, E] = evalc('ieeglab_preprocess(EEG, opt)');
tc.verifyGreaterThan(E.trials, 1);
tc.verifyEqual(E.nbchan, 20);
end

function test_crp_statistics(tc)
EEG = loadSeeg(tc, true);
E = ieeglab_preprocess(EEG, defaultOpt());
E = ieeglab_stats_subject(E, struct('crp_window',[15 400],'min_trials',8, ...
    'plot',false,'verbose',false));

tc.verifyTrue(isfield(E.ieeglab,'stats'));
T = E.ieeglab.stats.table;
tc.verifyNotEmpty(T, 'CRP produced no results.');
for c = {'site','channel','tR_ms','p','p_adj','significant'}
    tc.verifyTrue(ismember(c{1}, T.Properties.VariableNames), ...
        sprintf('Results table is missing the %s column.', c{1}));
end
% Adjusted p-values must never be smaller than raw ones
tc.verifyGreaterThanOrEqual(T.p_adj, T.p - 1e-12);
% Response durations must sit inside the requested window
tc.verifyGreaterThan(min(T.tR_ms), 0);
tc.verifyLessThanOrEqual(max(T.tR_ms), 400);
% A stimulated contact must not be reported as responding to its own site
for i = 1:height(T)
    tc.verifyFalse(any(strcmpi(T.channel{i}, split(T.site{i},'-'))), ...
        sprintf('%s reported as a response channel for its own stimulation site.', T.channel{i}));
end
% The per-pair CRP figure shows the values of the results table
[~, info] = ieeglab_plot_crp(E, T.site{1}, T.channel{1}, 'draw', false);
tc.verifyEqual(info.tR_ms, T.tR_ms(1), 'AbsTol', 1e-9);
tc.verifyEqual(mean(info.explained_var), T.explained_var(1), 'AbsTol', 1e-9);
tc.verifyEqual(numel(info.shape), numel(info.shape_ms));
tc.verifyEqual(size(info.trials, 2), T.n_trials(1));
[~, iMax] = max(info.profile);
tc.verifyEqual(info.profile_ms(iMax), info.tR_ms, 'AbsTol', 1e-9, 'tau_R must be where the projection profile peaks.');
tc.verifyError(@() ieeglab_plot_crp(E, T.site{1}, 'NOPE', 'draw', false), 'ieeglab_plot_crp:noPair');
tc.verifyEmpty(findall(0,'Type','figure'));
end

function test_car_methods_differ(tc)
% CARLA and the legacy fixed-fraction method must not silently be the same code
% path, and both must leave the data changed.
EEG = loadSeeg(tc, true);
E = ieeglab_preprocess(EEG, setfield(defaultOpt(), 'apply_car', false)); %#ok<SFLD>

Ea = ieeglab_car(E, struct('car_method','carla','car_nboot',20,'verbose',false));
Eb = ieeglab_car(E, struct('car_method','varsubset','verbose',false));

tc.verifyNotEqual(Ea.data, E.data, 'CARLA did not change the data.');
tc.verifyNotEqual(Eb.data, E.data, 'varsubset did not change the data.');
tc.verifyNotEqual(Ea.data, Eb.data, 'The two methods produced identical output.');
tc.verifyEqual(Ea.ref, 'CARLA');
end

function test_car_none_is_a_noop(tc)
EEG = loadSeeg(tc, true);
E = ieeglab_preprocess(EEG, setfield(defaultOpt(), 'apply_car', false)); %#ok<SFLD>
[E2, out] = ieeglab_car(E, struct('car_method','none','verbose',false));
tc.verifyEqual(E2.data, E.data);
tc.verifyEmpty(out);
end

function test_baseline_window_validation(tc)
% The window used to be snapped to the nearest samples with no validation, so a
% window outside the epoch collapsed onto an edge sample and a post-stimulus
% window silently baselined the response against itself.
EEG = loadSeeg(tc, true);
E = ieeglab_preprocess(EEG, setfield(defaultOpt(), 'apply_baseline', false)); %#ok<SFLD>

mk = @(win) setfield(E, 'ieeglab', setfield(E.ieeglab, 'opt', ...
    struct('baseline_period', win, 'baseline_mode', 'subtract', 'baseline_method', 'median'))); %#ok<SFLD>

tc.verifyError(@() ieeglab_rm_baseline(mk([-5000 -4000])), ...
    'ieeglab_rm_baseline:windowOutsideEpoch', ...
    'A window entirely outside the epoch must error, not snap to an edge sample.');

tc.verifyWarning(@() ieeglab_rm_baseline(mk([-200 300])), ...
    'ieeglab_rm_baseline:windowIncludesPostStim', ...
    'A window crossing t=0 must warn: the response would be subtracted from itself.');

tc.verifyWarningFree(@() ieeglab_rm_baseline(mk([-450 -50])));
end

function test_baseline_subtract_is_correct_and_per_trial(tc)
% Subtraction must zero the baseline, preserve the evoked response, and be
% computed per trial rather than pooled across trials.
srate = 1000; t = -500:1000/srate:1000; T = numel(t);
C = 6; N = 20;
rng(3);
X = randn(C,T,N)*20;
bump = 100*exp(-((t-80)/25).^2); bump(t<0) = 0;
X = X + repmat(reshape(bump,1,T,1), C, 1, N);
for k = 1:N, X(:,:,k) = X(:,:,k) + 200*k/N; end   % trial-varying DC drift

E = struct('data',X, 'times',t, 'srate',srate, 'nbchan',C, 'pnts',T, 'trials',N, ...
           'xmin',t(1)/1000, 'xmax',t(end)/1000);
E.chanlocs = struct('labels', arrayfun(@(k) sprintf('C%d',k), 1:C, 'UniformOutput', false));
E.ieeglab.opt = struct('baseline_period',[-500 -50], 'baseline_mode','subtract', 'baseline_method','median');

Eo = ieeglab_rm_baseline(E);
bi = t >= -500 & t <= -50;

tc.verifyLessThan(abs(mean(Eo.data(:,bi,:), 'all')), 1, 'Baseline was not removed.');
tc.verifyGreaterThan(max(mean(Eo.data,3), [], 'all'), 80, 'The evoked response was destroyed.');

% Per-trial: each trial's own baseline must be near zero despite the drift
perTrial = squeeze(mean(Eo.data(:,bi,:), [1 2]));
tc.verifyLessThan(max(abs(perTrial)), 5, ...
    'Baseline looks pooled across trials: per-trial drift survived correction.');
end

function test_divisive_baseline_not_offered(tc)
% Divisive baselining assumes a multiplicative model (Gyurkovics et al. 2021).
% On high-passed time-domain data the baseline mean is near zero, so the ratio
% explodes and its sign is arbitrary. Only subtraction is supported.
srate = 1000; t = -500:1000/srate:1000; T = numel(t);
rng(4);
X = randn(6,T,15)*20;                       % zero-mean, like high-passed data
E = struct('data',X, 'times',t, 'srate',srate, 'nbchan',6, 'pnts',T, 'trials',15, ...
           'xmin',t(1)/1000, 'xmax',t(end)/1000);
E.chanlocs = struct('labels', arrayfun(@(k) sprintf('C%d',k), 1:6, 'UniformOutput', false));
E.ieeglab.opt = struct('baseline_period',[-500 -50], 'baseline_mode','divide', 'baseline_method','median');

tc.verifyError(@() ieeglab_rm_baseline(E), 'ieeglab_rm_baseline:onlySubtract');
end

function test_aperiodic_baseline_is_gone(tc)
EEG = loadSeeg(tc, true);
E = ieeglab_preprocess(EEG, setfield(defaultOpt(), 'apply_baseline', false)); %#ok<SFLD>
E.ieeglab.opt.baseline_method = '1/f';
tc.verifyError(@() ieeglab_rm_baseline(E), 'ieeglab_rm_baseline:removedMethod');
end

function test_carla_refuses_non_ccep(tc)
% CARLA is defined for CCEP only; on other data it must fall back with a warning
% rather than return an uninterpretable reference.
EEG = loadEcogAsErp(tc);
tc.verifyEqual(ieeglab_detect_mode(EEG), 'erp', ...
    'A dataset with condition labels must not be classified as CCEP.');

opt = defaultOpt(); opt.epoch_window = [-200 600]; opt.car_nboot = 10;
[~, E] = evalc('ieeglab_preprocess(EEG, opt)');
tc.verifyEqual(E.ref, 'CAR', ...
    'CARLA must fall back to a plain CAR on non-CCEP data.');
end

function test_mode_detection(tc)
EEG = loadSeeg(tc, true);
tc.verifyEqual(ieeglab_detect_mode(EEG), 'ccep');
EEG.event = [];
tc.verifyEqual(ieeglab_detect_mode(EEG), 'continuous');
end

function [V, F] = local_tetra()
% a small closed surface (mm) and its 0-based faces, as FreeSurfer stores them
V = [0 0 0; 10 0 0; 0 10 0; 0 0 10];
F = [0 1 2; 0 1 3; 0 2 3; 1 2 3];
end

function local_write_fs(f, V, F)
fid = fopen(f, 'w', 'ieee-be');
fwrite(fid, [255 255 254], 'uint8'); fprintf(fid, 'created by test\n\n');
fwrite(fid, [size(V, 1) size(F, 1)], 'int32'); fwrite(fid, V', 'float32'); fwrite(fid, F', 'int32');
fclose(fid);
end

function local_coordsystem(folder, space)
delete(fullfile(folder, '*_coordsystem.json'));
fid = fopen(fullfile(folder, 'sub-01_ses-1_coordsystem.json'), 'w');
fprintf(fid, '{"iEEGCoordinateSystem": "%s"}', space);
fclose(fid);
end
