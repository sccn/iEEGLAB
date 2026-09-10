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

tc.TestData.seegDir = fullfile(tc.TestData.root, 'tutorial', 'dataset_seeg');
tc.TestData.ecogDir = fullfile(tc.TestData.root, 'tutorial', 'dataset_ecog');
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
[~, out] = ieeglab_car(E, struct('car_method','carla','car_nboot',20,'verbose',false));

tc.verifyNotEmpty(out, 'No stimulation-site groups were formed.');
for g = 1:numel(out)
    tc.verifyEmpty(intersect(out(g).excluded_channels, out(g).car_channels), ...
        sprintf('Site %s: a stimulated contact is inside its own CAR reference set.', out(g).group));
    tc.verifyNumElements(out(g).excluded_channels, 2, ...
        sprintf('Site %s should exclude exactly the stimulated pair.', out(g).group));
end
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

function test_ecog_non_ccep_pipeline(tc)
% The eCoG dataset has ordinary condition labels, not stimulation pairs.
d = tc.TestData.ecogDir;
EEG = pop_loadset('filename','sub-02_ses-01_task-visual_run-01_ieeg.set','filepath',d);
elecs = readtable(fullfile(d,'sub-02_ses-01_electrodes.tsv'),'FileType','text','Delimiter','\t');
ev    = readtable(fullfile(d,'sub-02_ses-01_task-visual_run-01_events.tsv'),'FileType','text','Delimiter','\t');
[~, EEG] = evalc('get_elec_coor(EEG, elecs)');
EEG.event = [];
for i = 1:height(ev)
    t = ev.trial_type(i);
    if iscell(t), t = t{1}; end
    if ~ischar(t), t = num2str(t); end
    EEG.event(i).type = t;
    EEG.event(i).latency = ev.onset(i)*EEG.srate + 1;
end
EEG = eeg_checkset(EEG,'eventconsistency');
EEG.ieeglab.opt.events = ev;

opt = defaultOpt();
opt.epoch_window = [-200 600];
opt.car_nboot = 15;
[~, E] = evalc('ieeglab_preprocess(EEG, opt)');
tc.verifyGreaterThan(E.trials, 1);
tc.verifyEqual(E.nbchan, 96);
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
