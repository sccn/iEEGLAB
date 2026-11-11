function tests = test_ieeglab_import_fieldtrip
% FieldTrip import (ieeglab_import_fieldtrip): units of positions and values,
% electrode structures, events, history, file information, and the FieldTrip
% coordinate system names read by ieeglab_brain_surfaces. Headless.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(mfilename('fullpath')));
tc.assumeNotEmpty(which('eeglab'), 'EEGLAB is not on the MATLAB path.');
if isempty(which('pop_loadset')), evalc('eeglab nogui'); end
addpath(root); addpath(fullfile(root, 'functions'));
end

function test_position_units(tc)
% elec.unit m, dm, cm, mm, um and µm all give mm; an unknown unit is an error;
% no unit: estimated from the size of the layout, with a warning naming it.
d = local_ft((-50:100) / 100, 3, 50);
pos = [10 20 30; 40 50 60; 70 80 90];                          % mm
units = {'m', 1e-3; 'dm', 1e-2; 'cm', 0.1; 'mm', 1; 'um', 1e3; [char(181) 'm'], 1e3};
for k = 1:size(units, 1)
    el = local_elec(pos * units{k, 2}, units{k, 1});
    EEG = tc.verifyWarningFree(@() ieeglab_import_fieldtrip(d, 'elec', el, 'verbose', false));
    tc.verifyEqual(local_xyz(EEG), pos, 'AbsTol', 1e-9, units{k, 1});
end
tc.verifyError(@() ieeglab_import_fieldtrip(d, 'elec', local_elec(pos, 'inch'), 'verbose', false), ...
    'ieeglab_import_fieldtrip:unit');
for k = [1 3 4]                                                % no unit, positions in m, cm, mm
    el = rmfield(local_elec(pos * units{k, 2}, units{k, 1}), 'unit');
    EEG = tc.verifyWarning(@() ieeglab_import_fieldtrip(d, 'elec', el, 'verbose', false), ...
        'ieeglab_import_fieldtrip:unitAssumed');
    tc.verifyEqual(local_xyz(EEG), pos, 'AbsTol', 1e-9, units{k, 1});
    [~, txt] = local_import(d, 'elec', el, 'verbose', false);
    tc.verifyTrue(contains(txt, ['assumed ' units{k, 1} ',']), txt);
end
end

function test_value_units(tc)
% chanunit V or mV on values of that size: converted to uV, and said so. A label
% that disagrees with the size (FieldTrip's 'V' on values in uV, as SubjectUCI29):
% values kept as they are, and a warning.
d = local_ft((-50:100) / 100, 3, 50);                          % uV: median absolute deviation ~34
ref = cat(3, d.trial{:});
dv = d; dv.trial = cellfun(@(x) x * 1e-6, d.trial, 'UniformOutput', false); dv.chanunit = repmat({'V'}, 3, 1);
[EEG, txt] = local_import(dv, 'verbose', true);
tc.verifyEqual(double(EEG.data), ref, 'RelTol', 1e-6);
tc.verifyTrue(contains(txt, 'data.chanunit is V (3 channels): values converted to uV'), txt);
tc.verifyFalse(contains(txt, 'Warning'), txt);
dm = d; dm.trial = cellfun(@(x) x * 1e-3, d.trial, 'UniformOutput', false);
dm.elec = local_elec([10 20 30; 40 50 60; 70 80 90], 'mm'); dm.elec.chanunit = repmat({'mV'}, 3, 1);
EEG = tc.verifyWarningFree(@() ieeglab_import_fieldtrip(dm, 'verbose', false));
tc.verifyEqual(double(EEG.data), ref, 'RelTol', 1e-6);
dw = d; dw.elec = dm.elec; dw.elec.chanunit = repmat({'V'}, 3, 1);   % label wrong
EEG = tc.verifyWarning(@() ieeglab_import_fieldtrip(dw, 'verbose', false), 'ieeglab_import_fieldtrip:chanunit');
tc.verifyEqual(EEG.data, single(ref));
[~, txt] = local_import(dw, 'verbose', false);
tc.verifyTrue(contains(txt, 'data.elec.chanunit is V (3 channels), but the values have the size of uV') && ...
    contains(txt, 'the label looks wrong'), txt);
du = d; du.chanunit = 'uV';
EEG = tc.verifyWarningFree(@() ieeglab_import_fieldtrip(du, 'verbose', false));
tc.verifyEqual(EEG.data, single(ref));
end

function test_elec_from_structure_without_positions(tc)
% An elec file or structure with label but no positions: its field elec is used
% (FieldTrip data); without one, a clear error. chanpos before elecpos.
folder = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
pos = [10 20 30; 40 50 60; 70 80 90];
d = local_ft((-50:100) / 100, 3, 50);
data = d; data.elec = local_elec(pos, 'mm');
save(fullfile(folder, 'ft_data.mat'), 'data');
EEG = ieeglab_import_fieldtrip(d, 'elec', fullfile(folder, 'ft_data.mat'), 'verbose', false);
tc.verifyEqual(local_xyz(EEG), pos, 'AbsTol', 1e-9);
EEG = ieeglab_import_fieldtrip(d, 'elec', data, 'verbose', false);
tc.verifyEqual(local_xyz(EEG), pos, 'AbsTol', 1e-9);
freq = struct('label', {d.label}, 'dimord', 'chan_freq', 'freq', 1:4, 'powspctrm', ones(3, 4));
save(fullfile(folder, 'ft_freq.mat'), 'freq');
tc.verifyError(@() ieeglab_import_fieldtrip(d, 'elec', fullfile(folder, 'ft_freq.mat'), 'verbose', false), ...
    'ieeglab_import_fieldtrip:noPositions');
tc.verifyError(@() ieeglab_import_fieldtrip(d, 'elec', freq, 'verbose', false), 'ieeglab_import_fieldtrip:noPositions');
el = local_elec(pos, 'mm'); el.elecpos = pos + 1;
EEG = ieeglab_import_fieldtrip(d, 'elec', el, 'verbose', false);
tc.verifyEqual(local_xyz(EEG), pos, 'AbsTol', 1e-9);
end

function test_epochs_must_contain_time_zero(tc)
% The event of each epoch goes at time 0: epochs that start after it or end
% before it are an error; time 0 on the first sample is fine.
tc.verifyError(@() ieeglab_import_fieldtrip(local_ft((10:100) / 100, 3, 50), 'verbose', false), ...
    'ieeglab_import_fieldtrip:noTimeZero');
tc.verifyError(@() ieeglab_import_fieldtrip(local_ft((-100:-10) / 100, 3, 50), 'verbose', false), ...
    'ieeglab_import_fieldtrip:noTimeZero');
EEG = ieeglab_import_fieldtrip(local_ft((0:100) / 100, 3, 50), 'verbose', false);
tc.verifyEqual([EEG.event.latency], (0:2) * 101 + 1);
tc.verifyEqual(local_m([EEG.epoch.eventlatency]), [0 0 0], 'AbsTol', 1e-9);
end

function test_time_zero_between_samples(tc)
% Time 0 between two samples: fractional latency, at 0 ms in every epoch, kept
% by eeg_checkset. On a sample: an integer latency.
d = local_ft((-50:100) / 100 + 0.004, 3, 50);                  % time 0 at sample 50.6
EEG = ieeglab_import_fieldtrip(d, 'verbose', false);
tc.verifyEqual([EEG.event.latency], (0:2) * 151 + 50.6, 'AbsTol', 1e-9);
tc.verifyEqual(local_m([EEG.epoch.eventlatency]), [0 0 0], 'AbsTol', 1e-6);
E2 = eeg_checkset(EEG, 'eventconsistency');
tc.verifyEqual([E2.event.latency], (0:2) * 151 + 50.6, 'AbsTol', 1e-9);
EEG = ieeglab_import_fieldtrip(local_ft((-50:100) / 100, 3, 50), 'verbose', false);
tc.verifyEqual([EEG.event.latency], (0:2) * 151 + 51);
end

function test_one_trial_is_continuous(tc)
% One trial is continuous data: no event, time from 0, data unchanged.
d = local_ft((0:999) / 100 + 5, 1, 50); d.trialinfo = 7;
[EEG, txt] = local_import(d, 'verbose', false);
tc.verifyEmpty(EEG.event);
tc.verifyEqual([EEG.trials EEG.pnts EEG.xmin], [1 1000 0]);
tc.verifyEqual(EEG.times(end), 9990, 'AbsTol', 1e-9);
tc.verifyEqual(EEG.data, single(d.trial{1}));
tc.verifyFalse(contains(txt, 'Warning'), txt);
end

function test_history_reruns(tc)
% Files by name, structures by the caller's variable name, numbers as numbers;
% re-running the history gives the same dataset. A structure that is not a
% variable cannot be written: no history.
folder = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
d = local_ft((-50:100) / 100, 3, 50);
data = d; file = fullfile(folder, 'ft_data.mat'); save(file, 'data');
el = local_elec([1 2 3; 4 5 6; 7 8 9], 'cm');
[EEG, com] = ieeglab_import_fieldtrip(file, 'elec', el, 'event_type', 4, 'verbose', false);
tc.verifyEqual(com, ['EEG = ieeglab_import_fieldtrip(''' file ''', ''elec'', el, ''event_type'', 4);']);
R = local_rerun(com, el);
tc.verifyEqual(local_xyz(R), local_xyz(EEG));
tc.verifyEqual(local_xyz(R), [10 20 30; 40 50 60; 70 80 90], 'AbsTol', 1e-9);
tc.verifyEqual(R.event(1).type, 4);
[~, com] = ieeglab_import_fieldtrip(d, 'types', {'CH*', 'ECOG'}, 'verbose', false);
tc.verifyEqual(com, 'EEG = ieeglab_import_fieldtrip(d, ''types'', {''CH*'', ''ECOG''});');
s.el = el;
[~, com] = ieeglab_import_fieldtrip(file, 'elec', s.el, 'verbose', false);
tc.verifyEmpty(com);
end

function test_type_patterns_ignore_case(tc)
% Type patterns match labels as positions do, case ignored; first match wins.
d = local_ft((-50:100) / 100, 3, 50);
EEG = ieeglab_import_fieldtrip(d, 'types', {'ch1', 'ecog'; 'c?2', 'seeg'; 'CH*', 'EEG'}, 'verbose', false);
tc.verifyEqual({EEG.chanlocs.type}, {'ECOG', 'SEEG', 'EEG'});
end

function test_missing_positions_listed(tc)
% Channels without a position (NaN in elec, or not in elec) are named.
d = local_ft((-50:100) / 100, 3, 50);
el = local_elec([10 20 30; nan(1, 3); 70 80 90], 'mm'); el.label{3} = 'OTHER';
tc.verifyWarning(@() ieeglab_import_fieldtrip(d, 'elec', el, 'verbose', false), 'ieeglab_import_fieldtrip:noPosition');
[EEG, txt] = local_import(d, 'elec', el, 'verbose', false);
tc.verifyTrue(contains(txt, '2 channel(s) have no position in elec: CH2 CH3'), txt);
tc.verifyEqual(local_xyz(EEG), [10 20 30; nan(2, 3)], 'AbsTol', 1e-9);
end

function test_filepath_from_data_file(tc)
% Data read from a file keep its folder and name (also for a relative name);
% no subject is made up. A structure leaves them empty.
folder = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
data = local_ft((-50:100) / 100, 3, 50);
save(fullfile(folder, 'ft_data.mat'), 'data');
EEG = ieeglab_import_fieldtrip(fullfile(folder, 'ft_data.mat'), 'verbose', false);
tc.verifyEqual(lower(fullfile(EEG.filepath, EEG.filename)), lower(fullfile(folder, 'ft_data.mat')));
tc.verifyEqual(EEG.filename, 'ft_data.mat');
tc.verifyEmpty(EEG.subject);
tc.applyFixture(matlab.unittest.fixtures.CurrentFolderFixture(folder));
EEG = ieeglab_import_fieldtrip('ft_data.mat', 'verbose', false);
tc.verifyEqual(lower(EEG.filepath), lower(folder));
EEG = ieeglab_import_fieldtrip(local_ft((-50:100) / 100, 3, 50), 'verbose', false);
tc.verifyEmpty(EEG.filepath); tc.verifyEmpty(EEG.filename);
end

function test_fieldtrip_coordsys_names(tc)
% FieldTrip's coordsys names as ieeglab_brain_surfaces reads them: acpc is the
% participant's space; mni, spm, mni152 and tal are templates; unknown is not
% given, so a coordsystem.json next to the data is read.
d = local_ft((-50:100) / 100, 3, 50);
EEG = ieeglab_import_fieldtrip(d, 'elec', local_elec([10 20 30; 40 50 60; 70 80 90], 'mm'), 'verbose', false);
cases = {'acpc', 'acpc', false; 'mni', 'MNI', true; 'spm', 'MNI', true; 'mni152', 'MNI152', true; ...
         'tal', 'Talairach', true; 'unknown', '', false};
for k = 1:size(cases, 1)
    EEG.ieeglab.coordsys = cases{k, 1};
    [~, info] = ieeglab_brain_surfaces(EEG, 'check_fit', false);
    tc.verifyEqual(info.space, cases{k, 2}, cases{k, 1});
    tc.verifyEqual(info.template_ok, cases{k, 3}, cases{k, 1});
end
EEG.ieeglab.coordsys = 'spm';
[~, info] = ieeglab_brain_surfaces(EEG, 'check_fit', false);
if ~isempty(ieeglab_fsaverage_surfaces()), tc.verifyEqual(info.kind, 'fsaverage_approx'); end
folder = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
fid = fopen(fullfile(folder, 'sub-01_coordsystem.json'), 'w');
fprintf(fid, '{"iEEGCoordinateSystem": "ACPC", "iEEGCoordinateUnits": "mm"}');
fclose(fid);
EEG.filepath = folder; EEG.ieeglab.coordsys = 'unknown';
[~, info] = ieeglab_brain_surfaces(EEG, 'check_fit', false);
tc.verifyEqual(info.space, 'ACPC');
end

% ===================== helpers =====================

function d = local_ft(time, nTr, amp)
% FieldTrip raw data: channels CH1-CH3 at 100 Hz, nTr trials on one time axis
rng(11);
d = struct('label', {{'CH1'; 'CH2'; 'CH3'}}, 'fsample', 100, ...
    'trial', {arrayfun(@(k) amp * randn(3, numel(time)), 1:nTr, 'UniformOutput', false)}, ...
    'time', {repmat({time}, 1, nTr)});
end

function el = local_elec(pos, unit)
el = struct('label', {{'CH1'; 'CH2'; 'CH3'}}, 'chanpos', pos, 'unit', unit);
end

function xyz = local_xyz(EEG)
xyz = nan(numel(EEG.chanlocs), 3);
for c = 1:numel(EEG.chanlocs)
    v = [EEG.chanlocs(c).X EEG.chanlocs(c).Y EEG.chanlocs(c).Z];
    if numel(v) == 3, xyz(c, :) = v; end
end
end

function [EEG, txt] = local_import(varargin)
% The import, with what it prints (warnings included)
EEG = [];
txt = evalc('EEG = ieeglab_import_fieldtrip(varargin{:});');
% one line, single spaces: MATLAB on Linux wraps long warnings at 80 columns and
% brackets them with backspaces
txt = regexprep(strrep(txt, char(8), ''), '\s+', ' ');
end

function EEG = local_rerun(com, el) %#ok<INUSD>
% The history line, run where its variable el exists
EEG = [];
eval(com);
end

function y = local_m(x)
if iscell(x), y = cell2mat(x); else, y = x; end
end
