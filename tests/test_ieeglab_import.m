function tests = test_ieeglab_import
% ieeglab_import: reading part of a raw recording (channel types, stimulation
% sites, time range) and loading its sidecars afterwards. Headless.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(mfilename('fullpath')));
tc.assumeNotEmpty(which('eeglab'), 'EEGLAB is not on the MATLAB path.');
if isempty(which('pop_loadset')), evalc('eeglab nogui'); end
addpath(root); addpath(fullfile(root, 'functions'));
tc.TestData.dir = fullfile(root, 'tutorial', 'seeg');
tc.TestData.set = fullfile(tc.TestData.dir, 'sub-02_ses-ieeg01_task-ccep_run-01_ieeg.set');
tc.TestData.O = pop_loadset('filename', 'sub-02_ses-ieeg01_task-ccep_run-01_ieeg.set', ...
    'filepath', tc.TestData.dir, 'verbose', 'off');
tc.TestData.tmp = tempname; mkdir(tc.TestData.tmp);
end

function teardownOnce(tc)
if isfolder(tc.TestData.tmp), rmdir(tc.TestData.tmp, 's'); end
end

function test_info_reads_no_data(tc)
info = ieeglab_import(tc.TestData.set, 'info', true);
tc.verifyEqual(info.labels, {tc.TestData.O.chanlocs.labels});
tc.verifyEqual(info.srate, tc.TestData.O.srate);
tc.verifyEqual(info.pnts, tc.TestData.O.pnts);
tc.verifyEqual(info.sites, {'RA2-RA3','RA3-RA4','RA4-RA5'});
tc.verifyEqual(sum(info.site_counts), 34);
tc.verifyTrue(all(strcmp(info.types, 'SEEG')));
end

function test_whole_recording_is_single_and_identical(tc)
O = tc.TestData.O;
[~, E, com] = evalc('ieeglab_import(tc.TestData.set)');
tc.verifyClass(E.data, 'single');
tc.verifyEqual(E.data, single(O.data));
tc.verifyEqual({E.chanlocs.labels}, {O.chanlocs.labels});
tc.verifyEqual(E.subject, 'sub-02');
tc.verifySubstring(com, 'ieeglab_import(');
EEG = [];
evalc(com);                                          % the history line runs
tc.verifyEqual(EEG.data, E.data);
end

function test_sites_read_only_their_stretch(tc)
% Only the stretch around the pulses of one site, and after loading the
% sidecars only that site's events, at the right samples.
O = tc.TestData.O;
[~, E] = evalc('ieeglab_import(tc.TestData.set, ''sites'', {''RA3-RA4''})');
s0 = E.etc.ieeglab_import.first_sample;
tc.verifyLessThan(E.pnts, O.pnts / 2);
tc.verifyEqual(E.data, single(O.data(:, s0 + (1:E.pnts))));
[~, L] = evalc('ieeglab_load(E, struct(''event_field'', ''electrical_stimulation_site''))');
ty = cellfun(@(x) char(string(x)), {L.event.type}, 'UniformOutput', false);
tc.verifyEqual(unique(ty), {'RA3-RA4'});
ev = readtable(fullfile(tc.TestData.dir, 'sub-02_ses-ieeg01_task-ccep_run-01_events.tsv'), 'FileType','text', 'Delimiter','\t');
on = ev.onset(strcmp(ev.electrical_stimulation_site, 'RA3-RA4') & ~strcmpi(ev.status, 'bad'));
tc.verifyEqual(numel(L.event), numel(on));
tc.verifyEqual(([L.event.latency] - 1) / L.srate + s0 / L.srate, on', 'AbsTol', 0.5 / L.srate);
end

function test_time_range_shifts_events(tc)
[~, E] = evalc('ieeglab_import(tc.TestData.set, ''time_range'', [30 90])');
tc.verifyEqual(E.pnts, 60 * E.srate);
tc.verifyEqual(E.etc.ieeglab_import.onset_shift_s, 30);
[~, L] = evalc('ieeglab_load(E, struct(''event_field'', ''electrical_stimulation_site''))');
ev = readtable(fullfile(tc.TestData.dir, 'sub-02_ses-ieeg01_task-ccep_run-01_events.tsv'), 'FileType','text', 'Delimiter','\t');
in = ev.onset >= 30 & ev.onset <= 90 - 1 / E.srate & ~strcmpi(ev.status, 'bad') & ~strcmpi(ev.electrical_stimulation_site, 'n/a');
tc.verifyEqual(numel(L.event), nnz(in));
tc.verifyEqual(([L.event.latency] - 1) / L.srate, ev.onset(in)' - 30, 'AbsTol', 0.5 / L.srate);
end

function test_channel_types_and_labels(tc)
% Brain channels only by default; 'all' and explicit labels otherwise.
d = fullfile(tc.TestData.tmp, 'types'); copyfile(tc.TestData.dir, d);
f = fullfile(d, 'sub-02_ses-ieeg01_task-ccep_run-01_channels.tsv');
ct = readtable(f, 'FileType','text', 'Delimiter','\t', 'TextType','char');
ct.type(2) = {'EEG'}; ct.type(5) = {'ECG'};
writetable(ct, f, 'FileType','text', 'Delimiter','\t');
setf = fullfile(d, 'sub-02_ses-ieeg01_task-ccep_run-01_ieeg.set');
[~, E] = evalc('ieeglab_import(setf, ''time_range'', [0 2])');
tc.verifyEqual({E.chanlocs.labels}, ct.name(setdiff(1:height(ct), [2 5]))');
[~, E] = evalc('ieeglab_import(setf, ''time_range'', [0 2], ''types'', ''all'')');
tc.verifyEqual(E.nbchan, height(ct));
[~, E] = evalc('ieeglab_import(setf, ''time_range'', [0 2], ''types'', {''EEG'',''ECG''})');
tc.verifyEqual({E.chanlocs.labels}, ct.name([2 5])');
[~, E] = evalc('ieeglab_import(setf, ''time_range'', [0 2], ''channels'', {''RA3'',''RA1''})');
tc.verifyEqual({E.chanlocs.labels}, {'RA3','RA1'});
tc.verifyEqual(E.data, single(tc.TestData.O.data([3 1], 1:2 * E.srate)));
tc.verifyError(@() ieeglab_import(setf, 'types', {'DBS'}, 'verbose', false), 'ieeglab_import:noType');
tc.verifyError(@() ieeglab_import(setf, 'sites', {'ZZ1-ZZ2'}, 'verbose', false), 'ieeglab_import:unknownSite');
end

function test_mef3_round_trip(tc)
% A small MEF3 session written here, read back by type and by site: values,
% channel order (channels.tsv, not the file's own) and window.
tc.assumeTrue(ieeglab_is_mex('read_mef_ts_data'), 'MEF3 reader not compiled for this platform.');
tc.assumeTrue(ieeglab_is_mex('write_mef_ts_segment_data'), 'MEF3 writer not compiled for this platform.');
tc.assumeNotEmpty(which('writeMef3'), 'MEF3 plugin not installed.');
d = fullfile(tc.TestData.tmp, 'mef'); mkdir(d);
fs = 1000; n = 30 * fs;
names = {'LB1','LA2','EKG','LA1'};                   % file order differs from channels.tsv below
x = round(1000 * sin((1:n) / 37));
data = [x; -x; 2 * x; x + 5];
sess = fullfile(d, 'sub-01_task-ccep_run-01_ieeg.mefd');
evalc('writeMef3(sess, data, fs, names)');
local_write(fullfile(d, 'sub-01_task-ccep_run-01_channels.tsv'), ...
    sprintf('name\ttype\tunits\nLA1\tSEEG\tuV\nLA2\tSEEG\tuV\nLB1\tSEEG\tuV\nEKG\tECG\tuV\n'));
local_write(fullfile(d, 'sub-01_task-ccep_run-01_events.tsv'), ...
    sprintf('onset\tduration\telectrical_stimulation_site\n5\t0\tLA1-LA2\n10\t0\tLA1-LA2\n20\t0\tLB1-LA1\n25\t0\tLB1-LA1\n'));
local_write(fullfile(d, 'sub-01_electrodes.tsv'), ...
    sprintf('name\tx\ty\tz\tsize\nLA1\t1\t0\t0\t1\nLA2\t2\t0\t0\t1\nLB1\t3\t0\t0\t1\n'));
info = ieeglab_import(sess, 'info', true);
tc.verifyEqual(info.labels, {'LA1','LA2','LB1','EKG'});
tc.verifyEqual(info.pnts, n);
tc.verifyEqual(info.sites, {'LA1-LA2','LA1-LB1'}, 'A site is named whatever the polarity.');
[~, E] = evalc('ieeglab_import(sess)');
tc.verifyEqual({E.chanlocs.labels}, {'LA1','LA2','LB1'});
tc.verifyEqual(double(E.data), data([4 2 1], :), 'AbsTol', 1e-3);
[~, E] = evalc('ieeglab_import(sess, ''sites'', {''LB1-LA1''}, ''pad'', [2 2])');
tc.verifyEqual(E.etc.ieeglab_import.first_sample, 18 * fs);
tc.verifyEqual(E.pnts, 9 * fs);
tc.verifyEqual(double(E.data), data([4 2 1], 18 * fs + (1:9 * fs)), 'AbsTol', 1e-3);
[~, L] = evalc('ieeglab_load(E, struct(''event_field'', ''electrical_stimulation_site''))');
tc.verifyEqual([L.event.latency], [2 7] * fs + 1);
tc.verifyEqual([L.chanlocs.X], [1 2 3]);
end

function local_write(f, txt)
fid = fopen(f, 'w'); fwrite(fid, txt, 'char'); fclose(fid);
end
