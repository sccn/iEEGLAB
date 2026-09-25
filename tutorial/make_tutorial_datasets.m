function out = make_tutorial_datasets(outRoot, src)
% make_tutorial_datasets() - Build the native-rate tutorial datasets from OpenNeuro.
%
% Usage:
%   make_tutorial_datasets(outRoot, src)
%
%   outRoot  folder that receives one small BIDS dataset per tutorial dataset
%   src      struct of local OpenNeuro copies, fields (either may be omitted):
%              .ds004696   root of ds004696 (sEEG, MEF3), e.g. the NEMAR download
%              .ds004080   root of ds004080 (ECoG, BrainVision)
%
% Each output is a valid BIDS dataset with one subject and one run: the recording
% at its native sampling rate, cut to a few contacts and a few stimulation sites,
% saved as an EEGLAB .set file, next to the original sidecars restricted to what
% was kept (channels, electrodes, events with onsets shifted to the new start,
% ieeg.json with the new duration). Nothing is filtered, re-referenced or
% resampled. The .set files carry no events and no coordinates on purpose: the
% tutorial loads them from the BIDS sidecars.
%
% The selection (contacts, sites, padding) is set in local_config below. Each .set
% stays near 30 MB so the datasets can live in the GitHub repository (GitHub
% warns above 50 MB per file and refuses 100 MB). Sites are
% stimulated in consecutive blocks in both datasets, so the kept stretch is one
% continuous window: no boundary events, filters behave as on the full file.
%
% Requires EEGLAB with the MEF3 plugin (ds004696) and bva-io (ds004080).
%
% Cedric Cannard, iEEGLAB, 2026

if nargin < 2, src = struct(); end
if ~exist(outRoot, 'dir'), mkdir(outRoot); end
cfgs = local_config();
out = struct([]);
for c = cfgs
    if ~isfield(src, c.dataset) || isempty(src.(c.dataset))
        fprintf('Skipping %s: no source folder given.\n', c.dataset);
        continue
    end
    r = local_make_one(c, char(src.(c.dataset)), outRoot);
    out = [out r]; %#ok<AGROW>
end
end

% -------------------------------------------------------------------------
function cfgs = local_config()
% sEEG: ds004696 sub-02 (Ojeda Valencia et al., 2023), the subject of the earlier
% tutorial. RA2-RA3, RA3-RA4 and RA4-RA5 are three of the 21 sites in the
% published CRP results (derivatives/stats), stimulated back to back: 34 pulses.
% The whole RA lead plus RB13-RB15 gives 18 contacts, 3 of them clinician-bad
% (RA15, RB14, RB15, which are also the three without coordinates).
cfgs(1) = struct( ...
    'dataset',  'ds004696', ...
    'name',     'iEEGLAB tutorial, sEEG: extract of ds004696 sub-02', ...
    'outdir',   'ieeglab_tutorial_seeg', ...
    'doi',      'doi:10.18112/openneuro.ds004696.v1.0.1', ...
    'sub', 'sub-02', 'ses', 'ses-ieeg01', 'task', 'ccep', 'run', 'run-01', ...
    'format',   'mef', ...
    'contacts', {{'^RA\d+$', '^RB1[3-5]$'}}, ...
    'sites',    {{'RA2-RA3','RA3-RA4','RA4-RA5'}}, ...
    'pad',      [5 4], ...
    'surfaces', {{'derivatives/freesurfer/sub-02/pial.L.surf.gii', ...
                  'derivatives/freesurfer/sub-02/pial.R.surf.gii'}});
% ECoG: ds004080 (van Blooijs et al., 2023), sub-ccepAgeUMCU02, a 2048 Hz run with
% 16 sites of 10 pulses. The 2x8 FB grid holds the three sites kept here (30
% pulses); F05-F08 of the 2x4 F grid add 2 clinician-bad contacts (F07, F08).
% Coordinates are in fsaverage space (no individual MRI is shared).
cfgs(2) = struct( ...
    'dataset',  'ds004080', ...
    'name',     'iEEGLAB tutorial, ECoG: extract of ds004080 sub-ccepAgeUMCU02', ...
    'outdir',   'ieeglab_tutorial_ecog', ...
    'doi',      'doi:10.18112/openneuro.ds004080', ...
    'sub', 'sub-ccepAgeUMCU02', 'ses', 'ses-1', 'task', 'SPESclin', 'run', 'run-041456', ...
    'format',   'brainvision', ...
    'contacts', {{'^FB\d+$', '^F0[5-8]$'}}, ...
    'sites',    {{'FB54-FB62','FB54-FB53','FB54-FB55'}}, ...
    'pad',      [5 4], ...
    'surfaces', {{}});
end

% -------------------------------------------------------------------------
function r = local_make_one(c, srcRoot, outRoot)
ieegDir = fullfile(srcRoot, c.sub, c.ses, 'ieeg');
stem    = sprintf('%s_%s_task-%s_%s', c.sub, c.ses, c.task, c.run);
sesStem = sprintf('%s_%s', c.sub, c.ses);
dstRoot = fullfile(outRoot, c.outdir);
dstIeeg = fullfile(dstRoot, c.sub, c.ses, 'ieeg');
if ~exist(dstIeeg, 'dir'), mkdir(dstIeeg); end
info = jsondecode(fileread(fullfile(ieegDir, [stem '_ieeg.json'])));
fs = info.SamplingFrequency;

% ---- channels matching the patterns, in file order
[chH, chR] = local_tsv_read(fullfile(ieegDir, [stem '_channels.tsv']));
names = chR(:, strcmp(chH, 'name'));
keepCh = false(size(names));
for p = c.contacts, keepCh = keepCh | ~cellfun(@isempty, regexp(names, p{1}, 'once')); end
if ~any(keepCh), error('make_tutorial_datasets:noChannels', 'No channel matches %s.', strjoin(c.contacts, ' ')); end
names = names(keepCh);

% ---- window around the chosen sites, on exact samples
[evH, evR] = local_tsv_read(fullfile(ieegDir, [stem '_events.tsv']));
onset = str2double(evR(:, strcmp(evH, 'onset')));
site  = evR(:, strcmp(evH, 'electrical_stimulation_site'));
isSite = ismember(site, c.sites);
missing = setdiff(c.sites, site(isSite));
if ~isempty(missing), error('make_tutorial_datasets:noSite', 'Site(s) not in events.tsv: %s', strjoin(missing, ', ')); end
s0 = floor((min(onset(isSite)) - c.pad(1)) * fs);   % first kept sample, 0-based
s1 = ceil((max(onset(isSite)) + c.pad(2)) * fs);    % one past the last kept sample
t0 = s0 / fs;

% ---- data
switch c.format
    case 'mef'
        f = fullfile(ieegDir, [stem '_ieeg.mefd']);
        [~, X] = ieeglab_load_mefd(f, [], names, 'samples', [s0 s1]);
    case 'brainvision'
        hdr = pop_loadbv(ieegDir, [stem '_ieeg.vhdr'], [], [], true);
        [tf, idx] = ismember(names, {hdr.chanlocs.labels});
        if ~all(tf), error('make_tutorial_datasets:noChannels', 'Not in the .vhdr: %s', strjoin(names(~tf), ', ')); end
        E = pop_loadbv(ieegDir, [stem '_ieeg.vhdr'], [s0+1 s1], idx(:)');
        X = double(E.data);
    otherwise
        error('make_tutorial_datasets:format', 'Unknown format %s.', c.format);
end
if size(X, 2) ~= s1 - s0, error('make_tutorial_datasets:short', 'Read %d samples, expected %d.', size(X,2), s1 - s0); end
if any(~isfinite(X(:))), error('make_tutorial_datasets:nan', 'The extract contains NaN or Inf.'); end

EEG = eeg_emptyset;
EEG.setname  = c.name;
EEG.data     = single(X);
EEG.srate    = fs;
EEG.nbchan   = size(X, 1);
EEG.pnts     = size(X, 2);
EEG.trials   = 1;
EEG.xmin     = 0;
chType = chR(keepCh, strcmp(chH, 'type'));
EEG.chanlocs = struct('labels', names(:)', 'type', chType(:)');
EEG.ref      = 'intracranial';
EEG.etc.source = struct('dataset', c.dataset, 'doi', c.doi, 'file', [stem '_ieeg'], ...
    'first_sample', s0, 'last_sample', s1 - 1, 'onset_shift_s', t0);
EEG = eeg_checkset(EEG);
pop_saveset(EEG, 'filename', [stem '_ieeg.set'], 'filepath', dstIeeg, 'savemode', 'onefile');

% ---- sidecars restricted to what was kept
local_tsv_write(fullfile(dstIeeg, [stem '_channels.tsv']), chH, chR(keepCh, :));
[elH, elR] = local_tsv_read(fullfile(ieegDir, [sesStem '_electrodes.tsv']));
local_tsv_write(fullfile(dstIeeg, [sesStem '_electrodes.tsv']), elH, elR(ismember(elR(:,1), names), :));
local_copy(fullfile(ieegDir, [sesStem '_electrodes.json']), dstIeeg);
local_copy(fullfile(ieegDir, [sesStem '_coordsystem.json']), dstIeeg);
local_copy(fullfile(ieegDir, [stem '_channels.json']), dstIeeg);
local_copy(fullfile(ieegDir, [stem '_events.json']), dstIeeg);

% events: rows whose onset falls in the window, times shifted to the new start
inWin = onset >= t0 & onset < s1 / fs;
evR = evR(inWin, :);
for col = {'onset', 'offset'}
    k = strcmp(evH, col{1});
    if any(k), evR(:, k) = local_shift(evR(:, k), -t0, '%.15g'); end
end
for col = {'sample_start', 'sample_end'}
    k = strcmp(evH, col{1});
    if any(k), evR(:, k) = local_shift(evR(:, k), -s0, '%d'); end
end
local_tsv_write(fullfile(dstIeeg, [stem '_events.tsv']), evH, evR);

info.RecordingDuration = (s1 - s0) / fs;
countField = struct('SEEG', 'SEEGChannelCount', 'ECOG', 'ECOGChannelCount', 'EEG', 'EEGChannelCount', ...
    'ECG', 'ECGChannelCount', 'EMG', 'EMGChannelCount', 'EOG', 'EOGChannelCount', ...
    'MISC', 'MiscChannelCount', 'TRIG', 'TriggerChannelCount');
for t = fieldnames(countField)'
    if isfield(info, countField.(t{1})), info.(countField.(t{1})) = nnz(strcmpi(chType, t{1})); end
end
local_json_write(fullfile(dstIeeg, [stem '_ieeg.json']), info);

% ---- dataset level files
dd = jsondecode(fileread(fullfile(srcRoot, 'dataset_description.json')));
desc = struct('Name', c.name, 'BIDSVersion', '1.9.0', 'DatasetType', 'raw', 'License', 'CC0');
if isfield(dd, 'Authors'), desc.Authors = dd.Authors; end
desc.HowToAcknowledge = sprintf('Please cite the source dataset (%s) and its paper.', c.doi);
desc.SourceDatasets = {struct('DOI', c.doi, 'URL', ['https://openneuro.org/datasets/' c.dataset])};
desc.GeneratedBy = {struct('Name', 'iEEGLAB make_tutorial_datasets', 'CodeURL', 'https://github.com/sccn/iEEGLAB')};
local_json_write(fullfile(dstRoot, 'dataset_description.json'), desc);
[pH, pR] = local_tsv_read(fullfile(srcRoot, 'participants.tsv'));
local_tsv_write(fullfile(dstRoot, 'participants.tsv'), pH, pR(strcmp(pR(:,1), c.sub), :));
nSite = arrayfun(@(k) nnz(strcmp(site(inWin), c.sites{k})), 1:numel(c.sites));
readme = sprintf(['# %s\n\nExtract of OpenNeuro %s (%s), made by iEEGLAB''s tutorial/make_tutorial_datasets.m ' ...
    'for the iEEGLAB tutorial. Native sampling rate (%g Hz), not filtered, not re-referenced.\n\n' ...
    '- Recording: %s, samples %d to %d of the original (0-based), %.1f s\n' ...
    '- Contacts: %s (%d)\n- Stimulation sites: %s (%s pulses)\n\n' ...
    'Event onsets are shifted by -%.6f s relative to the original events.tsv. ' ...
    'Cite the original dataset and paper when using these data.\n'], ...
    c.name, c.dataset, c.doi, fs, stem, s0, s1 - 1, (s1 - s0) / fs, strjoin(names', ' '), numel(names), ...
    strjoin(c.sites, ', '), strjoin(arrayfun(@num2str, nSite, 'uni', 0), '/'), t0);
local_text_write(fullfile(dstRoot, 'README.md'), readme);
for s = c.surfaces
    d = fullfile(dstRoot, fileparts(s{1}));
    if ~exist(d, 'dir'), mkdir(d); end
    copyfile(fullfile(srcRoot, s{1}), d);
end

setFile = dir(fullfile(dstIeeg, [stem '_ieeg.set']));
r = struct('dataset', c.dataset, 'folder', dstRoot, 'contacts', numel(names), 'seconds', (s1 - s0) / fs, ...
    'pulses', nnz(isSite & inWin), 'mb', setFile.bytes / 1e6);
fprintf('%s: %d contacts, %.1f s, %d pulses at %s, %.0f MB -> %s\n', c.dataset, r.contacts, r.seconds, ...
    r.pulses, strjoin(c.sites, ' '), r.mb, dstRoot);
end

% -------------------------------------------------------------------------
% TSV as text, so 'n/a' and number formats survive unchanged
function [h, rows] = local_tsv_read(f)
L = splitlines(strtrim(fileread(f)));
L = L(~cellfun(@isempty, L));
h = strsplit(L{1}, '\t');
rows = cell(numel(L) - 1, numel(h));
for i = 2:numel(L)
    v = strsplit(L{i}, '\t', 'CollapseDelimiters', false);
    v(end+1:numel(h)) = {'n/a'};
    rows(i-1, :) = v(1:numel(h));
end
end

function local_tsv_write(f, h, rows)
L = [{strjoin(h, '\t')}; arrayfun(@(i) strjoin(rows(i,:), '\t'), (1:size(rows,1))', 'uni', 0)];
local_text_write(f, [strjoin(L, newline) newline]);
end

function v = local_shift(v, d, fmt)
for i = 1:numel(v)
    x = str2double(v{i});
    if isfinite(x), v{i} = sprintf(fmt, x + d); end
end
end

function local_json_write(f, s)
local_text_write(f, jsonencode(s, 'PrettyPrint', true));
end

function local_text_write(f, txt)
fid = fopen(f, 'w');
if fid < 0, error('make_tutorial_datasets:write', 'Cannot write %s.', f); end
fwrite(fid, txt, 'char');
fclose(fid);
end

function local_copy(f, d)
if exist(f, 'file'), copyfile(f, d); end
end
