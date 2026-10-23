function [EEG, com] = pop_ieeglab_import(file, varargin)
% pop_ieeglab_import() - Choose a raw iEEG recording and what to read from it,
%                        then read it with ieeglab_import.
%
% Usage:
%   [EEG, com] = pop_ieeglab_import                 % dialogs
%   [EEG, com] = pop_ieeglab_import(folder)         % dialogs, starting from this folder or recording
%   [EEG, com] = pop_ieeglab_import(file, 'sites', {'RA2-RA3'})   % as ieeglab_import
%
% The dialogs: (1) the folder of the dataset or of one participant, (2) the
% recording, when the folder holds several, (3) the channel types and,
% optionally, the stimulation sites to read. Reading fewer channels or sites
% takes proportionally less memory and time (see ieeglab_import).
%
% Cedric Cannard, iEEGLAB, 2026

EEG = []; com = '';
if nargin >= 2
    [EEG, com] = ieeglab_import(file, varargin{:});
    return
end

% ---------- 1. where the data are ----------
if nargin < 1 || isempty(file)
    bar = repmat('=', 1, 78);
    fprintf(['\n%s\n iEEGLAB > Import raw recording\n%s\n' ...
        ' SELECT, in the window that opens, the folder of the dataset as downloaded\n' ...
        '   (for example ds004696), or the folder of one participant (sub-02).\n' ...
        ' The recordings found in it (MEF3 .mefd, BrainVision .vhdr, EDF, BDF, .set) are\n' ...
        '   listed next; then you choose the channels and the stimulation sites to read.\n%s\n\n'], bar, bar, bar);
    file = uigetdir(pwd, 'SELECT THE DATASET FOLDER (e.g. ds004696) OR A PARTICIPANT FOLDER (sub-02)');
    if isequal(file, 0), return; end
end
file = regexprep(char(file), '[\\/]+$', '');

% ---------- 2. which recording ----------
[~, ~, ext] = fileparts(file);
if isfolder(file) && ~strcmpi(ext, '.mefd')
    fprintf('[import] Looking for recordings in %s ...\n', file);
    rec = local_find_recordings(file);
    if isempty(rec)
        errordlg(sprintf(['No iEEG recording found in\n%s\n\nExpected a BIDS dataset or participant folder holding ' ...
            '*_ieeg.mefd, *_ieeg.vhdr, *_ieeg.edf, *_ieeg.bdf or *_ieeg.set.'], file), 'iEEGLAB');
        return
    end
    k = 1;
    if numel(rec) > 1
        shown = cellfun(@(r) r(numel(file) + 2:end), rec, 'UniformOutput', false);
        [k, ok] = listdlg('ListString', shown, 'SelectionMode', 'single', 'Name', 'iEEGLAB: import raw recording', ...
            'PromptString', sprintf('%d recordings found. Which one?', numel(rec)), 'ListSize', [620 300]);
        if ~ok || isempty(k), return; end
    end
    file = rec{k};
end

% ---------- 3. what to read ----------
info = ieeglab_import(file, 'info', true);
[~, name, ext] = fileparts(file);
known = ~cellfun(@isempty, info.types);
uTypes = unique(info.types(known), 'stable');
nType  = cellfun(@(t) nnz(strcmp(info.types, t)), uTypes);
isBrain = ismember(uTypes, {'SEEG','ECOG','DBS'});
if ~any(isBrain), isBrain(:) = true; end
gb = @(nCh, nSamp) nCh * nSamp * 4 / 2^30;
nBrain = sum(nType(isBrain)); if isempty(uTypes), nBrain = numel(info.labels); end

typeList = arrayfun(@(i) sprintf('%s (%d channels)', uTypes{i}, nType(i)), 1:numel(uTypes), 'UniformOutput', false);
siteList = arrayfun(@(i) sprintf('%s (%d pulses)', info.sites{i}, info.site_counts(i)), 1:numel(info.sites), 'UniformOutput', false);
hasTypes = ~isempty(uTypes); hasSites = ~isempty(info.sites);
if ~hasTypes, typeList = {'(no channels.tsv next to the recording: all channels are read)'}; end
if ~hasSites, siteList = {'(no stimulation sites in an events.tsv next to the recording)'}; end
cbSites = 'set(findobj(gcbf, ''tag'', ''sites''), ''enable'', fastif(get(gcbo, ''value''), ''on'', ''off''));';

uilist = { ...
    {'style' 'text' 'string' [name ext] 'fontweight' 'bold'}, ...
    {'style' 'text' 'string' sprintf('%d channels, %.1f minutes at %g Hz', numel(info.labels), info.pnts / info.srate / 60, info.srate)}, ...
    {'style' 'text' 'string' sprintf('The whole recording with the brain channels takes %.1f GB of memory. Fewer channel types or sites take less.', gb(nBrain, info.pnts))}, ...
    {}, ...
    {'style' 'text' 'string' 'Channel types to read (Ctrl+click for several):'}, ...
    {'style' 'listbox' 'string' typeList 'tag' 'types' 'max' 2 'value' find(isBrain) 'enable' fastif(hasTypes, 'on', 'off')}, ...
    {'style' 'text' 'string' ['The brain channels are selected. Ctrl+click ECG, EMG, EOG or others to read them too (for ' ...
        'example for heart-brain analyses with BrainBeats): they are kept in the dataset and left out of the CCEP analysis.']}, ...
    {}, ...
    {'style' 'checkbox' 'string' 'Read only these stimulation sites (Ctrl+click for several):' 'tag' 'usesites' 'value' 0 ...
        'callback' cbSites 'enable' fastif(hasSites, 'on', 'off')}, ...
    {'style' 'listbox' 'string' siteList 'tag' 'sites' 'max' 2 'value' [] 'enable' 'off'}, ...
    {'style' 'text' 'string' 'Coordinates, events and bad channels are loaded in the next step: iEEGLAB > Load electrode coordinates and events.'} };
geom = {1 1 1 1 1 1 1 1 1 1 1};
gv   = [1 1 1 0.3 1 max(2, min(4, numel(typeList))) 2 0.3 1 max(3, min(8, numel(siteList))) 1];
[res, ~, ~, out] = inputgui('geometry', geom, 'geomvert', gv, 'uilist', uilist, ...
    'helpcom', 'pophelp(''ieeglab_import'')', 'title', 'iEEGLAB: import raw recording', 'minwidth', 640);
if isempty(res), return; end

args = {};
if hasTypes
    sel = out.types;
    if isempty(sel), errordlg('Select at least one channel type.', 'iEEGLAB'); return; end
    brain = uTypes(ismember(uTypes, {'SEEG','ECOG','DBS'}));
    if ~isempty(brain) && isequal(sort(uTypes(sel)), sort(brain))
        % the default of ieeglab_import: nothing to pass
    elseif numel(sel) == numel(uTypes)
        args = [args {'types', 'all'}];
    else
        args = [args {'types', uTypes(sel)}];
    end
end
if hasSites && out.usesites
    if isempty(out.sites), errordlg('Tick the sites to read, or untick "Read only these stimulation sites".', 'iEEGLAB'); return; end
    args = [args {'sites', info.sites(out.sites)}];
end
[EEG, com] = ieeglab_import(file, args{:}, 'waitbar', true);
end

function rec = local_find_recordings(root)
% Raw iEEG recordings below root, derivatives left out
rec = {};
for pat = {'*_ieeg.mefd', '*_ieeg.vhdr', '*_ieeg.edf', '*_ieeg.bdf', '*_ieeg.set'}
    d = dir(fullfile(root, '**', pat{1}));
    for i = 1:numel(d)
        f = fullfile(d(i).folder, d(i).name);
        if contains(f, [filesep 'derivatives' filesep]), continue; end
        if strcmp(pat{1}, '*_ieeg.mefd') ~= d(i).isdir, continue; end   % .mefd is a folder, the others are files
        rec{end+1} = f; %#ok<AGROW>
    end
end
rec = unique(rec, 'stable');
end
