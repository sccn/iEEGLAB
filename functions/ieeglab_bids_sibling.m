function p = ieeglab_bids_sibling(EEG, suffix)
% ieeglab_bids_sibling() - The BIDS sidecar that belongs to this dataset file.
%
% Usage:
%   p = ieeglab_bids_sibling(EEG, 'electrodes')   % per session
%   p = ieeglab_bids_sibling(EEG, 'events')       % per run
%   p = ieeglab_bids_sibling(EEG, 'channels')     % per run
%
% Returns '' when nothing suitable is found.
%
% A sidecar is accepted only when its BIDS entities are a subset of the
% dataset's own (the BIDS inheritance rule): sub-02_ses-01_electrodes.tsv
% applies to sub-02_ses-01_task-ccep_run-01_ieeg.set, but a run-02 events or
% channels file never applies to run 01. The 'space' entity of electrodes
% files is ignored for matching, since it describes coordinates, not data.
% Earlier versions fell back to "any single file with this suffix", which
% injected one run's events - or one run's clinician bad-channel labels - into
% another run.
%
% When the dataset filename carries no BIDS entities at all, a single
% candidate in the folder is accepted, as before.
%
% Cedric Cannard, iEEGLAB, 2026

p = '';
% Where to look, in order: next to the dataset file; next to the raw BIDS file it
% was imported from, when EEGLAB's BIDS import recorded it (EEG.BIDS.sourcefile);
% and, for datasets saved by that import in its default output folder
% <root>/derivatives/<pipeline>/sub-..., the matching folder of the raw dataset.
% Without the last two, nothing is found after File > Import data > BIDS, since
% the imported .set lives in derivatives and the sidecars do not.
places = {};
if isfield(EEG,'filepath') && ~isempty(EEG.filepath)
    fname = '';
    if isfield(EEG,'filename'), fname = char(EEG.filename); end
    places(end+1,:) = {char(EEG.filepath), fname};
end
if isfield(EEG,'BIDS') && isstruct(EEG.BIDS) && isfield(EEG.BIDS,'sourcefile') && ~isempty(EEG.BIDS.sourcefile)
    [sd, sn, se] = fileparts(char(EEG.BIDS.sourcefile));
    places(end+1,:) = {sd, [sn se]};
end
if ~isempty(places)
    tok = regexp(places{1,1}, '^(.*)[\\/]derivatives[\\/][^\\/]+[\\/](sub-[^\\/]+.*)$', 'tokens', 'once');
    if ~isempty(tok), places(end+1,:) = {fullfile(tok{1}, tok{2}), places{1,2}}; end
end
for i = 1:size(places, 1)
    if isfolder(places{i,1})
        p = local_find(places{i,1}, places{i,2}, suffix);
        if ~isempty(p), return; end
    end
end
end

function p = local_find(d, fname, suffix)
p = '';
stem = '';
if ~isempty(fname)
    stem = regexprep(fname, '_(ieeg|eeg)\.[^.]+$', '');
    stem = regexprep(stem, '\.[^.]+$', '');
end
ents = local_entities(stem);

if ~isempty(stem)
    c = fullfile(d, [stem '_' suffix '.tsv']);
    if exist(c, 'file') == 2, p = c; return; end
end

g = dir(fullfile(d, ['*_' suffix '.tsv']));
if isempty(g), return; end

if isempty(fieldnames(ents))
    if numel(g) == 1, p = fullfile(g.folder, g.name); end
    return
end

ok = {}; nEnt = [];
for i = 1:numel(g)
    ce = local_entities(regexprep(g(i).name, ['_' suffix '\.tsv$'], ''));
    if isfield(ce, 'space'), ce = rmfield(ce, 'space'); end
    if local_is_subset(ce, ents)
        ok{end+1} = fullfile(g(i).folder, g(i).name); %#ok<AGROW>
        nEnt(end+1) = numel(fieldnames(ce)); %#ok<AGROW>
    end
end
if isempty(ok)
    warning('ieeglab_bids_sibling:noMatch', ...
        ['Found %d *_%s.tsv file(s) in %s, but none belongs to %s (BIDS entities differ, ' ...
         'e.g. another run). Not using any; pass the file explicitly if it is right.'], ...
        numel(g), suffix, d, fname);
    return
end
% Most specific match wins (run-level beats session-level); a tie is ambiguous.
best = find(nEnt == max(nEnt));
if numel(best) == 1
    p = ok{best};
else
    warning('ieeglab_bids_sibling:ambiguous', ...
        'Several *_%s.tsv files match %s equally well (%s). Pass the one to use explicitly.', ...
        suffix, fname, strjoin(cellfun(@local_name, ok(best), 'UniformOutput', false), ', '));
end
end

function e = local_entities(stem)
e = struct();
if isempty(stem), return; end
parts = strsplit(stem, '_');
for i = 1:numel(parts)
    kv = regexp(parts{i}, '^([a-zA-Z]+)-(.+)$', 'tokens', 'once');
    if ~isempty(kv) && isvarname(kv{1})
        e.(lower(kv{1})) = kv{2};
    end
end
end

function tf = local_is_subset(a, b)
fa = fieldnames(a);
tf = true;
for i = 1:numel(fa)
    if ~isfield(b, fa{i}) || ~strcmp(a.(fa{i}), b.(fa{i}))
        tf = false; return
    end
end
end

function n = local_name(p)
[~, a, b] = fileparts(p); n = [a b];
end
