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
if ~isfield(EEG,'filepath') || isempty(EEG.filepath) || ~isfolder(char(EEG.filepath)), return; end
d = char(EEG.filepath);
stem = '';
if isfield(EEG,'filename') && ~isempty(EEG.filename)
    stem = regexprep(char(EEG.filename), '_(ieeg|eeg)\.[^.]+$', '');
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
        numel(g), suffix, d, char(EEG.filename));
    return
end
% Most specific match wins (run-level beats session-level); a tie is ambiguous.
best = find(nEnt == max(nEnt));
if numel(best) == 1
    p = ok{best};
else
    warning('ieeglab_bids_sibling:ambiguous', ...
        'Several *_%s.tsv files match %s equally well (%s). Pass the one to use explicitly.', ...
        suffix, char(EEG.filename), strjoin(cellfun(@local_name, ok(best), 'UniformOutput', false), ', '));
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
