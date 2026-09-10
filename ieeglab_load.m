function [EEG, com] = ieeglab_load(EEG, opt)
% ieeglab_load() - Load BIDS electrode coordinates, events and channel status
%                  into an EEGLAB iEEG dataset, and select channels and events.
%
% Usage:
%   [EEG, com] = ieeglab_load(EEG)          % two dialogs
%   [EEG, com] = ieeglab_load(EEG, opt)     % headless, no dialogs
%
% Options (headless path; the dialogs fill the same fields):
%   .elec_tsv        *_electrodes.tsv path, or 'auto' (default) to find the BIDS
%                    sibling of the dataset file; '' to skip
%   .events_tsv      *_events.tsv path, 'auto' (default), or ''
%   .channels_tsv    *_channels.tsv path, 'auto' (default), or ''. Its 'status'
%                    column (good/bad) marks clinician/annotator bad channels.
%                    They are MARKED here and removed during preprocessing if
%                    requested, so the decision stays visible and reversible.
%   .event_field     events.tsv column naming the condition or stimulation site.
%                    Default: electrical_stimulation_site if present, else
%                    trial_type, else type.
%   .event_to_keep   values of event_field to keep. Default: all.
%   .chan_list       channel labels to keep. Default: all.
%   .drop_bad_events drop events whose BIDS 'status' column says 'bad'. Default true.
%
% Cancelling either dialog returns the dataset unmodified with com empty,
% which tells EEGLAB not to store a half-processed dataset.
%
% Example (headless):
%   EEG = ieeglab_load(EEG, struct('event_field','electrical_stimulation_site'));
%
% Cedric Cannard, iEEGLAB, 2025-2026

com = '';
EEG_in = EEG;
interactive = (nargin < 2 || isempty(opt));
if interactive, opt = struct(); end

if ndims(EEG.data) ~= 3 && (~isfield(EEG,'trials') || isempty(EEG.trials) || EEG.trials == 0)
    EEG.trials = 1;
end

% ---------- which files ----------
auto = struct('elec_tsv', local_bids_sibling(EEG, 'electrodes'), ...
              'events_tsv', local_bids_sibling(EEG, 'events'), ...
              'channels_tsv', local_bids_sibling(EEG, 'channels'));
if interactive
    [g, wasCanceled] = ieeglab_gui_load1(EEG.filepath, auto);
    if wasCanceled || isempty(g), return; end
    opt.elec_tsv = g.elec_tsv; opt.events_tsv = g.events_tsv; opt.channels_tsv = g.channels_tsv;
else
    for f = {'elec_tsv','events_tsv','channels_tsv'}
        if ~isfield(opt, f{1}) || strcmpi(opt.(f{1}), 'auto')
            opt.(f{1}) = auto.(f{1});
        end
    end
end
if ~isfield(opt,'drop_bad_events') || isempty(opt.drop_bad_events), opt.drop_bad_events = true; end

% ---------- electrodes ----------
elecs = [];
if ~isempty(opt.elec_tsv)
    elecs = readtable(opt.elec_tsv, 'FileType', 'text', 'Delimiter', '\t');
    % Properties.VariableNames, not fieldnames: fieldnames on a table also
    % returns Properties/Row/Variables.
    vn = elecs.Properties.VariableNames;
    missingCols = setdiff({'x','y','z'}, lower(vn));
    if ~isempty(missingCols)
        warning('ieeglab_load:noCoordinateColumns', ...
            ['The electrodes file has no %s column(s), so 3D coordinates cannot be loaded.\n' ...
             'Columns found: %s\nContinuing without electrode coordinates.'], ...
             strjoin(upper(missingCols), '/'), strjoin(vn, ', '));
        elecs = [];
    else
        xyz = elecs{:, ismember(lower(vn), {'x','y','z'})};
        if isempty(xyz) || all(~isfinite(xyz(:)))
            warning('ieeglab_load:emptyCoordinates', ...
                'The X/Y/Z columns of the electrodes file are empty or all non-finite.');
            elecs = [];
        end
    end
end

% ---------- events ----------
% events_from_tsv records the source: only a BIDS table has an 'onset' column
% in seconds, so only that case is re-injected into EEG.event.
opt.events_from_tsv = false;
if ~isempty(opt.events_tsv)
    events = readtable(opt.events_tsv, 'FileType', 'text', 'Delimiter', '\t');
    fprintf('%d events read from %s\n', height(events), local_short(opt.events_tsv));
    if ~ismember('onset', events.Properties.VariableNames)
        error('ieeglab_load:noOnset', ...
            ['The events file has no "onset" column, which BIDS requires and which is ' ...
             'needed to place events in the recording.\nColumns found: %s'], ...
             strjoin(events.Properties.VariableNames, ', '));
    end
    late = events.onset > EEG.xmax;
    if any(late)
        warning('ieeglab_load:eventsBeyondData', ...
            'Removing %d events whose onset is after the end of the recording (%.1f s).', nnz(late), EEG.xmax);
        events(late,:) = [];
    end
    % Clinician/annotator trial rejection: BIDS events 'status' good/bad
    if opt.drop_bad_events && ismember('status', events.Properties.VariableNames)
        st = lower(strtrim(string(events.status)));
        bad = st == "bad";
        if any(bad)
            fprintf('Dropping %d/%d events marked status=bad in the events file.\n', nnz(bad), numel(bad));
            events(bad,:) = [];
        end
    end
    opt.events = events;
    opt.events_from_tsv = true;
elseif isfield(EEG,'event') && ~isempty(EEG.event)
    fprintf('No events file; using the %d events already in the dataset.\n', numel(EEG.event));
    opt.events = struct2table(EEG.event, 'AsArray', true);
else
    fprintf('No events: continuous mode.\n');
end

% ---------- coordinates ----------
if ~isempty(elecs)
    EEG = get_elec_coor(EEG, elecs);
end
if EEG.nbchan == 0
    error('ieeglab_load:noChannels', 'No channels left in the dataset.');
end
if ~isfield(EEG.chanlocs, 'X') || isempty([EEG.chanlocs.X])
    error('ieeglab_load:noCoordinates', ...
        ['No electrode coordinates. Provide an electrodes.tsv with x, y and z columns ' ...
         '(BIDS: sub-XX_ses-YY_electrodes.tsv).']);
end

% ---------- clinician bad-channel labels: mark only ----------
if ~isempty(opt.channels_tsv) || ~isempty(elecs)
    [EEG, badT] = ieeglab_bad_channels(EEG, struct('channels_tsv', local_or(opt.channels_tsv, ''), ...
        'elec_tsv', local_or_table(elecs), 'action', 'mark', 'drop_bad_stim_sites', false, ...
        'auto_detect', false, 'verbose', true));
    opt.bad_labels  = cellstr(badT.label(badT.status == "bad"));
    opt.bad_reasons = cellstr(badT.reason(badT.status == "bad"));
end
opt.elec_labels = {EEG.chanlocs.labels};

% ---------- channel and event selection ----------
if interactive
    [opt, wasCanceled] = ieeglab_gui_load2(opt);
    if wasCanceled || isempty(opt)
        EEG = EEG_in;          % get_elec_coor already rewrote chanlocs
        return
    end
else
    if isfield(opt,'events') && istable(opt.events) && ~isempty(opt.events) && opt.events_from_tsv ...
            && (~isfield(opt,'event_field') || isempty(opt.event_field))
        opt.event_field = local_default_event_field(opt.events);
    end
end

% ---------- inject TSV events ----------
if isfield(opt,'events') && istable(opt.events) && ~isempty(opt.events) && opt.events_from_tsv
    col = opt.event_field;
    if ~ismember(col, opt.events.Properties.VariableNames)
        error('ieeglab_load:badEventField', ...
            'Event field "%s" is not a column of the events file. Available: %s', ...
            col, strjoin(opt.events.Properties.VariableNames, ', '));
    end
    ty = local_to_str(opt.events.(col));
    missing = ty == "" | lower(ty) == "n/a";
    if any(missing)
        fprintf('Dropping %d/%d events with no value in "%s".\n', nnz(missing), numel(missing), col);
        opt.events(missing,:) = [];
        ty(missing) = [];
    end
    lat = opt.events.onset;
    hasDur = ismember('duration', opt.events.Properties.VariableNames);

    % Replace existing events wholesale so no stale event or field survives.
    EEG.event = [];
    EEG.urevent = [];
    for iEv = 1:numel(ty)
        EEG.event(iEv).type    = char(ty(iEv));
        EEG.event(iEv).latency = lat(iEv) * EEG.srate + 1;     % samples, 1-based
        if hasDur
            d = opt.events.duration(iEv);
            if iscell(d), d = str2double(d); end
            if isnumeric(d) && isfinite(d), EEG.event(iEv).duration = d * EEG.srate; end
        end
    end
    fprintf('%d events placed using column "%s".\n', numel(ty), col);
end
EEG = eeg_checkset(EEG, 'eventconsistency');
EEG = eeg_checkset(EEG);

% ---------- keep only events of interest ----------
if isfield(opt,'event_to_keep') && ~isempty(opt.event_to_keep) && isfield(EEG,'event') && ~isempty(EEG.event)
    keep = ismember(string({EEG.event.type}), string(opt.event_to_keep));
    fprintf('Keeping %d/%d events of interest (%s).\n', nnz(keep), numel(keep), ...
        strjoin(cellstr(string(opt.event_to_keep)), ', '));
    if ~any(keep)
        error('ieeglab_load:noEventsKept', ...
            'None of the requested event types (%s) are present. Types available: %s', ...
            strjoin(cellstr(string(opt.event_to_keep)), ', '), strjoin(unique({EEG.event.type}), ', '));
    end
    EEG.event(~keep) = [];
    if isfield(opt,'events') && istable(opt.events) && height(opt.events) == numel(keep)
        opt.events(~keep,:) = [];
    end
    EEG = eeg_checkset(EEG, 'eventconsistency');
end

% ---------- keep only channels of interest ----------
if isfield(opt,'chan_list') && ~isempty(opt.chan_list)
    want = cellstr(opt.chan_list);
    allLabels = {EEG.chanlocs.labels};
    unknown = setdiff(want, allLabels);
    if ~isempty(unknown)
        warning('ieeglab_load:unknownChannels', 'Ignoring channel names not in the dataset: %s', strjoin(unknown, ', '));
    end
    removed = setdiff(allLabels, want, 'stable');
    if ~isempty(removed) && numel(removed) < numel(allLabels)
        fprintf('Removing %d/%d channels not selected: %s\n', numel(removed), numel(allLabels), strjoin(removed, ', '));
        EEG = pop_select(EEG, 'nochannel', removed);
        % CCEP: drop trials whose stimulated pair used a removed contact. Exact
        % token match - contains() also matched 'RA10-RA9' when removing 'RA1'.
        if isfield(EEG,'event') && ~isempty(EEG.event)
            drop = ieeglab_events_using(EEG.event, removed);
            if any(drop)
                fprintf('Removing %d/%d stimulation events whose pair uses a removed contact.\n', nnz(drop), numel(drop));
                EEG.event(drop) = [];
                if isfield(opt,'events') && istable(opt.events) && height(opt.events) == numel(drop)
                    opt.events(drop,:) = [];
                end
                EEG = eeg_checkset(EEG, 'eventconsistency');
            end
        end
    elseif numel(removed) == numel(allLabels)
        error('ieeglab_load:noChannelsSelected', 'None of the requested channels are in the dataset.');
    end
end

if ~isfield(EEG,'ieeglab') || isempty(EEG.ieeglab), EEG.ieeglab = struct(); end
if isfield(EEG.ieeglab,'opt') && isstruct(EEG.ieeglab.opt)
    base = EEG.ieeglab.opt; fns = fieldnames(opt);
    for k = 1:numel(fns), base.(fns{k}) = opt.(fns{k}); end
    opt = base;
end
EEG.ieeglab.opt = opt;
EEG = eeg_checkset(EEG);

com = sprintf('EEG = ieeglab_load(EEG, %s);', local_literal(opt));
end

% ======================= helpers =======================

function p = local_bids_sibling(EEG, suffix)
% BIDS sidecar next to the dataset file. electrodes.tsv is per session (no
% task/run entities); events and channels are per run.
p = '';
if ~isfield(EEG,'filepath') || isempty(EEG.filepath) || ~isfolder(EEG.filepath), return; end
d = EEG.filepath;
stem = '';
if isfield(EEG,'filename') && ~isempty(EEG.filename)
    stem = regexprep(EEG.filename, '_(ieeg|eeg)\.[^.]+$', '');
end
cands = {};
if ~isempty(stem)
    cands{end+1} = fullfile(d, [stem '_' suffix '.tsv']);
    ses = regexp(stem, '^(sub-[^_]+(_ses-[^_]+)?)', 'tokens', 'once');
    if ~isempty(ses), cands{end+1} = fullfile(d, [ses{1} '_' suffix '.tsv']); end
end
for i = 1:numel(cands)
    if exist(cands{i}, 'file') == 2, p = cands{i}; return; end
end
g = dir(fullfile(d, ['*_' suffix '.tsv']));
if numel(g) == 1, p = fullfile(g.folder, g.name); end
end

function f = local_default_event_field(T)
vn = T.Properties.VariableNames;
for c = {'electrical_stimulation_site','trial_type','type','value'}
    if ismember(c{1}, vn), f = c{1}; return; end
end
rest = setdiff(vn, {'onset','duration','sample','sample_start'}, 'stable');
f = rest{1};
end

function s = local_to_str(x)
if iscell(x)
    s = strings(numel(x),1);
    for i = 1:numel(x)
        v = x{i};
        if isnumeric(v) || islogical(v), v = num2str(v); end
        s(i) = strtrim(string(v));
    end
elseif isnumeric(x) || islogical(x)
    s = strings(numel(x),1);
    ok = ~isnan(double(x));
    s(ok) = string(x(ok));
else
    s = strtrim(string(x));
    s(ismissing(s)) = "";
end
s = s(:);
end

function v = local_or(x, d)
if isempty(x), v = d; else, v = x; end
end

function v = local_or_table(t)
if isempty(t), v = ''; else, v = t; end
end

function s = local_short(p)
[~, n, e] = fileparts(char(p)); s = [n e];
end

function s = local_literal(opt)
% Re-runnable struct literal for the EEGLAB history (file choices and selections).
keys = {'elec_tsv','events_tsv','channels_tsv','event_field','event_to_keep','chan_list','drop_bad_events'};
parts = {};
for k = 1:numel(keys)
    if ~isfield(opt, keys{k}), continue; end
    v = opt.(keys{k});
    if ischar(v) || isstring(v)
        parts{end+1} = sprintf('''%s'',''%s''', keys{k}, strrep(char(v), '''', '''''')); %#ok<AGROW>
    elseif iscell(v) && ~isempty(v)
        q = cellfun(@(x) ['''' strrep(char(string(x)), '''', '''''') ''''], v, 'UniformOutput', false);
        parts{end+1} = sprintf('''%s'',{{%s}}', keys{k}, strjoin(q, ',')); %#ok<AGROW>
    elseif islogical(v) || isnumeric(v)
        parts{end+1} = sprintf('''%s'',%s', keys{k}, mat2str(v)); %#ok<AGROW>
    end
end
s = ['struct(' strjoin(parts, ', ') ')'];
end
