function [EEG, com] = ieeglab_load(EEG, opt)
% ieeglab_load() - Load BIDS electrode coordinates, events and channel status
%                  into an EEGLAB iEEG dataset, and select channels and events.
%
% Usage:
%   [EEG, com] = ieeglab_load(EEG)          % two dialogs
%   [EEG, com] = ieeglab_load(EEG, opt)     % headless, no dialogs
%
% Options (headless path; the dialogs fill the same fields):
%   .elec_tsv        *_electrodes.tsv path, 'auto' (default: the BIDS sidecar of
%                    this dataset, see ieeglab_bids_sibling), or '' to skip
%   .events_tsv      *_events.tsv path, 'auto' (default), or ''
%   .channels_tsv    *_channels.tsv path, 'auto' (default), or ''. Channels with
%                    status 'bad' are MARKED here; whether they are removed is
%                    decided in preprocessing, so the decision stays visible.
%   .event_field     column naming the condition or stimulation site - of the
%                    events.tsv, or a field of EEG.event when the events are
%                    already in the dataset. Default: electrical_stimulation_site
%                    if present, else trial_type, else type.
%   .event_to_keep   values of event_field to keep. Default: all.
%   .chan_list       channel labels to keep. Default: all.
%   .drop_bad_events drop events whose BIDS 'status' is 'bad'. Default true.
%
% Provenance: every event placed from an events.tsv carries tsv_row, its row
% in that file, and so does every epoch made from it (EEG.epoch(k).eventtsv_row).
% EEG.ieeglab.opt.events holds exactly the rows still represented in EEG.event,
% linked by tsv_row rather than by position, so sorting, 'boundary' events and
% later deletions cannot misalign them.
%
% Cancelling either dialog returns the dataset unmodified with com empty.
%
% Cedric Cannard, iEEGLAB, 2025-2026

com = '';
EEG_in = EEG;
interactive = (nargin < 2 || isempty(opt));
if interactive, opt = struct(); end

if ndims(EEG.data) ~= 3 && (~isfield(EEG,'trials') || isempty(EEG.trials) || EEG.trials == 0)
    EEG.trials = 1;
end

% Labels with leading/trailing whitespace would match nothing downstream
lab = {EEG.chanlocs.labels};
trimmed = strtrim(lab);
if ~isequal(lab, trimmed)
    warning('ieeglab_load:labelWhitespace', ...
        'Trimming leading/trailing whitespace from %d channel label(s).', nnz(~strcmp(lab, trimmed)));
    [EEG.chanlocs.labels] = trimmed{:};
end

% ---------- which files ----------
auto = struct('elec_tsv', ieeglab_bids_sibling(EEG, 'electrodes'), ...
              'events_tsv', ieeglab_bids_sibling(EEG, 'events'), ...
              'channels_tsv', ieeglab_bids_sibling(EEG, 'channels'));
if interactive
    [g, wasCanceled] = ieeglab_gui_load1(EEG.filepath, auto);
    if wasCanceled || isempty(g), return; end
    opt.elec_tsv = g.elec_tsv; opt.events_tsv = g.events_tsv; opt.channels_tsv = g.channels_tsv;
else
    for f = {'elec_tsv','events_tsv','channels_tsv'}
        if ~isfield(opt, f{1}) || (ischar(opt.(f{1})) || isstring(opt.(f{1}))) && strcmpi(opt.(f{1}), 'auto')
            opt.(f{1}) = auto.(f{1});
        end
        opt.(f{1}) = char(opt.(f{1}));
    end
end
if ~isfield(opt,'drop_bad_events') || isempty(opt.drop_bad_events), opt.drop_bad_events = true; end

% ---------- electrodes ----------
elecs = [];
if ~isempty(opt.elec_tsv)
    elecs = readtable(opt.elec_tsv, 'FileType', 'text', 'Delimiter', '\t');
    vn = elecs.Properties.VariableNames;
    missingCols = setdiff({'x','y','z'}, lower(vn));
    if ~isempty(missingCols)
        warning('ieeglab_load:noCoordinateColumns', ...
            ['The electrodes file has no %s column(s), so 3D coordinates cannot be loaded.\n' ...
             'Columns found: %s\nContinuing without electrode coordinates.'], ...
             strjoin(upper(missingCols), '/'), strjoin(vn, ', '));
        elecs = [];
    else
        % Text columns (all 'n/a', comma decimals) are parsed, not assumed numeric
        xyz = local_numeric(elecs(:, ismember(lower(vn), {'x','y','z'})));
        if all(~isfinite(xyz(:)))
            warning('ieeglab_load:emptyCoordinates', ...
                'The x/y/z columns of %s hold no finite value (all n/a?). Continuing without them.', local_short(opt.elec_tsv));
            elecs = [];
        end
    end
end

% ---------- events ----------
opt.events_from_tsv = false;
if isfield(opt,'events'), opt = rmfield(opt, 'events'); end
if ~isempty(opt.events_tsv)
    events = readtable(opt.events_tsv, 'FileType', 'text', 'Delimiter', '\t');
    fprintf('%d events read from %s\n', height(events), local_short(opt.events_tsv));
    if ~ismember('onset', events.Properties.VariableNames)
        error('ieeglab_load:noOnset', ...
            ['The events file has no "onset" column, which BIDS requires and which is ' ...
             'needed to place events in the recording.\nColumns found: %s'], ...
             strjoin(events.Properties.VariableNames, ', '));
    end
    events.tsv_row = (1:height(events))';
    events.onset = local_numeric(events(:, 'onset'));
    % Only events that fall inside the recording can be placed. BIDS allows
    % negative onsets (before the first sample) and n/a onsets; drop them here,
    % explicitly, instead of letting eeg_checkset remove them from EEG.event
    % alone and leave the table one row out of step.
    lastOnset = (EEG.pnts - 1) / EEG.srate;
    bad = ~isfinite(events.onset) | events.onset < 0 | events.onset > lastOnset;
    if any(bad)
        warning('ieeglab_load:eventsOutsideData', ...
            'Dropping %d event(s) whose onset is n/a or outside the recording [0 %.2f] s (rows %s).', ...
            nnz(bad), lastOnset, mat2str(events.tsv_row(bad)'));
        events(bad,:) = [];
    end
    if opt.drop_bad_events && ismember('status', events.Properties.VariableNames)
        st = lower(strtrim(string(events.status)));
        isBadEv = st == "bad";
        if any(isBadEv)
            fprintf('Dropping %d/%d events marked status=bad in the events file.\n', nnz(isBadEv), numel(isBadEv));
            events(isBadEv,:) = [];
        end
    end
    opt.events = events;
    opt.events_from_tsv = true;
elseif isfield(EEG,'event') && ~isempty(EEG.event)
    fprintf('No events file; using the %d events already in the dataset.\n', numel(EEG.event));
    opt.events = struct2table(EEG.event, 'AsArray', true);   % lets the dialog list their fields
else
    fprintf('No events: continuous mode.\n');
end

% ---------- labels after EEGLAB's BIDS import ----------
% EEG-BIDS (bids_importchanlocs, up to at least 2026-05) takes channel labels
% from channels.tsv, then overwrites them row by row with electrodes.tsv. When
% electrodes.tsv lists fewer contacts or another order, which BIDS allows, data
% rows get the wrong labels (ds004696 sub-02: 27 of 30 wrong). channels.tsv rows
% follow the data order, so while every channel is still there, restore them.
if ~isempty(opt.channels_tsv) && isfield(EEG,'BIDS') && isstruct(EEG.BIDS)
    ct = readtable(opt.channels_tsv, 'FileType', 'text', 'Delimiter', '\t', 'TextType', 'char');
    if ismember('name', ct.Properties.VariableNames) && height(ct) == EEG.nbchan
        tsvNames = strtrim(cellstr(ct.name))';
        if ~isequal({EEG.chanlocs.labels}, tsvNames)
            nWrong = nnz(~strcmp({EEG.chanlocs.labels}, tsvNames));
            warning('ieeglab_load:bidsImportLabels', ...
                ['%d of %d channel labels differ from the row order of %s, which is the data order. ' ...
                 'EEGLAB''s BIDS import assigns electrodes.tsv rows to channels by position; restoring ' ...
                 'the labels from channels.tsv and matching coordinates by name.'], ...
                nWrong, EEG.nbchan, local_short(opt.channels_tsv));
            [EEG.chanlocs.labels] = tsvNames{:};
            for f = {'X','Y','Z','theta','radius','sph_theta','sph_phi','sph_radius'}
                if isfield(EEG.chanlocs, f{1}), [EEG.chanlocs.(f{1})] = deal([]); end
            end
        end
    end
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
    [EEG, badT] = ieeglab_bad_channels(EEG, struct('channels_tsv', opt.channels_tsv, ...
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
elseif isfield(opt,'events') && istable(opt.events) && ~isempty(opt.events) ...
        && (~isfield(opt,'event_field') || isempty(opt.event_field))
    [opt.event_field, opt.events] = local_default_event_field(opt.events, local_bids_event_map(EEG, opt.events_from_tsv));
end
if ~isfield(opt,'event_field'), opt.event_field = ''; end
opt.event_field = char(opt.event_field);

% ---------- place the events ----------
if opt.events_from_tsv && isfield(opt,'events') && ~isempty(opt.events)
    EEG = local_place_tsv_events(EEG, opt);
elseif ~opt.events_from_tsv && isfield(EEG,'event') && ~isempty(EEG.event) && ~isempty(opt.event_field)
    % EEGLAB's BIDS import moves the events.tsv column chosen as event type into
    % EEG.event.type and records the move in EEG.BIDS.eInfo: follow it.
    map = local_bids_event_map(EEG, false);
    if ~isfield(EEG.event, opt.event_field) && isfield(map, opt.event_field)
        fprintf('Event column "%s" was imported as EEG.event.%s by EEGLAB''s BIDS import; using it.\n', ...
            opt.event_field, map.(opt.event_field));
        opt.event_field = map.(opt.event_field);
    end
    EEG = local_retype_dataset_events(EEG, opt.event_field);
end
opt = local_align_table(EEG, opt);

% ---------- keep only events of interest ----------
if isfield(opt,'event_to_keep') && ~isempty(opt.event_to_keep) && isfield(EEG,'event') && ~isempty(EEG.event)
    want = string(opt.event_to_keep);
    types = string({EEG.event.type});
    isB = strcmpi(types, 'boundary');
    keep = ismember(types, want) | isB;
    fprintf('Keeping %d/%d events of interest (%s).\n', nnz(keep & ~isB), nnz(~isB), strjoin(cellstr(want), ', '));
    if ~any(keep & ~isB)
        error('ieeglab_load:noEventsKept', ...
            'None of the requested event types (%s) are present. Types available: %s', ...
            strjoin(cellstr(want), ', '), strjoin(unique(cellstr(types(~isB))), ', '));
    end
    EEG.event(~keep) = [];
    EEG = eeg_checkset(EEG, 'eventconsistency');
    opt = local_align_table(EEG, opt);
end

% ---------- keep only channels of interest ----------
removedNow = {};
if isfield(opt,'chan_list') && ~isempty(opt.chan_list)
    want = cellstr(string(opt.chan_list));
    allLabels = {EEG.chanlocs.labels};
    unknown = setdiff(want, allLabels);
    if ~isempty(unknown)
        warning('ieeglab_load:unknownChannels', 'Ignoring channel names not in the dataset: %s', strjoin(unknown, ', '));
    end
    removedNow = setdiff(allLabels, want, 'stable');
    if numel(removedNow) == numel(allLabels)
        error('ieeglab_load:noChannelsSelected', 'None of the requested channels are in the dataset.');
    elseif ~isempty(removedNow)
        fprintf('Removing %d/%d channels not selected: %s\n', numel(removedNow), numel(allLabels), strjoin(removedNow, ', '));
        EEG = pop_select(EEG, 'nochannel', removedNow);
    end
end
if ~isfield(EEG,'ieeglab') || isempty(EEG.ieeglab), EEG.ieeglab = struct(); end
prevRemoved = {};
if isfield(EEG.ieeglab,'removed_channels'), prevRemoved = EEG.ieeglab.removed_channels; end
allRemoved = unique([prevRemoved(:); removedNow(:)], 'stable');
if ~isempty(removedNow), EEG.ieeglab.removed_channels = allRemoved; end

% CCEP: a stimulation whose pair uses a contact no longer in the montage -
% removed now or in an earlier load - is dropped (exact token match).
if ~isempty(allRemoved) && isfield(EEG,'event') && ~isempty(EEG.event)
    drop = ieeglab_events_using(EEG.event, allRemoved);
    if any(drop)
        fprintf('Removing %d/%d stimulation events whose pair uses a removed contact.\n', nnz(drop), numel(drop));
        EEG.event(drop) = [];
        EEG = eeg_checkset(EEG, 'eventconsistency');
        opt = local_align_table(EEG, opt);
    end
end

% ---------- store ----------
% Options owned by the loader describe THIS call only; stale selections from a
% previous load are not carried over (the dataset would contradict them).
owned = {'elec_tsv','events_tsv','channels_tsv','event_field','event_to_keep','chan_list', ...
         'drop_bad_events','events','events_from_tsv','bad_labels','bad_reasons','elec_labels'};
base = struct();
if isfield(EEG.ieeglab,'opt') && isstruct(EEG.ieeglab.opt)
    base = rmfield(EEG.ieeglab.opt, intersect(fieldnames(EEG.ieeglab.opt), owned));
end
fns = fieldnames(opt);
for k = 1:numel(fns), base.(fns{k}) = opt.(fns{k}); end
EEG.ieeglab.opt = base;
EEG = eeg_checkset(EEG);

applied = struct('elec_tsv', opt.elec_tsv, 'events_tsv', opt.events_tsv, 'channels_tsv', opt.channels_tsv, ...
                 'event_field', opt.event_field, 'drop_bad_events', logical(opt.drop_bad_events));
if isfield(opt,'event_to_keep') && ~isempty(opt.event_to_keep), applied.event_to_keep = cellstr(string(opt.event_to_keep)); end
if isfield(opt,'chan_list') && ~isempty(opt.chan_list), applied.chan_list = cellstr(string(opt.chan_list)); end
com = sprintf('EEG = ieeglab_load(EEG, %s);', ieeglab_literal(applied));
end

% ======================= helpers =======================

function EEG = local_place_tsv_events(EEG, opt)
col = opt.event_field;
T = opt.events;
if ~ismember(col, T.Properties.VariableNames)
    error('ieeglab_load:badEventField', ...
        'Event field "%s" is not a column of the events file. Available: %s', ...
        col, strjoin(T.Properties.VariableNames, ', '));
end
ty = local_to_str(T.(col));
missing = ty == "" | lower(ty) == "n/a";
if any(missing)
    fprintf('Dropping %d/%d events with no value in "%s".\n', nnz(missing), numel(missing), col);
    T(missing,:) = []; ty(missing) = [];
end
hasDur = ismember('duration', T.Properties.VariableNames);
dur = [];
if hasDur, dur = local_numeric(T(:, 'duration')); end

% EEGLAB 'boundary' events mark discontinuities (EDF+D gaps, concatenated
% runs, removed segments). They are kept, so pop_epoch still rejects epochs
% that straddle one.
bnd = [];
if isfield(EEG,'event') && ~isempty(EEG.event) && isfield(EEG.event,'type')
    isB = strcmpi(cellfun(@(x) char(string(x)), {EEG.event.type}, 'UniformOutput', false), 'boundary');
    bnd = EEG.event(isB);
    if ~isempty(bnd) && isfield(bnd,'duration')
        d = [bnd.duration];
        if any(d(isfinite(d)) > 0)
            warning('ieeglab_load:boundaryWithRemovedData', ...
                ['The dataset has boundary events that mark removed data (duration > 0). ' ...
                 'events.tsv onsets are relative to the ORIGINAL recording and may no longer ' ...
                 'map to the right samples after data were removed.']);
        end
    end
end

ev = struct('type', cellstr(ty)', 'latency', num2cell(T.onset' * EEG.srate + 1), ...
            'tsv_row', num2cell(T.tsv_row'));
if hasDur
    for k = 1:numel(ev)
        if isfinite(dur(k)), ev(k).duration = dur(k) * EEG.srate; else, ev(k).duration = 0; end
    end
end
if ~isempty(bnd)
    ev = local_concat_events(ev, bnd);
end
EEG.event = ev;
EEG.urevent = [];
EEG = eeg_checkset(EEG, 'eventconsistency');
EEG = eeg_checkset(EEG, 'makeur');
fprintf('%d events placed using column "%s".\n', numel(ty), col);
end

function EEG = local_retype_dataset_events(EEG, field)
% Events already in the dataset (e.g. from EEGLAB's BIDS importer, which puts
% trial_type in EEG.event.type and every other column in its own field): use
% the chosen field as the event type, as the events.tsv path does.
if ~isfield(EEG.event, field)
    error('ieeglab_load:badEventField', 'Event field "%s" is not a field of EEG.event. Fields: %s', ...
        field, strjoin(fieldnames(EEG.event)', ', '));
end
n = numel(EEG.event); drop = false(1, n);
for k = 1:n
    t = char(string(EEG.event(k).type));
    if strcmpi(t, 'boundary'), continue; end
    v = EEG.event(k).(field);
    if isnumeric(v) || islogical(v), v = num2str(v); end
    v = strtrim(char(string(v)));
    if isempty(v) || strcmpi(v, 'n/a')
        drop(k) = true;
        continue
    end
    EEG.event(k).orig_type = t;
    EEG.event(k).type = v;
end
if any(drop)
    fprintf('Dropping %d/%d events with no value in field "%s".\n', nnz(drop), n, field);
    EEG.event(drop) = [];
end
EEG = eeg_checkset(EEG, 'eventconsistency');
fprintf('Event types taken from field "%s".\n', field);
end

function opt = local_align_table(EEG, opt)
% Keep opt.events to exactly the events.tsv rows still in EEG.event, in the
% same order, linked by tsv_row.
if ~opt.events_from_tsv || ~isfield(opt,'events') || ~istable(opt.events), return; end
if isempty(EEG.event) || ~isfield(EEG.event, 'tsv_row')
    opt.events = opt.events([],:);
    return
end
r = [EEG.event.tsv_row];
r = r(isfinite(r));
[tf, loc] = ismember(r, opt.events.tsv_row);
opt.events = opt.events(loc(tf), :);
end

function ev = local_concat_events(a, b)
% Concatenate two event struct arrays whose fields differ.
f = union(fieldnames(a), fieldnames(b), 'stable');
for i = 1:numel(f)
    if ~isfield(a, f{i}), [a.(f{i})] = deal([]); end
    if ~isfield(b, f{i}), [b.(f{i})] = deal([]); end
end
a = orderfields(a, f); b = orderfields(b, f);
for k = 1:numel(b), if isempty(b(k).tsv_row), b(k).tsv_row = NaN; end, end
ev = [a(:); b(:)]';
end

function [f, T] = local_default_event_field(T, map)
% map: BIDS column -> EEG.event field, for columns EEGLAB's BIDS import renamed.
vn = T.Properties.VariableNames;
for c = {'electrical_stimulation_site','trial_type','type','value'}
    if isfield(map, c{1}) && ismember(map.(c{1}), vn), f = map.(c{1}); return; end
    if ismember(c{1}, vn), f = c{1}; return; end
end
rest = setdiff(vn, {'onset','duration','sample','sample_start','tsv_row','latency','urevent'}, 'stable');
if isempty(rest)
    % BIDS only requires onset and duration; give every event one generic type
    warning('ieeglab_load:noEventColumn', ...
        'The events file has no column naming the events (only %s). Using the type ''event'' for all.', strjoin(vn, ', '));
    T.event_type = repmat({'event'}, height(T), 1);
    f = 'event_type';
else
    f = rest{1};
end
end

function map = local_bids_event_map(EEG, fromTsv)
% Columns that EEGLAB's BIDS import (EEG-BIDS) stored under another EEG.event
% field, from the {column, field} pairs it keeps in EEG.BIDS.eInfo. Only
% relevant for events already in the dataset, not for events read from a file.
map = struct();
if fromTsv || ~isfield(EEG,'BIDS') || ~isstruct(EEG.BIDS) || ~isfield(EEG.BIDS,'eInfo'), return; end
e = EEG.BIDS.eInfo;
if ~iscell(e) || size(e, 2) < 2, return; end
for k = 1:size(e, 1)
    if ischar(e{k,1}) && ischar(e{k,2}) && isvarname(e{k,1}) && ~strcmp(e{k,1}, e{k,2})
        map.(e{k,1}) = e{k,2};
    end
end
end

function x = local_numeric(T)
% Numeric matrix from table columns that may be text ('n/a', comma decimals).
x = nan(height(T), width(T));
for j = 1:width(T)
    v = T{:, j};
    if isnumeric(v) || islogical(v)
        x(:, j) = double(v);
    else
        x(:, j) = str2double(strrep(string(v), ',', '.'));
    end
end
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

function v = local_or_table(t)
if isempty(t), v = ''; else, v = t; end
end

function s = local_short(p)
[~, n, e] = fileparts(char(p)); s = [n e];
end
