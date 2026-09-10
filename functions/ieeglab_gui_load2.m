function [opt, wasCanceled] = ieeglab_gui_load2(opt)
% iEEGLAB GUI (channels + optional event filters)
%
% Required (preferred):
%   opt.elec_labels  - cellstr of channel labels
%
% If opt.elec_labels is missing/empty AND there is no opt.elec_tsv,
% this function will fall back to {EEG.chanlocs.labels}.
%
% Optional (preselects only):
%   opt.chan_idx     - previously selected channel indices (preselect)
%   opt.events       - table; if present enables event-based filters UI
%   opt.event_field  - (optional) preselect events column name
%   opt.event_to_keep - (optional) preselect values (cell of char) for that column
%
% Outputs set by GUI:
%   opt.chan_idx, opt.chan_list
%   opt.event_field          - char (present only when events exist)
%   opt.event_to_keep         - {} when "All"/none selected; otherwise list of values (cellstr)
%   wasCanceled              - logical, true if user canceled

wasCanceled = false;

% ---------- Resolve channel labels ----------
% Priority:
%   1) opt.elec_labels (if provided and non-empty)
%   2) If NO opt.elec_tsv, fall back to EEG.chanlocs.labels from workspace
if ~isfield(opt,'elec_labels') || isempty(opt.elec_labels)
    if ~isfield(opt,'elec_tsv') || isempty(opt.elec_tsv)
        EEG = [];
        try
            if evalin('base','exist(''EEG'',''var'')'); EEG = evalin('base','EEG'); end
        catch, end
        if isempty(EEG)
            try
                if evalin('caller','exist(''EEG'',''var'')'); EEG = evalin('caller','EEG'); end
            catch, end
        end
        if isempty(EEG) || ~isstruct(EEG) || ~isfield(EEG,'chanlocs') || isempty(EEG.chanlocs)
            error(['ieeglab_gui_load2: No opt.elec_labels and no opt.elec_tsv. ', ...
                   'EEG with non-empty .chanlocs not found in workspace.']);
        end
        lbls = {EEG.chanlocs.labels};
        if isempty(lbls) || any(cellfun(@isempty, lbls))
            error('ieeglab_gui_load2: EEG.chanlocs.labels is missing or contains empty labels.');
        end
        opt.elec_labels = lbls;
    else
        error(['ieeglab_gui_load2: opt.elec_labels is empty. Provide labels (from elec_tsv) ', ...
               'or omit opt.elec_tsv to fall back to EEG.chanlocs.labels.']);
    end
end

labels = opt.elec_labels(:);
nch    = numel(labels);

% ---------- Channels preselect ----------
% Clinician-marked bad channels (channels.tsv status) start UNselected and are
% labelled with the reason, so the default choice is the clinical one and
% overriding it is one click.
isBadCh = false(nch, 1);
badNote = repmat({''}, nch, 1);
if isfield(opt,'bad_labels') && ~isempty(opt.bad_labels)
    [isBadCh, loc] = ismember(upper(labels), upper(opt.bad_labels(:)));
    if isfield(opt,'bad_reasons') && numel(opt.bad_reasons) == numel(opt.bad_labels)
        badNote(isBadCh) = opt.bad_reasons(loc(isBadCh));
    end
end
if ~isfield(opt,'chan_idx') || isempty(opt.chan_idx)
    if any(isBadCh), preCh = find(~isBadCh)' + 1; else, preCh = 1; end
else
    sel  = opt.chan_idx(:)';
    sel  = sel(sel>=1 & sel<=nch);
    preCh = iff(isempty(sel) || numel(sel)==nch, 1, sel+1);
end
labels_shown = labels;
for k = find(isBadCh)'
    if isempty(badNote{k}), labels_shown{k} = [labels{k} '   [marked bad]'];
    else, labels_shown{k} = [labels{k} '   [bad: ' badNote{k} ']']; end
end
labels_disp = [{'All channels'}; labels_shown];

% ---------- Events presence & lists ----------
haveEvents = isfield(opt,'events') && istable(opt.events) && ~isempty(opt.events);

fields_disp = {};
all_value_lists = {};
if haveEvents
    vars = opt.events.Properties.VariableNames;
    exclude = {'onset','latency','duration','channel','bvmknum','bvtime','urevent','code','visible'};
    cand = setdiff(vars, intersect(vars, exclude), 'stable');
    if isempty(cand), cand = vars; end

    for k = 1:numel(cand)
        col  = cand{k};
        vals = getUniques(opt.events, col);
        fields_disp{end+1}     = col; %#ok<AGROW>
        all_value_lists{end+1} = [{'All'}; cellstr(string(vals(:)))]; %#ok<AGROW>
    end
    if isempty(fields_disp)
        fields_disp     = {'<no event fields>'};
        all_value_lists = {{'All'}};
    end

    if isfield(opt,'event_field') && ~isempty(opt.event_field)
        fidx = find(strcmp(fields_disp, opt.event_field), 1);
        preFieldIdx = iff(~isempty(fidx), fidx, 1);
    else
        fidx = find(strcmpi(fields_disp, 'type'), 1);
        if ~isempty(fidx)
            preFieldIdx = fidx;
        else
            normNames = lower(strrep(fields_disp, '_', ' '));
            fidx = find(contains(normNames, 'stimulation site'), 1);
            preFieldIdx = iff(~isempty(fidx), fidx, 1);
        end
    end
else
    fields_disp     = {'<no events table>'};
    all_value_lists = {{'All'}};
    preFieldIdx     = 1;
end

list_values = all_value_lists{preFieldIdx};

% Preselect values (from opt.event_to_keep)
preValIdx = 1;
if haveEvents && isfield(opt,'event_to_keep') && ~isempty(opt.event_to_keep)
    want = cellstr(opt.event_to_keep(:));
    idxs = [];
    for w = 1:numel(want)
        i = find(strcmp(list_values, want{w}), 1);
        if ~isempty(i), idxs(end+1) = i; end %#ok<AGROW>
    end
    if ~isempty(idxs), preValIdx = idxs; end
end

% ---------- Geometry ----------
uigeom = {
    [1]
    [1]
    [1]
    [1 1]
    [1 1]
};
geomvert = [1 6 1 1 2];
ev_enable = iff(haveEvents,'on','off');

% ---------- UI ----------
uilist = {
    {'style' 'text'    'string' 'Select channel(s):' 'fontweight' 'bold'}
    {'style' 'listbox' 'tag' 'sel_list' 'string' labels_disp ...
     'value' preCh 'min' 0 'max' 2 'ListboxTop' 1}

    {'style' 'text' 'string' 'Event filters' 'fontweight' 'bold' 'enable' ev_enable}

    {'style' 'text' 'tag' 'ev_field_label' 'string' 'Column holding the event/condition name:' 'enable' ev_enable}
    {'style' 'popupmenu' 'tag' 'event_field' 'string' fields_disp ...
     'value' preFieldIdx ...
     'userdata' all_value_lists 'enable' ev_enable ...
     'callback', @cb_field}

    {'style' 'text' 'tag' 'ev_values_label' 'string' 'Values to load:' 'enable' ev_enable}
    {'style' 'listbox' 'tag' 'event_val_list' 'string' list_values ...
     'value' preValIdx 'min' 0 'max' 2 'ListboxTop' 1 'enable' ev_enable}
};

[res,~,~,outstruct] = inputgui('geometry', uigeom, ...
                               'geomvert', geomvert, ...
                               'uilist', uilist, ...
                               'title','iEEGLAB: load options', ...
                               'minwidth', 600);
if isempty(res) || isempty(outstruct)
    fprintf('iEEGLAB load dialog canceled\n');
    wasCanceled = true;
    opt = [];
    return;
end

% ---------- Channels mapping ----------
idxCh = outstruct.sel_list(:);
if isempty(idxCh)
    error('ieeglab_gui_load2:noChannels', 'No channel selected. Select at least one channel, or "All channels".');
end
% 'All channels' is entry 1. Selecting it together with other entries used to
% fall into the else branch and index the mask at 0; it now means all.
if any(idxCh == 1)
    opt.chan_idx  = 1:nch;
    opt.chan_list = labels;
else
    mask = false(1, nch); % convert to a logical mask
    mask(idxCh-1) = true; % -1 to ignore 'All Channels' field
    opt.chan_idx  = mask;
    opt.chan_list = labels(mask);
end

% ---------- Event outputs ----------
if haveEvents
    if isfield(outstruct,'event_field') && ~isempty(outstruct.event_field)
        ef_idx = outstruct.event_field;
    else
        ef_idx = preFieldIdx;
    end
    ef_idx = max(1, min(ef_idx, numel(fields_disp)));
    opt.event_field = fields_disp{ef_idx};

    vals_for_field = all_value_lists{ef_idx};
    if isfield(outstruct,'event_val_list') && ~isempty(outstruct.event_val_list)
        opt.event_to_keep = mapOutEmptyOnAll(outstruct.event_val_list, vals_for_field);
    else
        opt.event_to_keep = {};
    end
else
    opt = rmfield_if_exists(opt, {'event_field','event_to_keep'});
end

end % === main ===

% -------------------- callbacks & helpers --------------------
function cb_field(src, ~)
    idx = get(src,'value');
    vallists = get(src,'userdata');
    lst = vallists{idx};
    hlist = findobj(gcbf,'tag','event_val_list');
    set(hlist,'string',lst,'value',1);
end

function out = mapOutEmptyOnAll(idx_list, list_disp)
    if isempty(idx_list) || any(idx_list==1)
        out = {};
    else
        out = list_disp(idx_list);
        out = out(:)'; % row cell
    end
end

function vals = getUniques(ev, colname)
    if ~ismember(colname, ev.Properties.VariableNames)
        vals = string.empty(0,1);
        return
    end
    x = ev.(colname);

    if isstring(x)
        vals = unique(x(:));
        vals = vals(~ismissing(vals));
        return
    elseif isnumeric(x) || islogical(x) || isdatetime(x) || isduration(x)
        vals = string(unique(x(:)));
        return
    elseif ischar(x)
        vals = string({x});
        return
    end

    if iscell(x)
        y = strings(size(x));
        for i = 1:numel(x)
            y(i) = local_scalar_to_string(x{i});
        end
        vals = unique(y(:));
        return
    end

    vals = string.empty(0,1);
end

function s = local_scalar_to_string(v)
    if isempty(v), s = ""; return; end
    if isstring(v)
        s = iff(isscalar(v), v, join(v(:), ","));
    elseif ischar(v)
        s = string(v);
    elseif isnumeric(v) || islogical(v)
        s = string(iff(isscalar(v), v, sprintf("[len %d]", numel(v))));
    elseif isdatetime(v) || isduration(v)
        s = string(v);
    else
        s = string(class(v));
    end
end

function y = iff(cond, a, b)
    if cond, y = a; else, y = b; end
end

function s = rmfield_if_exists(s, names)
    for i = 1:numel(names)
        if isfield(s, names{i}), s = rmfield(s, names{i}); end
    end
end
