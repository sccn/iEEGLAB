function [EEG, wasCanceled] = ieeglab_gui_preprocess(EEG)
% iEEGLAB: Preprocessing Options GUI
% - CAR kept
% - Single global filter type (noncausal vs minimum-phase) for ALL filters
% - If no events exist, Epoch/CAR/Baseline controls are disabled and values forced OFF

wasCanceled = false;

% ------------------------------------------------------------------------
%  Defaults
% -------------------------------------------------------------------------
ds_should_enable = EEG.srate > 512;
ds_default_rate  = 512;

choices = struct();
choices.event_filters     = struct();   % filled by sub-GUI (or empty)
choices.remove_rare_cond  = true;
choices.min_trials        = 5;
choices.remove_no_coords  = true;

choices.apply_ds          = ds_should_enable;

choices.apply_highpass    = true;
choices.highpass          = 0.1;        % Hz

choices.apply_notch       = true;
choices.notch             = 60;         % Hz (scalar or vector, e.g., [60 120])

choices.apply_lowpass     = false;
choices.lowpass           = 40;         % Hz (empty -> none)

% Global filter type for ALL filters
choices.filter_type       = 1;          % 1 = Noncausal zero-phase (default), 2 = Minimum-phase

% Epoching (enabled by default IF events exist)
choices.apply_epoch       = true;
choices.epoch_window      = [-500 900]; % ms

% CAR (enabled by default IF events exist)
choices.apply_acar        = true;
choices.acar_fraction     = 0.20;
choices.acar_timewin      = [15 500];   % ms

% Baseline (enabled by default IF events exist)
choices.apply_baseline    = true;
choices.baseline_method   = 1;          % 1=median, 2=mean, 3=trimmed mean, 4=1/f
choices.baseline_period   = [-500 -50]; % ms
choices.baseline_mode     = 1;          % 1=subtract, 2=divide

% Downsample default
if ds_should_enable
    choices.downsample    = ds_default_rate;
else
    choices.downsample    = max(1, EEG.srate);
end

% -------------------------------------------------------------------------
%  Events availability (prefer stored table; fallback to EEG.event)
% -------------------------------------------------------------------------
evT = table(); showCols = {};
if isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'opt') && ...
   isfield(EEG.ieeglab.opt,'events') && istable(EEG.ieeglab.opt.events) && ~isempty(EEG.ieeglab.opt.events)
    evT = EEG.ieeglab.opt.events;
elseif isfield(EEG,'event') && ~isempty(EEG.event)
    try
        if istable(EEG.event), evT = EEG.event;
        elseif isstruct(EEG.event), evT = struct2table(EEG.event,'AsArray',true);
        end
    catch
        evT = table();
    end
end
haveEv = istable(evT) && ~isempty(evT) && width(evT)>0;

if haveEv
    ignoreCols = {'onset','duration','electrodes_involved_onset','electrodes_involved_offset','sample_start'};
    allCols    = evT.Properties.VariableNames;
    showCols   = setdiff(allCols, ignoreCols, 'stable');
    mustInclude = intersect({'stim_type','stim_num'}, allCols, 'stable');
    showCols   = unique([showCols(:); mustInclude(:)], 'stable');
end

% If NO events: force these OFF and keep their controls disabled
if ~haveEv
    choices.apply_epoch    = false;
    choices.apply_acar     = false;
    choices.apply_baseline = false;
end

% -------------------------------------------------------------------------
%  Callbacks and on/off states
% -------------------------------------------------------------------------
if haveEv && ~isempty(showCols)
    cb_ev = @(h,~) local_open_event_selector(h, evT, showCols);
    btn_en = 'on';
else
    cb_ev = @(~,~) warndlg('No events available to select.','iEEGLAB');
    btn_en = 'off';
end

% generic togglers
cb_ds  = local_cb_toggle({'lbl_ds_rate','ds_rate'});
cb_hp  = local_cb_toggle({'lbl_highpass','highpass'});
cb_no  = local_cb_toggle({'lbl_notch','notch'});
cb_lp  = local_cb_toggle({'lbl_lowpass','lowpass'});

% these depend on events
if haveEv
    cb_seg = local_cb_toggle({'lbl_epoch','epoch_window'});
    cb_acar= local_cb_toggle({'acar_fraction','acar_timewin','lbl_acar_fraction','lbl_acar_timewin'});
    cb_bl  = local_cb_toggle({'lbl_bl_method','bl_method','lbl_bl_period','bl_period','lbl_bl_mode','bl_mode'});
else
    cb_seg = @(h,~) set(h,'value',0);
    cb_acar= @(h,~) set(h,'value',0);
    cb_bl  = @(h,~) set(h,'value',0);
end

% initial enable states
onoff_ds  = iff(choices.apply_ds,'on','off');
onoff_hp  = iff(choices.apply_highpass,'on','off');
onoff_no  = iff(choices.apply_notch,'on','off');
onoff_lp  = iff(choices.apply_lowpass,'on','off');

onoff_seg = iff(haveEv && choices.apply_epoch,'on','off');
onoff_car = iff(haveEv && choices.apply_acar,'on','off');
onoff_bl  = iff(haveEv && choices.apply_baseline,'on','off');

% -------------------------------------------------------------------------
%  Geometry (rows match UI controls order)
% -------------------------------------------------------------------------
uigeom = {
    1
    [0.70 0.30]
    [0.70 0.30]
    [0.06 0.64 0.30]
    0.01
    1
    [0.70 0.30]
    1

    % Downsample
    [0.70 0.30]
    [0.06 0.64 0.30]

    % High-pass
    [0.70 0.30]
    [0.06 0.64 0.30]

    % Notch
    [0.70 0.30]
    [0.06 0.64 0.30]

    % Low-pass
    [0.70 0.30]
    [0.06 0.64 0.30]

    % Global Filter Type
    [0.70 0.30]

    % Segmentation
    [0.70 0.30]
    [0.06 0.64 0.30]

    % CAR
    [0.70 0.30]
    [0.06 0.64 0.30]
    [0.06 0.64 0.30]

    % Baseline
    [0.70 0.30]
    [0.06 0.64 0.30]
    [0.06 0.64 0.30]
    [0.06 0.64 0.30]
};

% -------------------------------------------------------------------------
%  Build UI
% -------------------------------------------------------------------------
uilist = {};  add = @(c) assignin('caller','uilist',[uilist,{c}]); %#ok<NASGU>
uilist = {};

% helper append
    function append(c), uilist{end+1} = c; end

% Events
append({'style' 'text' 'string' 'Events' 'fontweight' 'bold' 'horizontalalignment' 'left'});
append({'style' 'text' 'string' 'Select events of interest:' 'horizontalalignment' 'left'});
append({'style' 'pushbutton' 'string' 'Open selector…' 'callback' cb_ev 'enable' btn_en});

append({'style' 'text' 'string' 'Remove rare conditions:' 'horizontalalignment' 'left'});
append({'style' 'checkbox' 'tag' 'remove_rare_cond' 'value' choices.remove_rare_cond 'string' 'Enable'});

append({'style' 'text' 'string' ''});
append({'style' 'text' 'string' 'Min trials:' 'horizontalalignment' 'left'});
append({'style' 'edit' 'tag' 'min_trials' 'string' num2str(choices.min_trials)});

% hidden event filter store
append({'style' 'edit' 'tag' 'evsel_json' 'string' '' 'visible' 'off'});

% Electrodes
append({'style' 'text' 'string' 'Electrodes' 'fontweight' 'bold' 'horizontalalignment' 'left'});
append({'style' 'text' 'string' 'Remove electrodes with no XYZ coordinates?' 'horizontalalignment' 'left'});
append({'style' 'checkbox' 'tag' 'rm_no_coords' 'value' choices.remove_no_coords 'string' ''});

% Signal Processing header
append({'style' 'text' 'string' 'Signal Processing' 'fontweight' 'bold' 'horizontalalignment' 'left'});

% Downsample
append({'style' 'text' 'string' 'Downsample:' 'horizontalalignment' 'left'});
append({'style' 'checkbox' 'tag' 'apply_ds' 'value' choices.apply_ds 'string' 'Enable' 'callback' cb_ds});
append({'style' 'text' 'string' ''});
append({'style' 'text' 'tag' 'lbl_ds_rate' 'string' 'Target rate (Hz):' 'horizontalalignment' 'left' 'enable' onoff_ds});
append({'style' 'edit' 'tag' 'ds_rate' 'string' num2str(choices.downsample) 'enable' onoff_ds});

% High-pass
append({'style' 'text' 'string' 'High-pass filter:' 'horizontalalignment' 'left'});
append({'style' 'checkbox' 'tag' 'apply_highpass' 'value' choices.apply_highpass 'string' 'Enable' 'callback' cb_hp});
append({'style' 'text' 'string' ''});
append({'style' 'text' 'tag' 'lbl_highpass' 'string' 'Cutoff (Hz):' 'horizontalalignment' 'left' 'enable' onoff_hp});
append({'style' 'edit' 'tag' 'highpass' 'string' num2str(choices.highpass) 'enable' onoff_hp});

% Notch
append({'style' 'text' 'string' 'Notch filter:' 'horizontalalignment' 'left'});
append({'style' 'checkbox' 'tag' 'apply_notch' 'value' choices.apply_notch 'string' 'Enable' 'callback' cb_no});
append({'style' 'text' 'string' ''});
append({'style' 'text' 'tag' 'lbl_notch' 'string' 'Center(s) (Hz):' 'horizontalalignment' 'left' 'enable' onoff_no});
append({'style' 'edit' 'tag' 'notch' 'string' num2str(choices.notch) 'enable' onoff_no});

% Low-pass
append({'style' 'text' 'string' 'Low-pass filter:' 'horizontalalignment' 'left'});
append({'style' 'checkbox' 'tag' 'apply_lowpass' 'value' choices.apply_lowpass 'string' 'Enable' 'callback' cb_lp});
append({'style' 'text' 'string' ''});
append({'style' 'text' 'tag' 'lbl_lowpass' 'string' 'Cutoff (Hz):' 'horizontalalignment' 'left' 'enable' onoff_lp});
append({'style' 'edit' 'tag' 'lowpass' 'string' iff(isempty(choices.lowpass),'',num2str(choices.lowpass)) 'enable' onoff_lp});

% Global filter type (applies to all filters)
append({'style' 'text'  'string' 'Filter type (all filters):' 'horizontalalignment' 'left'});
append({'style' 'popupmenu' 'tag' 'filter_type' ...
        'string' {'Noncausal zero-phase (default)' 'Minimum-phase'} ...
        'value' choices.filter_type});

% Segmentation
append({'style' 'text' 'string' 'Segmentation (epoching):' 'horizontalalignment' 'left'});
append({'style' 'checkbox' 'tag' 'apply_epoch' 'value' choices.apply_epoch ...
        'string' 'Enable' 'callback' cb_seg 'enable' iff(haveEv,'on','off')});
append({'style' 'text' 'string' ''});
append({'style' 'text' 'tag' 'lbl_epoch' 'string' 'Epoch window [ms] (start end):' 'horizontalalignment' 'left' 'enable' onoff_seg});
append({'style' 'edit' 'tag' 'epoch_window' 'string' sprintf('%d %d',choices.epoch_window) 'enable' onoff_seg});

% CAR
append({'style' 'text' 'string' 'Adjusted Average Reference:' 'horizontalalignment' 'left'});
append({'style' 'checkbox' 'tag' 'apply_acar' 'value' choices.apply_acar ...
        'string' 'Enable' 'callback' cb_acar 'enable' iff(haveEv,'on','off')});
append({'style' 'text' 'string' ''});
append({'style' 'text' 'tag' 'lbl_acar_fraction' 'string' 'Fraction of channels (0–1):' 'horizontalalignment' 'left' 'enable' onoff_car});
append({'style' 'edit' 'tag' 'acar_fraction' 'string' num2str(choices.acar_fraction) 'enable' onoff_car});
append({'style' 'text' 'string' ''});
append({'style' 'text' 'tag' 'lbl_acar_timewin' 'string' 'CAR window [ms] (start end):' 'horizontalalignment' 'left' 'enable' onoff_car});
append({'style' 'edit' 'tag' 'acar_timewin' 'string' sprintf('%d %d',choices.acar_timewin) 'enable' onoff_car});

% Baseline
append({'style' 'text' 'string' 'Remove baseline:' 'horizontalalignment' 'left'});
append({'style' 'checkbox' 'tag' 'apply_baseline' 'value' choices.apply_baseline ...
        'string' 'Enable' 'callback' cb_bl 'enable' iff(haveEv,'on','off')});
append({'style' 'text' 'string' ''});
append({'style' 'text' 'tag' 'lbl_bl_method' 'string' 'Method:' 'horizontalalignment' 'left' 'enable' onoff_bl});
append({'style' 'popupmenu' 'tag' 'bl_method' 'string' {'Median (default)' 'Mean' 'Trimmed mean' '1/f'} 'value' choices.baseline_method 'enable' onoff_bl});
append({'style' 'text' 'string' ''});
append({'style' 'text' 'tag' 'lbl_bl_period' 'string' 'Baseline period [ms] (start end):' 'horizontalalignment' 'left' 'enable' onoff_bl});
append({'style' 'edit' 'tag' 'bl_period' 'string' sprintf('%d %d',choices.baseline_period) 'enable' onoff_bl});
append({'style' 'text' 'string' ''});
append({'style' 'text' 'tag' 'lbl_bl_mode' 'string' 'Mode:' 'horizontalalignment' 'left' 'enable' onoff_bl});
append({'style' 'popupmenu' 'tag' 'bl_mode' 'string' {'Subtract (default)' 'Divide'} 'value' choices.baseline_mode 'enable' onoff_bl});

% -------------------------------------------------------------------------
%  Launch
% -------------------------------------------------------------------------
[res, ~, ~, out] = inputgui('geometry', uigeom, 'uilist', uilist, ...
                            'title', 'iEEGLAB: Preprocessing Options', ...
                            'minwidth', 440);
if isempty(res) || isempty(out)
    fprintf('iEEGLAB preprocess dialog canceled\n');
    wasCanceled = true;
    return
end

% -------------------------------------------------------------------------
%  Parse outputs (with guards) + finalize consistency with event availability
% -------------------------------------------------------------------------
% event filters
if isfield(out,'evsel_json') && ~isempty(out.evsel_json)
    try, choices.event_filters = jsondecode(out.evsel_json); catch, choices.event_filters = struct(); end
end

% numbers & toggles
choices.remove_rare_cond = isfield(out,'remove_rare_cond') && logical(out.remove_rare_cond);
choices.remove_no_coords = isfield(out,'rm_no_coords') && logical(out.rm_no_coords);

if isfield(out,'min_trials'),  mt = str2double(out.min_trials);  if isfinite(mt) && mt>=0, choices.min_trials = mt; end, end
choices.apply_ds = isfield(out,'apply_ds') && logical(out.apply_ds);
if isfield(out,'ds_rate'),     dr = str2double(out.ds_rate);     if isfinite(dr) && dr>0, choices.downsample = dr; end, end

choices.apply_highpass = isfield(out,'apply_highpass') && logical(out.apply_highpass);
if isfield(out,'highpass'),    hp = str2double(out.highpass);    if isfinite(hp) && hp>0, choices.highpass = hp; end, end

choices.apply_notch = isfield(out,'apply_notch') && logical(out.apply_notch);
if isfield(out,'notch') && ~isempty(out.notch)
    nv = str2num(out.notch); %#ok<ST2NM>
    if isempty(nv), nv = str2double(out.notch); end
    if all(isfinite(nv)), choices.notch = nv; end
end

choices.apply_lowpass = isfield(out,'apply_lowpass') && logical(out.apply_lowpass);
if isfield(out,'lowpass') && ~isempty(out.lowpass)
    lp = str2double(out.lowpass); if isfinite(lp) && lp>0, choices.lowpass = lp; else, choices.lowpass = []; end
else
    if ~choices.apply_lowpass, choices.lowpass = []; end
end

% global filter type
if isfield(out,'filter_type') && ~isempty(out.filter_type)
    choices.filter_type = out.filter_type;   % 1 or 2
end
idx = min(max(round(choices.filter_type),1),2);
ft_labels = {'noncausal-zero-phase','minimum-phase'};
choices.filter_type_label = ft_labels{idx};

% epoch / CAR / baseline (respect event availability)
choices.apply_epoch = isfield(out,'apply_epoch') && logical(out.apply_epoch);
if isfield(out,'epoch_window') && ~isempty(out.epoch_window)
    tw = sscanf(out.epoch_window,'%f'); if numel(tw)>=2, choices.epoch_window = tw(1:2).'; end
end

choices.apply_acar = isfield(out,'apply_acar') && logical(out.apply_acar);
if isfield(out,'acar_fraction') && ~isempty(out.acar_fraction)
    choices.acar_fraction = max(0, min(1, str2double(out.acar_fraction)));
end
if isfield(out,'acar_timewin') && ~isempty(out.acar_timewin)
    tw = sscanf(out.acar_timewin,'%f'); if numel(tw)>=2, choices.acar_timewin = tw(1:2).'; end
end

choices.apply_baseline = isfield(out,'apply_baseline') && logical(out.apply_baseline);
bl_labels = {'median','mean','trimmed mean','1/f'};
if isfield(out,'bl_method') && ~isempty(out.bl_method), bl_idx = out.bl_method; else, bl_idx = choices.baseline_method; end
choices.baseline_method = bl_labels{min(max(bl_idx,1),4)};
if isfield(out,'bl_period') && ~isempty(out.bl_period)
    bp = sscanf(out.bl_period,'%f'); if numel(bp)>=2, choices.baseline_period = bp(1:2).'; end
end
bm_labels = {'subtract','divide'};
if isfield(out,'bl_mode') && ~isempty(out.bl_mode), bm_idx = out.bl_mode; else, bm_idx = choices.baseline_mode; end
choices.baseline_mode = bm_labels{min(max(bm_idx,1),2)};

% FINAL consistency: if no events, force these OFF (ignores any GUI value)
if ~haveEv
    choices.apply_epoch    = false;
    choices.apply_acar     = false;
    choices.apply_baseline = false;
end

%% ------------------------------------------------------------------------
%  Save back into EEG
% -------------------------------------------------------------------------
if ~isfield(EEG,'ieeglab'), EEG.ieeglab = struct(); end
if ~isfield(EEG.ieeglab,'opt') || ~isstruct(EEG.ieeglab.opt), EEG.ieeglab.opt = struct(); end
EEG.ieeglab.opt = merge_structs(EEG.ieeglab.opt, choices);

end

%% ========================================================================
%  Sub-GUI: event selector
% ========================================================================
function local_open_event_selector(srcBtn, evT, showCols)
lists = cell(1,numel(showCols)); tags = cell(1,numel(showCols));
for k = 1:numel(showCols)
    col = showCols{k}; vals = evT.(col);
    if iscell(vals), vals = string(vals); end
    if iscategorical(vals), vals = string(vals); end
    if isnumeric(vals), vals = string(vals); end
    u = unique(vals(~ismissing(vals)));
    lists{k} = [{'All'}; cellstr(u(:))];
    tags{k}  = ['fld_' regexprep(col,'[^\w]','_')];
end
nRows = 1+numel(showCols); uigeom = cell(nRows,1); gvert = ones(1,nRows);
uilist = cell(0,1);
uilist{end+1} = {'style' 'text' 'string' 'Select events to KEEP (leave at "All")' 'fontweight' 'bold'};
uigeom{1} = [1]; gvert(1)=1; ii=2;
for k=1:numel(showCols)
    uigeom{ii}=[1 1]; gvert(ii)=1.4;
    uilist{end+1}={'style' 'text' 'string' showCols{k} 'horizontalalignment' 'left'};
    uilist{end+1}={'style' 'listbox' 'tag' tags{k} 'string' lists{k} 'value' 1 'min' 0 'max' 2};
    ii=ii+1;
end
[res,~,~,out]=inputgui('geometry',uigeom,'geomvert',gvert,'uilist',uilist,'title','Select events of interest','minwidth',440);
if isempty(res)||isempty(out), return; end
filters=struct();
for k=1:numel(showCols)
    tg=tags{k}; selIdx=1;
    if isfield(out,tg)&&~isempty(out.(tg)), selIdx=out.(tg); end
    if any(selIdx==1), filters.(showCols{k})={}; else, filters.(showCols{k})=lists{k}(selIdx); end
end
hFig=ancestor(srcBtn,'figure'); hStore=findobj(hFig,'tag','evsel_json');
if ~isempty(hStore)&&ishghandle(hStore)
    try, set(hStore,'string',jsonencode(filters));
    catch
        keys=fieldnames(filters); fr=strings(0);
        for i=1:numel(keys)
            v=string(filters.(keys{i})); if isempty(v), v="<ALL>"; end
            fr(end+1)=keys{i}+"="+strjoin(v,"|"); %#ok<AGROW>
        end
        set(hStore,'string',char(strjoin(fr,";")));
    end
end
end

%%  Helpers

function out = iff(cond,a,b), if cond, out=a; else, out=b; end, end
function s = merge_structs(a,b), s=a; f=fieldnames(b); for i=1:numel(f), s.(f{i})=b.(f{i}); end, end
function cb = local_cb_toggle(tag_list), cb = @(h,~) local_set_enable(ancestor(h,'figure'), tag_list, get(h,'value')~=0); end
function local_set_enable(hFig, tag_list, is_on)
    onoff = iff(is_on,'on','off');
    for k = 1:numel(tag_list)
        h = findobj(hFig,'tag',tag_list{k});
        if ~isempty(h) && ishghandle(h), set(h,'enable',onoff); end
    end
end
