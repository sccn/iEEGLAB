function [EEG, com] = ieeglab_preprocess(EEG, opt)
% ieeglab_preprocess() - Preprocess iEEG data: event selection, resampling,
%                        filtering, epoching, re-referencing, baseline removal.
%
% Usage:
%   EEG        = ieeglab_preprocess(EEG)        % opens the GUI
%   [EEG, com] = ieeglab_preprocess(EEG, opt)   % headless, no dialog
%
% Passing opt skips the dialog entirely, which is what makes the pipeline
% scriptable and testable. Any field left unset falls back to the same default
% the GUI would have shown. See ieeglab_gui_preprocess for the full list; the
% commonly used ones are:
%
%   .apply_highpass / .highpass        logical / Hz
%   .apply_notch    / .notch          logical / Hz (scalar or vector)
%   .apply_lowpass  / .lowpass        logical / Hz
%   .filter_type                      1 = noncausal zero-phase, 2 = minimum phase
%   .apply_epoch    / .epoch_window   logical / [start stop] ms
%   .apply_car      / .car_method     logical / 'carla' | 'varsubset' | 'car'
%   .apply_baseline / .baseline_period logical / [start stop] ms
%   .remove_rare_cond / .min_trials   logical / integer
%   .remove_no_coords                 logical
%   .plot                             logical, draw before/after figures.
%                                     Defaults to false in the headless path.
%
% Example (headless):
%   EEG = ieeglab_preprocess(EEG, struct('apply_highpass',true,'highpass',0.5, ...
%           'apply_epoch',true,'epoch_window',[-500 1000], ...
%           'apply_car',true,'car_method','carla','apply_baseline',true));
%
% Cedric Cannard, iEEGLAB, 2025-2026

com = '';
interactive = (nargin < 2 || isempty(opt));

if interactive
    % GUI to get user choices
    [EEG, wasCancelled] = ieeglab_gui_preprocess(EEG);
    if wasCancelled
        return  % user aborted, exit gracefully
    end
    opt = EEG.ieeglab.opt;
else
    % Merge the supplied options over anything already on the dataset, so a
    % scripted call can override just the fields it cares about.
    if isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'opt') && isstruct(EEG.ieeglab.opt)
        base = EEG.ieeglab.opt;
        f = fieldnames(opt);
        for ii = 1:numel(f), base.(f{ii}) = opt.(f{ii}); end
        opt = base;
    end
    if ~isfield(opt,'plot') || isempty(opt.plot), opt.plot = false; end
    EEG.ieeglab.opt = opt;
end
if ~isfield(opt,'plot') || isempty(opt.plot), opt.plot = interactive; end
if ~isfield(opt,'verbose') || isempty(opt.verbose), opt.verbose = true; end

% Accept the pre-1.0 field names for re-referencing
if isfield(opt,'apply_acar') && ~isfield(opt,'apply_car'), opt.apply_car = opt.apply_acar; end

% % --- Preconditions ---
% if ~isfield(opt,'events') || ~istable(opt.events) || isempty(opt.events) ...
%         || ~isfield(opt,'event_filters') || isempty(fieldnames(opt.event_filters))
%     % nothing to do
%     return;
% end

% Event selection. Filters name either 'type' (the EEGLAB event type) or a
% column of the BIDS events table; both branches now actually delete the
% non-matching events, and opt.events is kept row-aligned with EEG.event.
if isfield(opt,'event_filters') && isstruct(opt.event_filters) && ~isempty(fieldnames(opt.event_filters))

    ev_choices = opt.event_filters;
    vars       = fieldnames(ev_choices);
    haveTbl    = isfield(opt,'events') && istable(opt.events) && ~isempty(opt.events);

    if haveTbl && height(opt.events) ~= numel(EEG.event)
        warning('ieeglab_preprocess:eventMismatch', ...
            ['The BIDS events table has %d rows but the dataset has %d events, so table ' ...
             'columns cannot be used for filtering. Filtering on event type only.'], ...
            height(opt.events), numel(EEG.event));
        haveTbl = false;
    end

    for iVar = 1:numel(vars)
        varName = vars{iVar};
        valsToKeep = ev_choices.(varName);
        if isempty(valsToKeep), continue; end
        if isempty(EEG.event), break; end

        if strcmpi(varName, 'type') || strcmpi(varName, 'var_type')
            ev_values = string({EEG.event.type});
        elseif haveTbl && ismember(varName, opt.events.Properties.VariableNames)
            ev_values = string(opt.events.(varName));
        else
            warning('ieeglab_preprocess:unknownFilterField', ...
                'Event filter field "%s" is neither the event type nor a column of the events table; ignoring it.', varName);
            continue
        end

        % string/ismember rather than ismissing: type-safe for numeric columns,
        % which previously threw when the GUI stringified the choice list.
        toRemove = ~ismember(ev_values(:)', string(valsToKeep(:))');
        if ~any(toRemove), continue; end

        fprintf('Removing %d/%d events whose "%s" is not one of: %s\n', ...
            sum(toRemove), numel(toRemove), varName, strjoin(cellstr(string(valsToKeep(:))'), ', '));
        EEG.event(toRemove) = [];
        if haveTbl
            opt.events(toRemove,:) = [];
        end
        EEG = eeg_checkset(EEG, 'eventconsistency');
    end
    EEG.ieeglab.opt = opt;   % persist, so a second run sees the aligned table
end

% %  Drop heavy event table from options (to save memory) 
% if isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'opt') && isfield(EEG.ieeglab.opt,'events')
%     EEG.ieeglab.opt = rmfield(EEG.ieeglab.opt, 'events');
% end


% Remove electrodes with no coordinates (optional)
if isfield(opt, 'remove_no_coords') && opt.remove_no_coords && isfield(EEG,'chanlocs') && ~isempty(EEG.chanlocs)
    hasXYZ = arrayfun(@(c) isfield(c,'X') && isfield(c,'Y') && isfield(c,'Z') && ...
        ~isempty(c.X) && ~isempty(c.Y) && ~isempty(c.Z) && ...
        all(isfinite([c.X c.Y c.Z])), EEG.chanlocs);

    removed_elecs = {EEG.chanlocs(~hasXYZ).labels};  % for logging
    if any(~hasXYZ)
        warning('Removing %d channels with no 3D (XYZ) coordinates: %s\n', ...
            nnz(~hasXYZ), strjoin(removed_elecs, ', '));
        keepChanIdx = find(hasXYZ);
        EEG = pop_select(EEG, 'channel', keepChanIdx);
        EEG = eeg_checkset(EEG);
    else
        removed_elecs = {};
        disp("All electrodes have 3D (XYZ) coordinates.")
    end

    % Remove events that reference a removed electrode (CCEP stimulation sites).
    % Guarded: continuous datasets legitimately have no events at all.
    if ~isempty(removed_elecs) && isfield(EEG,'event') && ~isempty(EEG.event)
        % Exact token match: contains() also matched 'RA10-RA9' when removing 'RA1'
        trials_to_rem = ieeglab_events_using(EEG.event, removed_elecs);
        if any(trials_to_rem)
            fprintf('Removing %d/%d events that reference an electrode with no 3D coordinates.\n', ...
                sum(trials_to_rem), numel(trials_to_rem));
            EEG.event(trials_to_rem) = [];
            if isfield(opt,'events') && istable(opt.events) && height(opt.events) == numel(trials_to_rem)
                opt.events(trials_to_rem,:) = [];
            end
            EEG = eeg_checkset(EEG, 'eventconsistency');
            EEG = eeg_checkset(EEG);
            EEG.ieeglab.opt = opt;
        end
    end
end

% Bad channels: clinician labels (BIDS channels.tsv status), seizure-zone labels,
% an explicit list, or automatic detection. ieeglab_load MARKS them; this step
% removes them when asked. Channels left marked are still kept out of the CAR
% reference and out of the statistics.
doBad = local_opt(opt,'remove_bad_channels',false) || local_opt(opt,'exclude_soz',false) || ...
        local_opt(opt,'exclude_irritative',false) || local_opt(opt,'auto_bad_channels',false) || ...
        ~isempty(local_opt(opt,'bad_channels',[]));
if doBad
    if local_opt(opt,'remove_bad_channels',false), act = 'remove'; else, act = 'mark'; end
    EEG = ieeglab_bad_channels(EEG, struct( ...
        'channels_tsv',       local_opt(opt,'channels_tsv',''), ...
        'elec_tsv',           local_opt(opt,'elec_tsv',''), ...
        'bad_channels',       local_opt(opt,'bad_channels',[]), ...
        'exclude_soz',        local_opt(opt,'exclude_soz',false), ...
        'exclude_irritative', local_opt(opt,'exclude_irritative',false), ...
        'auto_detect',        local_opt(opt,'auto_bad_channels',false), ...
        'action',             act, ...
        'drop_bad_stim_sites',true, ...
        'verbose',            opt.verbose));
    if isfield(opt,'events') && istable(opt.events) && isfield(EEG.ieeglab,'opt') && isfield(EEG.ieeglab.opt,'events')
        opt.events = EEG.ieeglab.opt.events;     % keep the local copy aligned
    end
end

% No events is a legitimate state: the README advertises continuous mode for
% epilepsy and clinical monitoring. Filtering still applies; the event-dependent
% steps are skipped rather than treated as an error.
continuousMode = ~isfield(EEG,'event') || isempty(EEG.event);
if continuousMode
    if isfield(opt,'apply_epoch') && opt.apply_epoch || ...
       isfield(opt,'apply_car') && opt.apply_car || ...
       isfield(opt,'apply_baseline') && opt.apply_baseline
        warning('ieeglab_preprocess:continuousMode', ...
            ['No events in the dataset - running in CONTINUOUS mode. Filtering and ' ...
             'resampling will be applied; epoching, re-referencing and baseline ' ...
             'correction are skipped because they need events.']);
    end
    opt.apply_epoch    = false;
    opt.apply_car      = false;
    opt.apply_acar     = false;
    opt.apply_baseline = false;
end

% Stimulation-artifact blanking, BEFORE any filtering. A zero-phase FIR rings
% symmetrically about the artifact step, smearing it forward into the early
% response and backward into the baseline (issue #10). Blanking first means the
% filter never sees the discontinuity.
if isfield(opt,'apply_blank') && opt.apply_blank && ~continuousMode
    if strcmp(ieeglab_detect_mode(EEG), 'ccep')
        EEG = ieeglab_blank_stim(EEG, opt);
    elseif opt.verbose
        fprintf('[preprocess] Skipping stimulation blanking: not CCEP data.\n');
    end
end

% Downsample. Gated on the GUI's own checkbox (apply_ds); previously the flag
% was ignored, so unchecking "Downsample" still resampled the data.
do_ds = (~isfield(opt,'apply_ds') || isempty(opt.apply_ds) || opt.apply_ds);
if do_ds && isfield(opt, 'downsample') && ~isempty(opt.downsample) && opt.downsample<EEG.srate
    fprintf("Downsampling iEEG data to %g Hz... \n", opt.downsample)
    EEG = pop_resample(EEG, opt.downsample);
end

% Global filter type -> minphase flag for pop_eegfiltnew
minphase = false;  % default = noncausal zero-phase
if isfield(opt,'filter_type') && ~isempty(opt.filter_type)
    idx = round(double(opt.filter_type));
    minphase = (idx == 2);   % 1=noncausal, 2=minimum-phase
elseif isfield(opt,'filter_type_label') && ~isempty(opt.filter_type_label)
    lbl = lower(string(opt.filter_type_label));
    if contains(lbl,'minimum') || contains(lbl,'causal')
        minphase = true;
    elseif contains(lbl,'noncausal') || contains(lbl,'zero')
        minphase = false;
    end
end

% High-pass filter
if isfield(opt,'apply_highpass') && opt.apply_highpass && isfield(opt,'highpass') ...
        && ~isempty(opt.highpass) && opt.highpass > 0
    EEG = pop_eegfiltnew(EEG, 'locutoff', double(opt.highpass), 'usefftfilt', 1, 'minphase', minphase);
end

% Notch filter
if isfield(opt,'apply_notch') && opt.apply_notch && isfield(opt,'notch') && ~isempty(opt.notch)
    nyq = EEG.srate/2;
    centers = double(opt.notch(:))';
    BW = 2;                     % total bandwidth (Hz)
    for f0 = centers
        if ~isfinite(f0) || f0<=0 || f0>=nyq, continue; end
        lo = max(0, f0 - BW/2);
        hi = min(nyq-1e-6, f0 + BW/2);
        if hi <= lo, continue; end
        EEG = pop_eegfiltnew(EEG,'locutoff',lo,'hicutoff',hi,'usefftfilt',1,'revfilt',1, 'minphase', minphase);
    end
end

% Low-pass filter
if isfield(opt,'apply_lowpass') && opt.apply_lowpass && isfield(opt,'lowpass') ...
        && ~isempty(opt.lowpass) && opt.lowpass > 0
    lp = double(opt.lowpass);
    nyq = EEG.srate/2;
    if isfinite(lp) && lp > 0 && lp < nyq
        EEG = pop_eegfiltnew(EEG,'hicutoff',lp,'usefftfilt',1, 'minphase', minphase);
    end
end

% Remove conditions with too few trials
if isfield(opt,'remove_rare_cond') && opt.remove_rare_cond && isfield(opt,'min_trials')
    minN = max(0, round(double(opt.min_trials)));
    if minN>0 && isfield(EEG,'event') && ~isempty(EEG.event)
        types = string({EEG.event.type});
        u = unique(types);
        cnt = arrayfun(@(x) sum(types==x), u);
        rm = u(cnt < minN);
        if ~isempty(rm)
            fprintf('Removing %d condition(s) with < %d trials: %s\n', numel(rm), minN, strjoin(cellstr(rm), ', '));
            keep = ~ismember(types, rm);
            EEG.event = EEG.event(keep);
            try EEG = eeg_checkset(EEG,'makeur'); catch, end
        else
            fprintf('No conditions below %d trials.\n', minN);
        end
    end
end


% Epoching (uses ALL current event types if none chosen in GUI)
if isfield(opt,'apply_epoch') && opt.apply_epoch && ...
        isfield(opt,'epoch_window') && ...
        isfield(EEG,'event') && ~isempty(EEG.event)

    % window
    t_ms = double(opt.epoch_window(:))';
    if numel(t_ms)<2 || ~all(isfinite(t_ms(1:2))) || t_ms(2)<=t_ms(1)
        warning('Invalid epoch window; skipping epoching.');
    else
        % collect unique event types as strings
        if istable(EEG.event), evTypes = EEG.event.type; else, evTypes = {EEG.event.type}'; end
        if isnumeric(evTypes), evTypes = string(evTypes); end
        if iscell(evTypes),   evTypes = string(evTypes); end
        evtTypes = unique(cellstr(string(evTypes(:))));  % cellstr

        if isempty(evtTypes)
            warning('No event types found to epoch around; skipping.');
        else
            fprintf('Epoching around %d event types, window [%g %g] ms\n', numel(evtTypes), t_ms(1), t_ms(2));
            EEG = pop_epoch(EEG, evtTypes, t_ms/1000, 'epochinfo','yes', 'newname','iEEGLAB epochs');
            EEG = eeg_checkset(EEG);
        end
    end

    % event_filters were only for selection; drop from opts
    if isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'opt') && isfield(EEG.ieeglab.opt,'event_filters')
        EEG.ieeglab.opt = rmfield(EEG.ieeglab.opt, 'event_filters');
    end
end


% Re-referencing (CARLA by default; see ieeglab_car for the alternatives)
if isfield(opt,'apply_car') && opt.apply_car && isfield(EEG,'trials') && EEG.trials > 1

    if opt.plot
        respBefore = local_response_trace(EEG);
    end

    EEG = ieeglab_car(EEG, opt);

    if opt.plot
        respAfter = local_response_trace(EEG);
        local_plot_car_comparison(EEG.times(:), respBefore, respAfter, EEG.ref);
    end
end


% Baseline correction
if isfield(opt,'apply_baseline') && opt.apply_baseline
    EEG = ieeglab_rm_baseline(EEG);
end

% Persist the (possibly modified) options and report success to eeglab_new.
EEG.ieeglab.opt = opt;
EEG = eeg_checkset(EEG);
com = 'EEG = ieeglab_preprocess(EEG);';

end

% ========================== local helpers ==========================

function resp = local_response_trace(EEG)
% Mean response trace used for the before/after re-referencing figure. For CCEP
% data we follow the channel named after the '-' in the event type; for other
% designs we take a trimmed mean across channels. Trials whose target label does
% not match any channel are left as NaN instead of raising a size error.
nTr  = size(EEG.data, 3);
resp = nan(numel(EEG.times), nTr);
labels = {EEG.chanlocs.labels};
for iTrial = 1:nTr
    ty = '';
    if isfield(EEG,'epoch') && numel(EEG.epoch) >= iTrial && isfield(EEG.epoch,'eventtype')
        ty = EEG.epoch(iTrial).eventtype;
        if iscell(ty) && ~isempty(ty), ty = ty{1}; end
    end
    if ischar(ty) && contains(ty, '-')
        respIdx = find(strcmpi(labels, extractAfter(ty, '-')), 1);
        if isempty(respIdx), continue; end          % label not in this montage
        resp(:,iTrial) = squeeze(EEG.data(respIdx,:,iTrial));
    else
        resp(:,iTrial) = squeeze(trimmean(EEG.data(:,:,iTrial), 20, 1));
    end
end
end

function local_plot_car_comparison(x, before, after, refName)
mu = @(R) trimmean(R, 20, 2);
sem = @(R) std(R, 0, 2, 'omitmissing') ./ sqrt(max(sum(~isnan(R), 2), 1));
figure('color','w'); hold on
ax = gca; co = ax.ColorOrder;
band = @(m,s,c) fill([x; flipud(x)]', [m+s; flipud(m-s)]', c, 'FaceAlpha',0.2, 'EdgeColor', c);
m1 = mu(before); m2 = mu(after);
band(m1, sem(before), co(1,:));
h1 = plot(x, m1, 'LineWidth',2, 'Color', co(1,:), 'DisplayName','Before re-referencing');
band(m2, sem(after), co(2,:));
h2 = plot(x, m2, 'LineWidth',2, 'Color', co(2,:), 'DisplayName', sprintf('After %s', refName));
legend([h1 h2], 'Location','best'); box on
xlabel('Time (ms)'); ylabel('Amplitude (\muV)');
% Concatenated rather than sprintf'd: sprintf would warn on the TeX '\pm'.
title(['Mean \pm 1 SEM (before vs after ' char(refName) ')']);
end

function v = local_opt(opt, name, default)
% Field of opt, or default when absent or empty.
if isfield(opt, name) && ~isempty(opt.(name)), v = opt.(name); else, v = default; end
end
