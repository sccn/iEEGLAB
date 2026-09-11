function [EEG, com] = ieeglab_preprocess(EEG, opt)
% ieeglab_preprocess() - Preprocess iEEG data: event selection, resampling,
%                        filtering, epoching, re-referencing, baseline removal.
%
% Usage:
%   EEG        = ieeglab_preprocess(EEG)        % opens the GUI
%   [EEG, com] = ieeglab_preprocess(EEG, opt)   % headless, no dialog
%
% Passing opt skips the dialog entirely, which is what makes the pipeline
% scriptable and testable. The commonly used fields are:
%
%   .apply_highpass / .highpass        logical / Hz
%   .apply_notch    / .notch          logical / Hz (scalar or vector)
%   .apply_lowpass  / .lowpass        logical / Hz
%   .filter_type                      1 = noncausal zero-phase, 2 = minimum phase
%   .apply_ds       / .downsample     logical / Hz
%   .event_filters                    struct: field = event 'type' or an events.tsv
%                                     column, value = the values to keep
%   .apply_epoch    / .epoch_window   logical / [start stop] ms
%   .apply_car      / .car_method     logical / 'carla' | 'varsubset' | 'car'
%   .apply_baseline / .baseline_period logical / [start stop] ms
%   .remove_rare_cond / .min_trials   logical / integer, counted per stimulation site
%   .remove_no_coords                 logical
%   .remove_bad_channels              remove clinician-marked channels (else they are
%                                     only marked, and kept out of reference and stats)
%   .exclude_soz / .exclude_irritative  remove seizure-onset / irritative-zone contacts
%   .bad_channels                     labels or indices (of the dataset passed in)
%   .apply_blank / .blank_window      stimulation-artifact blanking (CCEP)
%   .reject_trials / .reject_z        automatic outlier-trial rejection (ieeglab_reject_trials)
%   .plot                             draw before/after figures. Default false headless.
%
% Which steps run is decided by the call: in a headless call, step switches
% (apply_*, remove_*, exclude_*, event_filters, bad_channels, reject_trials,
% downsample, plot) come only from opt. Parameters such as file paths and
% filter settings saved on the dataset by an earlier run are reused, but a
% saved switch never re-applies its step - calling ieeglab_preprocess again to
% change the baseline does not high-pass or re-reference the data a second time.
%
% The returned com (EEGLAB history) carries every option used, so it replays
% the same run headlessly.
%
% Example (headless):
%   EEG = ieeglab_preprocess(EEG, struct('apply_highpass',true,'highpass',0.5, ...
%           'apply_epoch',true,'epoch_window',[-500 1000], ...
%           'apply_car',true,'car_method','carla','apply_baseline',true));
%
% Cedric Cannard, iEEGLAB, 2025-2026

com = '';
interactive = (nargin < 2 || isempty(opt));

stepFields = {'apply_ds','downsample','apply_highpass','apply_notch','apply_lowpass', ...
    'apply_epoch','apply_car','apply_acar','apply_baseline','apply_blank', ...
    'remove_rare_cond','remove_no_coords','remove_bad_channels','exclude_soz', ...
    'exclude_irritative','auto_bad_channels','bad_channels','reject_trials', ...
    'event_filters','plot'};

stored = struct();
if isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'opt') && isstruct(EEG.ieeglab.opt)
    stored = local_aliases(EEG.ieeglab.opt);
    stored = rmfield(stored, intersect(fieldnames(stored), stepFields));
end

if interactive
    % The dialog starts from its own defaults; hand it only the saved parameters.
    EEGd = EEG;
    EEGd.ieeglab.opt = stored;
    [EEGd, wasCancelled] = ieeglab_gui_preprocess(EEGd);
    if wasCancelled
        return  % user aborted, exit gracefully
    end
    EEG = EEGd;
    opt = local_aliases(EEG.ieeglab.opt);
    if ~isfield(opt,'plot') || isempty(opt.plot), opt.plot = true; end
else
    base = stored;
    userOpt = local_aliases(opt);
    f = fieldnames(userOpt);
    for ii = 1:numel(f), base.(f{ii}) = userOpt.(f{ii}); end
    opt = base;
    if ~isfield(opt,'plot') || isempty(opt.plot), opt.plot = false; end
end
if ~isfield(opt,'verbose') || isempty(opt.verbose), opt.verbose = true; end

% Channel indices refer to the dataset as passed in. Resolve them to labels now,
% before any step removes a channel and shifts the indices.
opt.bad_channels = local_index_to_labels(local_opt(opt,'bad_channels',[]), {EEG.chanlocs.labels});
EEG.ieeglab.opt = local_persist(opt);

% ------------------------------------------------------------------
% Event selection. Filters name either 'type' (the EEGLAB event type) or a
% column of the BIDS events table, linked to each event by its tsv_row.
% 'boundary' markers are never removed: pop_epoch needs them to reject
% epochs that span removed data.
% ------------------------------------------------------------------
if isfield(opt,'event_filters') && isstruct(opt.event_filters) && ~isempty(fieldnames(opt.event_filters)) ...
        && isfield(EEG,'event') && ~isempty(EEG.event)

    ev_choices = opt.event_filters;
    vars       = fieldnames(ev_choices);
    haveTbl    = isfield(opt,'events') && istable(opt.events) && ~isempty(opt.events);
    nReal0     = nnz(~strcmpi(local_types(EEG), 'boundary'));

    for iVar = 1:numel(vars)
        varName = vars{iVar};
        valsToKeep = ev_choices.(varName);
        if ischar(valsToKeep) || isstring(valsToKeep), valsToKeep = cellstr(valsToKeep); end
        if isnumeric(valsToKeep) || islogical(valsToKeep), valsToKeep = num2cell(valsToKeep); end
        if isempty(valsToKeep), continue; end
        if isempty(EEG.event), break; end

        types = local_types(EEG);
        isB = strcmpi(types, 'boundary');
        if strcmpi(varName, 'type') || strcmpi(varName, 'var_type')
            ev_values = types;
        elseif haveTbl && ismember(varName, opt.events.Properties.VariableNames)
            ev_values = local_table_values(EEG, opt.events, varName);
            if all(ismissing(ev_values(~isB)))
                warning('ieeglab_preprocess:eventLink', ...
                    ['Cannot link the events table to the dataset''s events (no tsv_row and different ' ...
                     'row counts), so the "%s" filter is ignored. Reload with ieeglab_load.'], varName);
                continue
            end
        else
            warning('ieeglab_preprocess:unknownFilterField', ...
                'Event filter field "%s" is neither the event type nor a column of the events table; ignoring it.', varName);
            continue
        end

        keepStr = string(cellfun(@(v) char(string(v)), valsToKeep(:)', 'UniformOutput', false));
        toRemove = ~ismember(ev_values(:)', keepStr) & ~isB;
        if ~any(toRemove), continue; end

        fprintf('Removing %d/%d events whose "%s" is not one of: %s\n', ...
            sum(toRemove), nnz(~isB), varName, strjoin(cellstr(keepStr), ', '));
        EEG.event(toRemove) = [];
        EEG = eeg_checkset(EEG, 'eventconsistency');
        opt = local_sync_events(EEG, opt, toRemove);
    end
    if nReal0 > 0 && (isempty(EEG.event) || all(strcmpi(local_types(EEG), 'boundary')))
        error('ieeglab_preprocess:allEventsFiltered', ...
            'The event filters removed every event (%s). Check the values to keep.', ...
            strjoin(vars', ', '));
    end
end

% ------------------------------------------------------------------
% Remove electrodes with no coordinates (optional)
% ------------------------------------------------------------------
if local_opt(opt,'remove_no_coords',false) && isfield(EEG,'chanlocs') && ~isempty(EEG.chanlocs)
    hasXYZ = arrayfun(@(c) isfield(c,'X') && isfield(c,'Y') && isfield(c,'Z') && ...
        ~isempty(c.X) && ~isempty(c.Y) && ~isempty(c.Z) && ...
        all(isfinite([c.X c.Y c.Z])), EEG.chanlocs);

    removed_elecs = {EEG.chanlocs(~hasXYZ).labels};
    if all(~hasXYZ)
        warning('ieeglab_preprocess:noCoordinates', ...
            'No channel has 3D coordinates; remove_no_coords is ignored.');
        removed_elecs = {};
    elseif any(~hasXYZ)
        warning('ieeglab_preprocess:removedNoCoords', ...
            'Removing %d channels with no 3D (XYZ) coordinates: %s', ...
            nnz(~hasXYZ), strjoin(removed_elecs, ', '));
        EEG = pop_select(EEG, 'channel', find(hasXYZ));
        EEG = eeg_checkset(EEG);
        prev = {};
        if isfield(EEG.ieeglab,'removed_channels'), prev = EEG.ieeglab.removed_channels; end
        EEG.ieeglab.removed_channels = unique([prev(:); removed_elecs(:)], 'stable');
    elseif opt.verbose
        disp("All electrodes have 3D (XYZ) coordinates.")
    end

    % Remove events that stimulate a removed electrode (exact token match).
    if ~isempty(removed_elecs) && isfield(EEG,'event') && ~isempty(EEG.event)
        trials_to_rem = ieeglab_events_using(EEG.event, removed_elecs);
        if any(trials_to_rem)
            fprintf('Removing %d/%d events that reference an electrode with no 3D coordinates.\n', ...
                sum(trials_to_rem), numel(trials_to_rem));
            EEG.event(trials_to_rem) = [];
            EEG = eeg_checkset(EEG, 'eventconsistency');
            opt = local_sync_events(EEG, opt, trials_to_rem);
        end
    end
end

% ------------------------------------------------------------------
% Bad channels: clinician labels (BIDS channels.tsv status), seizure-zone
% labels, an explicit list, or automatic detection. ieeglab_load MARKS them.
% Whenever any channel is marked bad - whatever option marked it - the trials
% that stimulate it are dropped, so the analysed trials never depend on an
% unrelated switch. remove_bad_channels removes every marked channel;
% exclude_soz / exclude_irritative remove those contacts.
% ------------------------------------------------------------------
removeBad = local_opt(opt,'remove_bad_channels',false);
exclSOZ   = local_opt(opt,'exclude_soz',false);
exclIRR   = local_opt(opt,'exclude_irritative',false);
anyMarked = isfield(EEG.chanlocs,'status') && ...
    any(cellfun(@(x) ~isempty(x) && strcmpi(char(x),'bad'), {EEG.chanlocs.status}));
doBad = removeBad || exclSOZ || exclIRR || local_opt(opt,'auto_bad_channels',false) || ...
        ~isempty(opt.bad_channels) || anyMarked;
if doBad
    if removeBad, act = 'remove'; else, act = 'mark'; end
    EEG = ieeglab_bad_channels(EEG, struct( ...
        'channels_tsv',       local_opt(opt,'channels_tsv',''), ...
        'elec_tsv',           local_opt(opt,'elec_tsv',''), ...
        'bad_channels',       {opt.bad_channels}, ...
        'exclude_soz',        exclSOZ, ...
        'exclude_irritative', exclIRR, ...
        'auto_detect',        local_opt(opt,'auto_bad_channels',false), ...
        'action',             act, ...
        'drop_bad_stim_sites',true, ...
        'verbose',            opt.verbose));
    if ~removeBad && (exclSOZ || exclIRR)
        Tb = EEG.ieeglab.bad_channels;
        zoneLabels = cellstr(Tb.label(Tb.source == "seizure_zone"));
        rm = ismember({EEG.chanlocs.labels}, zoneLabels);
        if any(rm)
            fprintf('[preprocess] Removing %d seizure-zone contact(s): %s\n', nnz(rm), strjoin(zoneLabels, ', '));
            EEG = pop_select(EEG, 'nochannel', find(rm));
            prev = {};
            if isfield(EEG.ieeglab,'removed_channels'), prev = EEG.ieeglab.removed_channels; end
            EEG.ieeglab.removed_channels = unique([prev(:); zoneLabels(:)], 'stable');
        end
    end
    opt = local_sync_events(EEG, opt, []);
end

% No events is a legitimate state (continuous mode for epilepsy and clinical
% monitoring). Filtering still applies; the event-dependent steps are skipped.
continuousMode = ~isfield(EEG,'event') || isempty(EEG.event) || all(strcmpi(local_types(EEG), 'boundary'));
if continuousMode
    if local_opt(opt,'apply_epoch',false) || local_opt(opt,'apply_car',false) || local_opt(opt,'apply_baseline',false)
        warning('ieeglab_preprocess:continuousMode', ...
            ['No events in the dataset - running in CONTINUOUS mode. Filtering and ' ...
             'resampling will be applied; epoching, re-referencing and baseline ' ...
             'correction are skipped because they need events.']);
    end
    opt.apply_epoch    = false;
    opt.apply_car      = false;
    opt.apply_baseline = false;
end

% Global filter type -> minphase flag for pop_eegfiltnew
minphase = false;  % default = noncausal zero-phase
if isfield(opt,'filter_type') && ~isempty(opt.filter_type)
    minphase = (round(double(opt.filter_type)) == 2);   % 1=noncausal, 2=minimum-phase
elseif isfield(opt,'filter_type_label') && ~isempty(opt.filter_type_label)
    lbl = lower(string(opt.filter_type_label));
    minphase = contains(lbl,'minimum') || contains(lbl,'causal') && ~contains(lbl,'noncausal');
end
doHP = local_opt(opt,'apply_highpass',false) && local_opt(opt,'highpass',0) > 0;
doNotch = local_opt(opt,'apply_notch',false) && ~isempty(local_opt(opt,'notch',[]));
doLP = local_opt(opt,'apply_lowpass',false) && local_opt(opt,'lowpass',0) > 0;
do_ds = local_opt(opt,'apply_ds', isfield(opt,'downsample')) && ...
        ~isempty(local_opt(opt,'downsample',[])) && opt.downsample < EEG.srate;

% ------------------------------------------------------------------
% Stimulation-artifact blanking, BEFORE any filtering. A zero-phase FIR rings
% symmetrically about the artifact step, smearing it forward into the early
% response and backward into the baseline (issue #10). Blanking first means
% the filter never sees the discontinuity. blank_method 'nan' cannot go
% through a filter or the resampler (NaN would spread over the recording), so
% the gap is bridged linearly for them and set back to NaN afterwards.
% ------------------------------------------------------------------
reNaN = false;
if local_opt(opt,'apply_blank',false) && ~continuousMode
    if strcmp(ieeglab_detect_mode(EEG), 'ccep')
        bopt = opt;
        if strcmpi(local_opt(opt,'blank_method','pchip'), 'nan') && (doHP || doNotch || doLP || do_ds)
            bopt.blank_method = 'linear';
            reNaN = true;
        end
        EEG = ieeglab_blank_stim(EEG, bopt);
    elseif opt.verbose
        fprintf('[preprocess] Skipping stimulation blanking: not CCEP data.\n');
    end
end

% Downsample
if do_ds
    fprintf("Downsampling iEEG data to %g Hz... \n", opt.downsample)
    EEG = pop_resample(EEG, opt.downsample);
end

% High-pass filter
if doHP
    EEG = pop_eegfiltnew(EEG, 'locutoff', double(opt.highpass), 'usefftfilt', 1, 'minphase', minphase);
end

% Notch filter
if doNotch
    nyq = EEG.srate/2;
    BW = 2;                     % total bandwidth (Hz)
    for f0 = double(opt.notch(:))'
        if ~isfinite(f0) || f0<=0 || f0>=nyq, continue; end
        lo = max(0, f0 - BW/2);
        hi = min(nyq-1e-6, f0 + BW/2);
        if hi <= lo, continue; end
        EEG = pop_eegfiltnew(EEG,'locutoff',lo,'hicutoff',hi,'usefftfilt',1,'revfilt',1, 'minphase', minphase);
    end
end

% Low-pass filter
if doLP
    lp = double(opt.lowpass);
    if isfinite(lp) && lp < EEG.srate/2
        EEG = pop_eegfiltnew(EEG,'hicutoff',lp,'usefftfilt',1, 'minphase', minphase);
    end
end

if reNaN
    nopt = opt; nopt.blank_method = 'nan'; nopt.verbose = false;
    EEG = ieeglab_blank_stim(EEG, nopt);
    if opt.verbose, fprintf('[blank] Blanked samples set back to NaN after filtering.\n'); end
end

% ------------------------------------------------------------------
% Remove stimulation sites / conditions with too few trials. Counted per
% site as every later step defines it (ROP2-ROP4 and ROP4-ROP2 are one site);
% 'boundary' markers are not a condition.
% ------------------------------------------------------------------
if local_opt(opt,'remove_rare_cond',false) && ~continuousMode
    minN = max(0, round(double(local_opt(opt,'min_trials',0))));
    if minN > 0
        types = local_types(EEG);
        isB = strcmpi(types, 'boundary');
        site = ieeglab_canonical_site(cellstr(types), {EEG.chanlocs.labels});
        [u, ~, g] = unique(site(~isB));
        cnt = accumarray(g(:), 1);
        rmSites = u(cnt < minN);
        if ~isempty(rmSites)
            drop = ismember(site, rmSites) & ~isB;
            fprintf('Removing %d condition(s) with < %d trials: %s\n', numel(rmSites), minN, strjoin(cellstr(rmSites), ', '));
            EEG.event(drop) = [];
            EEG = eeg_checkset(EEG, 'eventconsistency');
            opt = local_sync_events(EEG, opt, drop);
        elseif opt.verbose
            fprintf('No conditions below %d trials.\n', minN);
        end
    end
end

% ------------------------------------------------------------------
% Epoching around every remaining event type except 'boundary'
% ------------------------------------------------------------------
if local_opt(opt,'apply_epoch',false) && isfield(opt,'epoch_window') && ~continuousMode
    t_ms = double(opt.epoch_window(:))';
    if numel(t_ms)<2 || ~all(isfinite(t_ms(1:2))) || t_ms(2)<=t_ms(1)
        warning('ieeglab_preprocess:badEpochWindow', 'Invalid epoch window; skipping epoching.');
    else
        types = local_types(EEG);
        evtTypes = unique(cellstr(types(~strcmpi(types,'boundary'))));
        if isempty(evtTypes)
            warning('ieeglab_preprocess:noEventTypes', 'No event types found to epoch around; skipping.');
        else
            fprintf('Epoching around %d event types, window [%g %g] ms\n', numel(evtTypes), t_ms(1), t_ms(2));
            EEG = pop_epoch(EEG, evtTypes, t_ms/1000, 'epochinfo','yes', 'newname','iEEGLAB epochs');
            EEG = eeg_checkset(EEG);
        end
    end
end

% Automatic outlier-trial rejection (off by default). Before re-referencing, so
% artifact trials cannot distort CARLA's covariance ranking.
if local_opt(opt,'reject_trials',false) && isfield(EEG,'trials') && EEG.trials > 1
    EEG = ieeglab_reject_trials(EEG, struct('z', local_opt(opt,'reject_z',5), ...
        'frac', local_opt(opt,'reject_frac',0.25), 'action', 'remove', 'verbose', opt.verbose));
end

% Re-referencing (CARLA by default; see ieeglab_car for the alternatives).
% Bad channels reach it through chanlocs.status, which survives channel removal.
if local_opt(opt,'apply_car',false) && isfield(EEG,'trials') && EEG.trials > 1
    carOpt = opt;
    carOpt.bad_channels = [];
    if opt.plot
        respBefore = local_response_trace(EEG);
    end

    EEG = ieeglab_car(EEG, carOpt);

    if opt.plot
        respAfter = local_response_trace(EEG);
        local_plot_car_comparison(EEG.times(:), respBefore, respAfter, EEG.ref);
    end
end

% Baseline correction (reads its options from EEG.ieeglab.opt)
if local_opt(opt,'apply_baseline',false)
    EEG.ieeglab.opt = local_persist(opt);
    EEG = ieeglab_rm_baseline(EEG);
end

% Persist the options (never the plot switch) and return a replayable history line.
EEG.ieeglab.opt = local_persist(opt);
if isfield(EEG.ieeglab.opt,'event_filters') && EEG.trials > 1
    EEG.ieeglab.opt = rmfield(EEG.ieeglab.opt, 'event_filters');   % applied; not to be re-applied
end
EEG = eeg_checkset(EEG);
replay = rmfield(opt, intersect(fieldnames(opt), {'events','plot','events_from_tsv'}));
replay.plot = false;
com = sprintf('EEG = ieeglab_preprocess(EEG, %s);', ieeglab_literal(replay));

end

% ========================== local helpers ==========================

function o = local_aliases(o)
% Pre-1.0 names -> current names. A current name, when present, wins.
if ~isstruct(o), o = struct(); return; end
map = {'apply_acar','apply_car'; 'acar_timewin','car_timewin'; 'acar_fraction','car_fraction'};
for k = 1:size(map,1)
    if isfield(o, map{k,1})
        if ~isfield(o, map{k,2}) || isempty(o.(map{k,2}))
            o.(map{k,2}) = o.(map{k,1});
        end
        o = rmfield(o, map{k,1});
    end
end
end

function o = local_persist(o)
if isfield(o,'plot'), o = rmfield(o,'plot'); end
end

function b = local_index_to_labels(b, labels)
if isempty(b) || ~(isnumeric(b) || islogical(b)), return; end
if islogical(b), b = find(b); end
b = double(b(:)');
bad = b < 1 | b > numel(labels) | b ~= round(b);
if any(bad)
    warning('ieeglab_preprocess:badChannelIndex', ...
        'bad_channels indices outside 1..%d are ignored: %s', numel(labels), mat2str(b(bad)));
end
b = labels(b(~bad));
end

function t = local_types(EEG)
t = strings(1, numel(EEG.event));
for i = 1:numel(EEG.event)
    v = EEG.event(i).type;
    if isnumeric(v) || islogical(v), v = num2str(v); end
    t(i) = string(v);
end
end

function v = local_table_values(EEG, T, varName)
% Value of events-table column varName for each event, linked by tsv_row
% (falling back to position when both have the same number of rows).
n = numel(EEG.event);
v = strings(1, n); v(:) = missing;
col = string(T.(varName));
if isfield(EEG.event,'tsv_row') && ismember('tsv_row', T.Properties.VariableNames)
    r = local_rows(EEG);
    [tf, loc] = ismember(r, T.tsv_row);
    v(tf) = col(loc(tf));
elseif height(T) == n
    v = col(:)';
end
end

function r = local_rows(EEG)
r = nan(1, numel(EEG.event));
for i = 1:numel(EEG.event)
    x = EEG.event(i).tsv_row;
    if ~isempty(x) && isnumeric(x), r(i) = double(x(1)); end
end
end

function opt = local_sync_events(EEG, opt, removed)
% Keep opt.events to the rows still represented in EEG.event.
if ~isfield(opt,'events') || ~istable(opt.events), return; end
T = opt.events;
if ismember('tsv_row', T.Properties.VariableNames) && isfield(EEG,'event') && isfield(EEG.event,'tsv_row')
    r = local_rows(EEG);
    r = r(isfinite(r));
    [tf, loc] = ismember(r, T.tsv_row);
    opt.events = T(loc(tf), :);
elseif ~isempty(removed) && height(T) == numel(removed)
    opt.events(removed, :) = [];
end
end

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
