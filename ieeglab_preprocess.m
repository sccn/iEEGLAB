function EEG = ieeglab_preprocess(EEG)

% GUI to get user choices
if nargin < 2
    [EEG, wasCancelled] = ieeglab_gui_preprocess(EEG);
    if wasCancelled
        return  % user aborted, exit gracefully
    end
end

opt = EEG.ieeglab.opt;

% % --- Preconditions ---
% if ~isfield(opt,'events') || ~istable(opt.events) || isempty(opt.events) ...
%         || ~isfield(opt,'event_filters') || isempty(fieldnames(opt.event_filters))
%     % nothing to do
%     return;
% end

if isfield(opt, 'event_filters') && ~isempty(opt.event_filters)
    ev_tbl0     = opt.events;              % original TSV table (row-aligned to how EEG.event was created)
    ev_choices  = opt.event_filters;
    vars        = fieldnames(ev_choices);
    
    % check EEG.event and opt.events still match length
    if size(ev_tbl0,1) ~= length(EEG.event)
        error("Events from tsv and EEGLAB dataset mismatch! ")
    end
    
    for iVar = 1:length(vars)
        varName = vars{iVar};
        if ~isempty(ev_choices.(varName))
            if strcmpi(varName, 'var_type')
                trialsToKeep = ismissing({EEG.event.type}, ev_choices.(varName));
                warning("Per request, removing %g/%g events that are not of type: ", sum(trialsToKeep), length(EEG.event))
                disp(ev_choices.(varName))
                EEG.event(~trialsToKeep) = [];
            else
                
                ev_values = ev_tbl0.(vars{iVar});
                valsToKeep = ev_choices.(varName);
                trialsToRem = ~ismissing(ev_values, valsToKeep);
                if any(trialsToRem)
                    warning("Removing %g/%g events with field '%s' that do not have value: ", sum(trialsToRem), length(ev_values), vars{iVar})
                    disp(valsToKeep)
                end
            end
        end
    end    
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

    % Remove corresponding events whose TYPE includes any removed label (if events exist)
    trials_to_rem = contains({EEG.event.type}, removed_elecs);
    if any(trials_to_rem)
        warning("Keeping %g/%g events containing electrodes that had no 3D coordinates.", sum(trials_to_rem), length({EEG.event.type}))
        EEG.event(trials_to_rem) = [];
        opt.events(trials_to_rem,:) = [];
        EEG = eeg_checkset(EEG, 'eventconsistency');
        EEG = eeg_checkset(EEG);
    end
end

% Sanity check that we still have some events left
if isempty(EEG.event)
    error("No events left after event filtering!")
end

% Downsample
if isfield(opt, 'downsample') && ~isempty(opt.downsample) && opt.downsample<EEG.srate
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


% Apply aCAR (Huang et al., 2024) for epoched data
if isfield(opt,'apply_acar') && opt.apply_acar && isfield(EEG,'trials') && EEG.trials>1

    % PLOT BEFORE CAR (OPTIONAL)
    respData = nan(numel(EEG.times), numel(EEG.epoch));
    for iTrial = 1:numel(EEG.epoch)
        if ~isnumeric(EEG.event(1).type) && contains(EEG.event(1).type,'-') % CCEP: pick target chan from event
            respElec   = extractAfter(EEG.event(iTrial).type,'-');
            respIdx    = strcmpi({EEG.chanlocs.labels}, respElec);
            respData(:,iTrial) = squeeze(EEG.data(respIdx,:,iTrial));
        else
            respData(:,iTrial) = squeeze(trimmean(EEG.data(:,:,iTrial),20,1));
        end
    end
    mu1  = trimmean(respData,20,2);
    n1   = sum(~isnan(respData),2);
    sem1 = std(respData,0,2,'omitmissing') ./ sqrt(max(n1,1));
    figure('color','w'); hold on 
    ax = gca; co = ax.ColorOrder; 
    col1 = co(1,:); col2 = co(2,:);
    x = EEG.times(:);
    fill([x; flipud(x)]', [mu1+sem1; flipud(mu1-sem1)]', col1, 'FaceAlpha',0.2, 'EdgeColor', col1);
    h1 = plot(x, mu1, 'LineWidth',2, 'Color', col1, 'DisplayName','Before CAR');

    % APPLY CAR 
    EEG = ieeglab_car(EEG);

    % PLOT AFTER CAR (OPTIONAL)
    respData = nan(numel(EEG.times), numel(EEG.epoch));
    for iTrial = 1:numel(EEG.epoch)
        if ~isnumeric(EEG.event(1).type) && contains(EEG.event(1).type,'-')
            respElec   = extractAfter(EEG.event(iTrial).type,'-');
            respIdx    = strcmpi({EEG.chanlocs.labels}, respElec);
            respData(:,iTrial) = squeeze(EEG.data(respIdx,:,iTrial));
        else
            respData(:,iTrial) = squeeze(trimmean(EEG.data(:,:,iTrial),20,1));
        end
    end
    mu2  = trimmean(respData,20,2);
    n2   = sum(~isnan(respData),2);
    sem2 = std(respData,0,2,'omitmissing') ./ sqrt(max(n2,1));
    fill([x; flipud(x)]', [mu2+sem2; flipud(mu2-sem2)]', col2, 'FaceAlpha',0.2, 'EdgeColor', col2);
    h2 = plot(x, mu2, 'LineWidth',2, 'Color', col2, 'DisplayName','After adjusted CAR');
    legend([h1 h2], 'Location','best'); box on
    xlabel('Time (ms)'); ylabel('Amplitude (\muV)');
    title('Mean \pm 1 SEM (Before vs After Ajusted CAR)');
end


% Baseline correction
if isfield(opt,'apply_baseline') && opt.apply_baseline
    EEG = ieeglab_rm_baseline(EEG);
end

end



