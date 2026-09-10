function [EEG, com] = ieeglab_load(EEG)
% ieeglab_load() - Load BIDS electrode coordinates and events into an EEGLAB
%                  iEEG dataset, and select the channels/events of interest.
%
% Usage:
%   [EEG, com] = ieeglab_load(EEG)
%
% Opens two dialogs: the first picks the electrodes.tsv / events.tsv files and
% general parameters, the second selects channels and event types. Cancelling
% either one returns the dataset unmodified.
%
% com is the EEGLAB history string; it is empty when the user cancels, which is
% what tells eeglab_new not to store a half-processed dataset.
%
% Cedric Cannard, iEEGLAB, 2025-2026

com = '';
EEG_in = EEG;   % restored if the user cancels the second dialog

% check if dataset is already epoched
if ndims(EEG.data) == 3
    % epoched = true;
    % if ~isfield(EEG, 'trials')
    %     EEG.trials = size(EEG.data,3);
    % end
else
    % epoched = false;
    if ~isfield(EEG, 'trials') || isempty(EEG.trials) || EEG.trials == 0
        EEG.trials = 1;
    end
end

% gui_load1 to files to load and general parameters
[opt, wasCanceled] = ieeglab_gui_load1(EEG.filepath);
if wasCanceled || isempty(opt)
    return;  % user aborted, exit gracefully
end

% Pull electrode labels from TSV file
if ~isempty(opt.elec_tsv)
    elecs = readtable(opt.elec_tsv, 'FileType', 'text', 'Delimiter', '\t');
    % opt.elecs = elecs;
    opt.elec_labels = elecs.name;
    opt.elec_labels = strtrim(regexprep(opt.elec_labels, '''', '')); % remove apostrophes inside the labels in some datasets

    % Check the TSV actually has Cartesian (XYZ) coordinate columns.
    % Use Properties.VariableNames, not fieldnames: fieldnames on a table also
    % returns Properties/Row/Variables, and the previous '&&' made the guard
    % unsatisfiable so it never fired in either direction.
    vn = elecs.Properties.VariableNames;
    missingCols = setdiff({'x','y','z'}, lower(vn));
    if ~isempty(missingCols)
        warning('ieeglab_load:noCoordinateColumns', ...
            ['The electrodes file has no %s column(s), so 3D coordinates cannot be loaded.\n' ...
             'Columns found: %s\nContinuing without electrode coordinates.'], ...
             strjoin(upper(missingCols), '/'), strjoin(vn, ', '));
        opt.elec_tsv = '';
    else
        xyz = elecs{:, ismember(lower(vn), {'x','y','z'})};
        if isempty(xyz) || all(~isfinite(xyz(:)))
            warning('ieeglab_load:emptyCoordinates', ...
                'The X/Y/Z columns of the electrodes file are empty or all non-finite. Continuing without electrode coordinates.');
            opt.elec_tsv = '';
        end
    end
end

% Pull events from TSV file or EEGLAB dataset.
% events_from_tsv records which of the two happened: only a BIDS table has an
% 'onset' column in seconds, so only that case may be re-injected into
% EEG.event. Events already on the dataset are left exactly as they are.
opt.events_from_tsv = false;
if ~isempty(opt.events_tsv)
    events = readtable(opt.events_tsv, 'FileType', 'text', 'Delimiter', '\t');
    fprintf("%g total events were imported from the .tsv event file.\n", size(events,1))
    if ~ismember('onset', events.Properties.VariableNames)
        error('ieeglab_load:noOnset', ...
            ['The events file has no "onset" column, which BIDS requires and which is ' ...
             'needed to place events in the recording.\nColumns found: %s'], ...
             strjoin(events.Properties.VariableNames, ', '));
    end
    opt.events = events;
    opt.events_from_tsv = true;

    % sometimes loaded dataset (e.g., mefd) is a portion of the file, so
    % some events no longer match the data --> let's remove them
    if any(opt.events.onset > EEG.xmax)  % in s
        idx = opt.events.onset > EEG.xmax;
        warning("Removing %g events from tsv table that have latency onset later than the iEEG data file length!", sum(idx))
        opt.events(idx,:) = [];
    end

else
    disp("You did not select a .tsv event file. Checking if your EEGLAB dataset already contains events...")
    if ~isempty(EEG.event)
        fprintf("EEGLAB dataset already contains %d events; using them as-is.\n", numel(EEG.event));
        opt.events = struct2table(EEG.event);
    else
        warning("No events detected in the EEGLAB dataset you loaded --> Analysis mode: CONTINUOUS.")
    end
end

% Electrode XYZ coordinates from TSV elec data
if ~isempty(opt.elec_tsv)
    EEG = get_elec_coor(EEG, elecs);
    opt.elec_labels = {EEG.chanlocs.labels};
end

% Abort if no channels left
if EEG.nbchan == 0
    error("No electrodes left in dataset.")
end

% Abort if there are no electrode locations
if ~isfield(EEG.chanlocs, 'X') || isempty([EEG.chanlocs.X])
    error("No electrode locations. Cannot proceed with analysis. You must load a valid .tsv file containing the electrodes' XYZ coordinates.")
end


% gui_load2 to select channels of interest
[opt, wasCanceled] = ieeglab_gui_load2(opt);
if wasCanceled || isempty(opt)
    % Restore the input dataset: get_elec_coor has already rewritten chanlocs
    % by this point, and storing that half-processed state as if the user had
    % confirmed was wrong.
    EEG = EEG_in;
    return
end


% Integrate event types and latencies from TSV event data into EEGLAB
% dataset
if isfield(opt,'events') && ~isempty(opt.events) && opt.events_from_tsv
    fprintf("Converting %g TSV events into EEGLAB format...\n", size(opt.events,1))
    col = opt.event_field;
    if ~ismember(col, opt.events.Properties.VariableNames)
        error('ieeglab_load:badEventField', ...
            'Event field "%s" is not a column of the events file. Available: %s', ...
            col, strjoin(opt.events.Properties.VariableNames, ', '));
    end
    ev_types = opt.events.(col);
    ev_lats  = opt.events.onset;

    % Replace any pre-existing events wholesale. Previously the loop only
    % overwrote slots 1..N, so a dataset that already carried more events kept
    % the tail, and stale duration/urevent fields survived on overwritten slots.
    EEG.event   = [];
    EEG.urevent = [];
    for iEv = 1:numel(ev_types)
        if iscell(ev_types(iEv))
            EEG.event(iEv).type = ev_types{iEv};
        else
            EEG.event(iEv).type = num2str(ev_types(iEv));
        end
        % EEGLAB latencies are in samples and 1-based: t=0 is sample 1.
        EEG.event(iEv).latency = ev_lats(iEv) * EEG.srate + 1;
        if ismember('duration', opt.events.Properties.VariableNames)
            d = opt.events.duration(iEv);
            if isnumeric(d) && isfinite(d), EEG.event(iEv).duration = d * EEG.srate; end
        end
    end
end
EEG = eeg_checkset(EEG, 'eventconsistency');
EEG = eeg_checkset(EEG);

% % Load events from ,tsv into EEGLAB
% if ~isempty(opt.events_tsv)
%     EEG = load_events(EEG,opt);
% end

% Keep only events of interest (if selected by user)
if isfield(opt,'event_to_keep') && ~isempty(opt.event_to_keep) ...
        && isfield(EEG,'event') && ~isempty(EEG.event)
    % ismember rather than ismissing: same result for cellstr but type-safe,
    % and the intent is membership, not missing-value detection.
    trials_to_rem = ~ismember(string({EEG.event.type}), string(opt.event_to_keep));
    nKeep = sum(~trials_to_rem);
    fprintf('Keeping %d/%d events of interest (%s).\n', ...
        nKeep, numel(trials_to_rem), strjoin(cellstr(string(opt.event_to_keep)), ', '));
    if nKeep == 0
        error('ieeglab_load:noEventsKept', ...
            ['None of the requested event types (%s) are present in the dataset.\n' ...
             'Types available: %s'], ...
             strjoin(cellstr(string(opt.event_to_keep)), ', '), ...
             strjoin(unique({EEG.event.type}), ', '));
    end
    EEG.event(trials_to_rem) = [];
    if isfield(opt,'events') && istable(opt.events) && height(opt.events) == numel(trials_to_rem)
        opt.events(trials_to_rem,:) = [];  % keep the table row-aligned
    end
    EEG = eeg_checkset(EEG, 'eventconsistency');
    EEG = eeg_checkset(EEG);
end

% Keep only electrodes of interest.
% Parenthesised: '&&' binds tighter than '||', so the unparenthesised version
% could take the chan_list branch while chan_idx was absent and then crash
% reading opt.chan_idx. Both fields are now required before entering.
haveIdx  = isfield(opt,'chan_idx')  && ~isempty(opt.chan_idx);
haveList = isfield(opt,'chan_list') && ~isempty(opt.chan_list);
subsetIdx  = haveIdx  && (numel(opt.chan_idx)  ~= EEG.nbchan || ~all(opt.chan_idx));
subsetList = haveList && (numel(opt.chan_list) ~= EEG.nbchan);
if (subsetIdx || subsetList) && haveList

    opt.chan_list = cellstr(opt.chan_list);  % ensure cellstr
    % Derive the removed labels from the kept list rather than from chan_idx,
    % which is not guaranteed to be a full-length logical mask.
    allLabels     = {EEG.chanlocs.labels};
    removed_elecs = setdiff(allLabels, opt.chan_list);
    if isempty(removed_elecs)
        fprintf('All %d electrodes selected for analysis.\n', EEG.nbchan);
    else
        fprintf('Removing %d/%d electrodes not selected for analysis: %s\n', ...
            numel(removed_elecs), EEG.nbchan, strjoin(removed_elecs, ', '));
        EEG = pop_select(EEG, 'channel', opt.chan_list);
    end

    % For CCEP data, drop events whose stimulation site references a removed
    % electrode. Guarded for continuous datasets, which have no events.
    if ~isempty(removed_elecs) && isfield(EEG,'event') && ~isempty(EEG.event) ...
            && any(contains({EEG.event.type}, '-'))
        trials_to_rem = contains({EEG.event.type}, removed_elecs);
        if any(trials_to_rem)
            fprintf('Removing %d/%d events whose stimulation site uses a removed electrode.\n', ...
                sum(trials_to_rem), numel(trials_to_rem));
            EEG.event(trials_to_rem) = [];
            if isfield(opt,'events') && istable(opt.events) && height(opt.events) == numel(trials_to_rem)
                opt.events(trials_to_rem,:) = [];
            end
            EEG = eeg_checkset(EEG, 'eventconsistency');
            EEG = eeg_checkset(EEG);
        end
    end

    % opt = rmfield(opt, 'chan_idx');
    % opt = rmfield(opt, 'chan_list');
end


% % Clear vars
% opt = rmfield(opt, 'elec_tsv');
% opt = rmfield(opt, 'elecs');
% opt = rmfield(opt, 'elec_labels');
% opt = rmfield(opt, 'events');
% opt = rmfield(opt, 'event_field');
% opt = rmfield(opt, 'event_values');
% opt = rmfield(opt, 'events_tsv');

EEG.ieeglab.opt = opt;

% Final check
EEG = eeg_checkset(EEG);

% Non-empty history string signals success to eeglab_new, which then stores the
% dataset and refreshes the main EEGLAB window.
com = 'EEG = ieeglab_load(EEG);';
