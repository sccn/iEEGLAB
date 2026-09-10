%% iEEGLAB tutorial
%
% Two shipped datasets are used:
%   dataset_seeg - stereo-EEG, single-pulse stimulation (CCEP), 16 channels
%                  subject 2 of https://nemar.org/dataexplorer/detail?dataset_id=ds004696
%   dataset_ecog - electrocorticography, visual task, 96 channels
%
% Both are reduced .set files that carry NO events and NO electrode
% coordinates: those come from the BIDS .tsv sidecars, which is what step 3
% demonstrates.
%
% Sections 1-8 are the interactive walkthrough (matching the wiki).
% Section 9 is the same pipeline run entirely from the command line - use that
% one for batch processing and as the regression test.
%
% Cedric Cannard, iEEGLAB, 2025-2026

%% Step 1: launch EEGLAB and locate the plugin

clear; close all; clc
eeglab; close

plugin_path = fileparts(which('eegplugin_ieeglab'));
cd(plugin_path)

% Confirm the dependencies are in place before going further
ieeglab_check_install


%% Step 2: load one of the two datasets
% Run ONE of these two blocks, not both.

% --- (a) sEEG / CCEP dataset ---
filepath = fullfile(plugin_path, 'tutorial', 'dataset_seeg');
filename = 'sub-02_ses-ieeg01_task-ccep_run-01_ieeg.set';
EEG = pop_loadset('filename', filename, 'filepath', filepath);

% --- (b) eCoG / visual task dataset ---
% filepath = fullfile(plugin_path, 'tutorial', 'dataset_ecog');
% filename = 'sub-02_ses-01_task-visual_run-01_ieeg.set';
% EEG = pop_loadset('filename', filename, 'filepath', filepath);


%% Step 3: load electrode coordinates and events from the .tsv sidecars
% Menu: iEEGLAB > Load electrode coordinates and events
%
% The first dialog picks electrodes.tsv and events.tsv; the second selects the
% channels and event types to keep. For the sEEG set, choose
% 'electrical_stimulation_site' as the column holding the event name - that is
% what makes the stimulated contact pair known to the rest of the pipeline.

EEG = ieeglab_load(EEG);


%% Step 4: visualise raw data and electrodes
% Menu: iEEGLAB > Visualize electrodes

% Raw time series
pop_eegplot(EEG, 1, 1, 1);

% 3D glass brain. Needs vistasoft on the path:
%   addpath(genpath('<your path>/vistasoft'))
% You may select several meshes at once in the file dialog.
EEG = ieeglab_vis_elec(EEG);

% Or name the surfaces directly. The sEEG dataset ships pial, white and
% inflated surfaces for both hemispheres.
EEG = ieeglab_vis_elec(EEG, struct('surf_files', ...
    {{'pial.L.surf.gii', 'pial.R.surf.gii'}}));


%% Step 5: preprocess
% Menu: iEEGLAB > Preprocess iEEG data
%
% Filtering and resampling run on the continuous data, then epoching,
% re-referencing (CARLA) and baseline correction.

EEG = ieeglab_preprocess(EEG);


%% Step 6: standard EEGLAB plots

% Epoched time series
pop_eegplot(EEG, 1, 1, 1);

% ERP image for a few channels
for iChan = 1:3:EEG.nbchan
    figure
    pop_erpimage(EEG, 1, iChan, [], EEG.chanlocs(iChan).labels, 10, 1, {}, [], '', ...
        'yerplabel','\muV', 'erp','on', 'cbar','on');
    pause(2)
    close(gcf)
end

% Butterfly plot over the actual epoch window (not a hardcoded one)
figure;
plottopo(mean(EEG.data,3), 'frames', EEG.pnts, ...
    'limits', [EEG.times(1) EEG.times(end) 0 0], ...
    'chans', 1:EEG.nbchan, 'ydir', 1);

% Power spectrum
figure; pop_spectopo(EEG, 1, [EEG.times(1) EEG.times(end)], 'EEG', ...
    'freqrange', [1 min(100, EEG.srate/2 - 1)], 'electrodes', 'off');


%% Step 7: CCEP plots
% Only meaningful for the sEEG / CCEP dataset.

% Heatmap of the trimmed mean across all channels
plot_ccep(mean(EEG.data,3), EEG.times, {EEG.chanlocs.labels}, 'all', [], 0.20);

% All trials for one channel, with the trimmed mean on top.
% Pick the trials belonging to one stimulation site.
site   = 'RA1-RA2';                                  % change to a site in your data
trials = strcmpi(local_epoch_sites(EEG), site);
fprintf('%d trials at site %s\n', sum(trials), site);
if any(trials)
    channel = 2;
    plot_ccep(EEG.data(:,:,trials), EEG.times, {EEG.chanlocs.labels}, 'single', channel, 0.20);
end

% For a conventional two-condition design (the eCoG dataset), compare stimuli:
% types = local_epoch_sites(EEG);
% u     = unique(types);
% idx1  = strcmp(types, u{1});
% idx2  = strcmp(types, u{2});
% plot_ccep(EEG.data(:,:,idx1), EEG.times, {EEG.chanlocs.labels}, 'single', 2, 0.20);


%% Step 8: within-subject statistics (CRP)
% Menu: iEEGLAB > Within-subject statistics (CRP)
%
% Fits the canonical response shape for every stimulation-site x channel pair
% and reports response duration (tau_R), explained variance and significance.

EEG = ieeglab_stats_subject(EEG);

T = EEG.ieeglab.stats.table;
disp(T(T.significant, :))

% Most responsive pair
if any(T.significant)
    [~, k] = max(T.explained_var);
    fprintf('Strongest response: %s -> %s, tau_R = %.0f ms, R2 = %.2f\n', ...
        T.site{k}, T.channel{k}, T.tR_ms(k), T.explained_var(k));
end

% Save
% pop_saveset(EEG, 'filename', 'sub-02_processed.set', 'filepath', filepath);


%% Step 9: the same pipeline with no GUI at all
% This is the scriptable form. Every step takes an options struct, so this
% block runs unattended - which is what makes batch processing and automated
% testing possible.

filepath = fullfile(plugin_path, 'tutorial', 'dataset_seeg');
EEG = pop_loadset('filename','sub-02_ses-ieeg01_task-ccep_run-01_ieeg.set', 'filepath',filepath);

% Electrode coordinates from the BIDS sidecar
elecs = readtable(fullfile(filepath,'sub-02_ses-ieeg01_electrodes.tsv'), ...
    'FileType','text', 'Delimiter','\t');
EEG = get_elec_coor(EEG, elecs);

% Events from the BIDS sidecar, using the stimulation site as the event type
events = readtable(fullfile(filepath,'sub-02_ses-ieeg01_task-ccep_run-01_events.tsv'), ...
    'FileType','text', 'Delimiter','\t');
EEG.event = [];
for iEv = 1:height(events)
    EEG.event(iEv).type    = events.electrical_stimulation_site{iEv};
    EEG.event(iEv).latency = events.onset(iEv) * EEG.srate + 1;   % samples, 1-based
end
EEG = eeg_checkset(EEG, 'eventconsistency');
EEG.ieeglab.opt.events = events;

% Preprocess
EEG = ieeglab_preprocess(EEG, struct( ...
    'apply_highpass', true,  'highpass', 0.5, ...
    'apply_notch',    true,  'notch', [60 120 180], ...
    'apply_lowpass',  false, ...
    'apply_epoch',    true,  'epoch_window', [-500 1000], ...
    'apply_car',      true,  'car_method', 'carla', ...
    'apply_baseline', true,  'baseline_period', [-500 -50], ...
    'plot', false));

% Statistics
EEG = ieeglab_stats_subject(EEG, struct('crp_window',[15 400], 'plot',false));

T = EEG.ieeglab.stats.table;
fprintf('\n%d of %d site-channel pairs show a significant response.\n', ...
    sum(T.significant), height(T));


%% local helper used by step 7
function sites = local_epoch_sites(EEG)
% The event type each epoch is locked to, as a cellstr.
sites = repmat({''}, 1, EEG.trials);
if ~isfield(EEG,'epoch') || isempty(EEG.epoch), return; end
for i = 1:numel(EEG.epoch)
    t = EEG.epoch(i).eventtype;
    if iscell(t)
        l = EEG.epoch(i).eventlatency;
        if iscell(l) && numel(l) == numel(t)
            [~, k] = min(cellfun(@(x) abs(double(x(1))), l));
        else
            k = 1;
        end
        t = t{k};
    end
    if isnumeric(t), t = num2str(t); end
    sites{i} = char(t);
end
end
