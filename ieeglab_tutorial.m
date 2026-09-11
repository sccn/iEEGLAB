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

% Heatmap of the 20% trimmed mean over trials, all channels
plot_ccep(EEG.data, EEG.times, {EEG.chanlocs.labels}, 'all', [], 0.20);

% All trials for one channel, with the trimmed mean on top.
% Pick the trials belonging to one stimulation site.
site   = 'RA1-RA2';                                  % change to a site in your data
sites  = ieeglab_epoch_sites(EEG);                   % canonical, order-independent
[cs, stimmed] = ieeglab_canonical_site({site}, {EEG.chanlocs.labels});
trials = sites == cs(1);
fprintf('%d trials at site %s\n', sum(trials), site);
if any(trials)
    % a contact NOT stimulated at this site (a stimulated one shows only the artifact)
    channel = find(~ismember(1:EEG.nbchan, stimmed{1}), 1);
    plot_ccep(EEG.data(:,:,trials), EEG.times, {EEG.chanlocs.labels}, 'single', channel, 0.20);
end

% For a conventional two-condition design (the eCoG dataset), compare stimuli:
% types = cellstr(ieeglab_epoch_sites(EEG));
% u     = unique(types);
% idx1  = strcmp(types, u{1});
% idx2  = strcmp(types, u{2});
% plot_ccep(EEG.data(:,:,idx1), EEG.times, {EEG.chanlocs.labels}, 'single', 2, 0.20);


%% Step 8: CCEP analysis - N1, CRP and the connectivity matrix
% Menu: iEEGLAB > CCEP analysis (N1, CRP, connectivity)
%
% Each stage can be switched on or off in the dialog. N1 gives amplitude and
% latency of the early response; CRP gives its canonical shape and duration;
% the connectivity matrix summarises which stimulation site drives which contact.

EEG = ieeglab_stats_subject(EEG);

% The connectivity matrix: rows are stimulation sites, columns contacts.
% Grey = not measured (the stimulated contacts), white = no response.
ieeglab_plot_ccep_matrix(EEG.ieeglab.ccep_matrix);                    % responses
ieeglab_plot_ccep_matrix(EEG.ieeglab.ccep_matrix, 'latency_ms');      % N1 latency

% The same results on the brain: how many sites evoke a response at each
% contact, and the N1 amplitude for one stimulation site.
% Menu: iEEGLAB > Plot electrode values on brain (also under EEGLAB's Plot menu)
ieeglab_topoplot(EEG, 'in_degree');
ieeglab_topoplot(EEG, 'n1_amplitude', 'site', EEG.ieeglab.ccep_matrix.sites{1});

% Export everything as TSV / JSON / MAT.  Menu: iEEGLAB > Export results
ieeglab_export(EEG, fullfile(filepath, 'derivatives', 'ieeglab'));


%% Step 9: the same pipeline with no dialogs at all
% Every step takes an options struct, so this block runs unattended - which is
% what makes batch processing and automated testing possible. ieeglab_load
% finds the BIDS sidecar files next to the dataset on its own.

filepath = fullfile(plugin_path, 'tutorial', 'dataset_seeg');
EEG = pop_loadset('filename','sub-02_ses-ieeg01_task-ccep_run-01_ieeg.set', 'filepath',filepath);

% Coordinates, events (events marked status=bad are dropped) and clinician
% channel labels, from the _electrodes.tsv / _events.tsv / _channels.tsv sidecars
EEG = ieeglab_load(EEG, struct('event_field', 'electrical_stimulation_site'));

% Preprocess. Blanking needs data at its native rate; this extract is 128 Hz,
% so it is left off here (the function warns about exactly this).
EEG = ieeglab_preprocess(EEG, struct( ...
    'remove_bad_channels', true, ...
    'apply_highpass', true,  'highpass', 0.5, ...
    'apply_notch',    true,  'notch', [60 120 180], ...
    'apply_lowpass',  false, ...
    'apply_epoch',    true,  'epoch_window', [-500 1000], ...
    'apply_car',      true,  'car_method', 'carla', ...
    'apply_baseline', true,  'baseline_period', [-500 -50], ...
    'plot', false));

% N1, CRP, connectivity matrix and export, in one call
outdir = fullfile(tempdir, 'ieeglab_tutorial_results');
EEG = ieeglab_stats_subject(EEG, struct('run_n1', true, 'run_crp', true, 'run_matrix', true, ...
    'crp_window', [15 400], 'n1_window', [15 100], 'export_dir', outdir, 'plot', false));

M = EEG.ieeglab.ccep_matrix;
fprintf('\n%d of %d tested site-contact pairs respond (density %.2f). Results in %s\n', ...
    M.n_significant, M.n_tested, M.density, outdir);
