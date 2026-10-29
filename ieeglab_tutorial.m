%% iEEGLAB tutorial: CCEP analysis of sEEG and ECoG data
%
% Runs the whole analysis on the two tutorial datasets shipped with the plugin,
% from the command line. Each step names the menu item that does the same in
% the EEGLAB window, so the tutorial can be followed either way.
%
%   tutorial/seeg  sEEG, OpenNeuro ds004696 sub-02 (Ojeda Valencia et al., 2023):
%                  17 contacts, 3 stimulation sites, 34 pulses, 2048 Hz, with the
%                  pial surfaces
%   tutorial/ecog  ECoG, OpenNeuro ds004080 sub-ccepAgeUMCU02 (van Blooijs et al.,
%                  2023): 20 contacts, 3 sites, 30 pulses, 2048 Hz, coordinates in
%                  fsaverage
%
% Both are small extracts of the original recordings at their native sampling
% rate, one folder each: the .set file next to its BIDS sidecar files
% (tutorial/make_tutorial_datasets.m). On sEEG the CCEP measure is
% the Canonical Response Parameterization (CRP); N1 detection is used on ECoG.
%
% Set show_figures = false to run without figures.
%
% Cedric Cannard, iEEGLAB, 2025-2026

show_figures = true;

%% 1. Start EEGLAB
% Needs EEGLAB with the firfilt and EEG-BIDS plugins (File > Manage EEGLAB
% extensions). iEEGLAB > Check installation lists anything missing.
eeglab
plugin_path = fileparts(which('eegplugin_ieeglab'));
addpath(plugin_path, fullfile(plugin_path, 'functions'));   % if iEEGLAB is not in eeglab/plugins
ieeglab_check_install;
outroot = fullfile(tempdir, 'ieeglab_tutorial_results');

%% ===== Part A: sEEG =====================================================

%% A2. Open the dataset
% Menu: File > Load existing dataset, and select the .set file in tutorial/seeg.
% (For a full recording downloaded from OpenNeuro, use iEEGLAB > Import raw
% recording instead: see iEEGLAB_GUI_walkthrough.html.)
tut_seeg = fullfile(plugin_path, 'tutorial', 'seeg');
EEG = pop_loadset('filename', 'sub-02_ses-ieeg01_task-ccep_run-01_ieeg.set', 'filepath', tut_seeg);

%% A3. Electrode coordinates, events and clinician channel labels
% Menu: iEEGLAB > Load electrode coordinates and events.
% Reads the electrodes, events and channels files next to the dataset. The
% event type is the stimulation site (e.g. RA3-RA4). Contacts marked bad in
% channels.tsv (here RA15, RB14, RB15) are marked, and removed in step A5.
EEG = ieeglab_load(EEG, struct('event_field', 'electrical_stimulation_site'));

%% A4. Electrodes on the brain
% Menu: iEEGLAB > Choose brain surfaces (Use these: the participant's pial
% surfaces next to the dataset), then iEEGLAB > Visualize electrodes. Chosen
% once, the surfaces are used by every brain figure without asking.
surf_seeg = fullfile(tut_seeg, {'pial.L.surf.gii', 'pial.R.surf.gii'});
EEG = pop_ieeglab_surfaces(EEG, 'surf_files', surf_seeg);
if show_figures
    EEG = ieeglab_vis_elec(EEG);
end

%% A5. Preprocess
% Menu: iEEGLAB > Preprocess iEEG data.
% Blank the stimulation artifact, remove the line noise (60 Hz, United States, and
% its harmonics; sines fitted away from the pulses, no filter), cut epochs around each pulse,
% drop sites with fewer than 5 pulses, outlier trials and clinician-bad contacts,
% and remove the baseline.
EEG = ieeglab_preprocess(EEG, struct( ...
    'apply_blank',         true, ...
    'apply_notch',         true,  'notch', [60 120 180], ...
    'apply_epoch',         true,  'epoch_window', [-1000 1000], ...
    'reject_trials',       true, ...
    'remove_rare_cond',    true,  'min_trials', 5, ...   % as the window's defaults
    'remove_no_coords',    true, ...
    'remove_bad_channels', true, ...
    'apply_baseline',      true,  'baseline_period', [-500 -50], ...
    'plot', show_figures));

%% A6. Re-reference
% Menu: iEEGLAB > iEEG re-referencing.
% CARLA chooses, for each stimulation site, the channels whose common average
% does not contain the evoked response (Huang et al., 2024). As its authors run
% it on real data, the stimulated pair, 2 contacts on each side of it and the
% bad contacts are left out, with the sensitive cutoff. A plain common average and ICA
% re-referencing (Michelmann et al., 2018: components spread uniformly over the
% contacts are removed) are kept for comparison.
EEG_car = pop_ieeglab_reref(EEG, 'method', 'car');
EEG_ica = pop_ieeglab_reref(EEG, 'method', 'ica');
EEG     = pop_ieeglab_reref(EEG, 'method', 'carla');

%% A7. CCEP analysis: CRP and the connectivity matrix
% Menu: iEEGLAB > CCEP analysis (N1, CRP, connectivity).
% CRP (Miller et al., 2023) gives, for every site x contact, the response
% duration (tau_R), its shape and the variance it explains; iEEGLAB tests each
% pair with a permutation test (the published t-test at tau_R is an option).
EEG = ieeglab_stats_subject(EEG, struct('run_n1', false, 'run_crp', true, ...
    'run_matrix', true, 'matrix_source', 'crp', 'plot', show_figures));
M = EEG.ieeglab.ccep_matrix;
fprintf('sEEG: %d of %d tested site-contact pairs respond (CRP).\n', M.n_significant, M.n_tested);

%% A8. Results on the brain
% Menu: iEEGLAB > Plot connectivity matrix, and > Plot electrode values on brain
% (on the surfaces chosen in A4).
if show_figures
    ieeglab_plot_ccep_matrix(M);
    ieeglab_topoplot(EEG, 'in_degree');
    ieeglab_topoplot(EEG, 'crp_explained_var', 'site', 'RA3-RA4');
end

%% A9. Export
% Menu: iEEGLAB > Export results. TSV tables, a JSON record of every option
% used, and a MAT file.
ieeglab_export(EEG, fullfile(outroot, 'seeg'));

%% ===== Part B: ECoG =====================================================

%% B2-B3. Open the dataset and load the sidecars
tut_ecog = fullfile(plugin_path, 'tutorial', 'ecog');
EEG = pop_loadset('filename', 'sub-ccepAgeUMCU02_ses-1_task-SPESclin_run-041456_ieeg.set', 'filepath', tut_ecog);
EEG = ieeglab_load(EEG, struct('event_field', 'electrical_stimulation_site'));

%% B4. Electrodes
% Coordinates are in fsaverage space and no individual MRI is shared. Without a
% choice in Choose brain surfaces, every brain figure is drawn on the brain
% that matches the coordinates (ieeglab_brain_surfaces): fsaverage's pial
% surfaces (FreeSurfer or MNE's fsaverage), else the template brain.
if show_figures
    EEG = ieeglab_vis_elec(EEG);
end

%% B5. Preprocess (line noise at 50 Hz: recorded in the Netherlands)
EEG = ieeglab_preprocess(EEG, struct( ...
    'apply_blank',         true, ...
    'apply_notch',         true,  'notch', 50, ...
    'apply_epoch',         true,  'epoch_window', [-1000 1000], ...
    'reject_trials',       true, ...
    'remove_rare_cond',    true,  'min_trials', 5, ...   % as the window's defaults
    'remove_no_coords',    true, ...
    'remove_bad_channels', true, ...
    'apply_baseline',      true,  'baseline_period', [-500 -50], ...
    'plot', show_figures));

%% B6. Re-reference: common average of the good contacts
EEG = pop_ieeglab_reref(EEG, 'method', 'car');

%% B7. CCEP analysis: N1 and the connectivity matrix
% N1 (van Blooijs et al., 2018): amplitude and latency of the largest negative
% peak 10-100 ms after the pulse; iEEGLAB tests it with a permutation test (the
% 3.4 SD rule of erdetect is an option).
EEG = ieeglab_stats_subject(EEG, struct('run_n1', true, 'run_crp', false, ...
    'run_matrix', true, 'matrix_source', 'n1', 'plot', show_figures));
M = EEG.ieeglab.ccep_matrix;
fprintf('ECoG: %d of %d tested site-contact pairs respond (N1).\n', M.n_significant, M.n_tested);

%% B8-B9. Plots and export
if show_figures
    ieeglab_plot_ccep_matrix(M, 'latency_ms');
    ieeglab_topoplot(EEG, 'n1_amplitude', 'site', 'FB54-FB62');
end
ieeglab_export(EEG, fullfile(outroot, 'ecog'));
fprintf('Results in %s\n', outroot);
