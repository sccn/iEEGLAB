%% iEEGLAB plugin - prep tutorial datasets

clear; close all; clc

% Long sEEG .mefd dataset (with pial surface)
% dataset source is subject 2 from: https://nemar.org/dataexplorer/detail?dataset_id=ds004696
filepath = '/Users/cedriccannard/Documents/dataset3';
filename = 'sub-02_ses-ieeg01_task-ccep_run-01_ieeg.mefd';
% EEG = pop_MEF3(fullfile(filepath, filename));

% [metadata, signal] = ieeglab_load_mefd(fileName);
metadata     = readMef3(fullfile(filepath, filename));
fs           = double(metadata.time_series_metadata.section_2.sampling_frequency);
num_samples  = double(metadata.time_series_metadata.section_2.number_of_samples);
startSample  = 369 * fs;     % import starting at 1st event
endSample    = 1509 * fs;    % 19 min of data

p = gcp('nocreate');
if isempty(p), parpool; end
tic
[metadata, signal] = ieeglab_load_mefd(fullfile(filepath, filename), [], [], ...
    'samples', [startSample endSample]);
toc

% Convert to EEGLAB format
EEG = eeg_emptyset;
EEG.filepath = filepath;
EEG.data = signal;
EEG.nbchan = size(EEG.data,1);
EEG.pnts = size(EEG.data,2);
EEG.srate = fs;
for iChan = 1:EEG.nbchan
    EEG.chanlocs(iChan).labels = metadata.time_series_channels(iChan).name;
end
EEG = eeg_checkset(EEG);

% Downsample
EEG = pop_resample(EEG, 128);


% --> Adjust events onset base from tsv file to match truncated EEG data
events = readtable('/Users/cedriccannard/Documents/seeg/sub-02_ses-ieeg01_task-ccep_run-01_events.tsv', 'FileType', 'text', 'Delimiter', '\t');
events.onset = events.onset - 369;
events.onset(1)

% remove events that are now out of bounds relative to truncated EEG
% dataset
idx = events.onset > EEG.xmax;
warning("Removing %g/%g events that are out of bounds.", sum(idx), size(events,1))
events(idx,:) = [];
writetable(events, '/Users/cedriccannard/Documents/seeg/sub-02_ses-ieeg01_task-ccep_run-01_events_ADJ.tsv','FileType', 'text', 'Delimiter', '\t')

% Remove unnecessary electrodes from removed events
unique_ev = unique(events.electrical_stimulation_site);
elecs = [];
for iEv = 1:length(unique_ev)
    elecs = [elecs extractBefore(unique_ev(iEv), '-')];
    elecs = [elecs extractAfter(unique_ev(iEv), '-')];
end
EEG = pop_select(EEG, 'channel', elecs);
 
% save EEGLAB dataset
pop_saveset(EEG, 'filename', sprintf('%s.set', extractBefore(filename, '.')),'filepath','/Users/cedriccannard/Documents/seeg/');


%% eCoG .vhdr data (with pial surface)
% dataset source: https://openneuro.org/datasets/ds005953/versions/1.0.0

filepath = '/Users/cedriccannard/Documents/ecog/sub-02/ses-01/ieeg/';
filename = 'sub-02_ses-01_task-visual_run-01_ieeg.vhdr';
EEG = pop_loadbv(filepath,filename);
EEG.filepath = filepath;
EEG.srate = round(EEG.srate);
EEG = eeg_checkset(EEG);

% Downsample
EEG = pop_resample(EEG, 128);

% save EEGLAB dataset
pop_saveset(EEG, 'filename', sprintf('%s.set', extractBefore(filename, '.')),'filepath',filepath);

