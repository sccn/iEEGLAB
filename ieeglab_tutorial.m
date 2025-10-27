%% step 1: launch EEGLAB

clear; close all; clc
eeglab; close

plugin_path = fileparts(which('eegplugin_ieeglab'));
cd(plugin_path)

%% Reduced sEEG .set dataset (with pial surface)

filepath = fullfile(plugin_path, 'tutorial', 'seeg_set');
filename = 'sub-02_ses-ieeg01_task-ccep_run-01_ieeg.set';
EEG = pop_loadset('filename', filename, 'filepath', filepath);

%% Reduced eCoG .set dataset (with pial surface)

filepath = fullfile(plugin_path, 'tutorial', 'ecog');
filename = 'sub-02_ses-01_task-visual_run-01_ieeg.set';
EEG = pop_loadset('filename', filename, 'filepath', filepath);


%% Load electrode coordinates and events of interest from .tsv files

EEG = ieeglab_load(EEG);


%% Visualize 

% Raw data
pop_eegplot(EEG,1,1,1);


% Electrodes in 3D glass brain 
addpath(genpath('/Users/cedriccannard/Documents/MATLAB/vistasoft'))

EEG = ieeglab_vis_elec(EEG);

% Or you can define the surface files by command line directly
EEG.ieeglab.opt.surf_files = {'pial_desc-qsiprep.L.surf.gii' 'pial_desc-qsiprep.R.surf.gii'};
ieeglab_vis_elec(EEG);

%% Preprocess channels

EEG = ieeglab_preprocess(EEG); 


%% Some EEGLAB plots

% Plot Epoched time series
pop_eegplot(EEG,1,1,1);

% ERP image
for iChan = 1:25:EEG.nbchan
    figure
    pop_erpimage(EEG,1, iChan,[],EEG.chanlocs(iChan).labels,10,1,{},[],'' ,'yerplabel','\muV','erp','on','cbar','on');
    pause(1)
    close(gcf)
end

figure; 
pop_timtopo(EEG, [EEG.times(1) EEG.times(end)], [], ...
    'CCEP - all electrodes','verbose','off');

% --> NEED TO GET THE NEW 3D BRAIN PLOT INSTEAD
% figure; pop_timtopo(EEG, [-500  990], NaN, 'ERP data and scalp maps');

% figure; pop_plottopo(EEG, 1:EEG.nbchan , '', 0, 'ydir',1);
figure; plottopo( trimmean(EEG.data,20,3), 'frames', EEG.pnts, 'limits', [-500 990 0 0], 'chans', 1:EEG.nbchan, 'ydir', 1);
% figure; plottopo( EEG.data, 'frames', EEG.pnts, 'limits', [-500 990 0 0], 'chans', 1:EEG.nbchan, 'ydir', 1);

% PSD
figure; pop_spectopo(EEG, 1, [], 'EEG' , 'freq', [], 'freqrange',[0.5 50],'electrodes','off');

% Heatmap of the average across trials (all channels x time)
avg_pair = squeeze(trimmean(EEG.data, 10, 3));             % [chan x time]
plot_ccep(avg_pair, EEG.times, {EEG.event.type}, 'all', [], 0.20);

% single-channel overlay of all trials + mean (CCEP experiment)
channel = 2;   
% trials = contains({EEG.event.type}, EEG.chanlocs(channel).labels);  % for CCEP trials
sum(trials)
plot_ccep(EEG.data(:,:,trials), EEG.times, {EEG.chanlocs.labels}, 'single', channel, []);

% single-channel overlay of all trials + mean (traditional experiment
% comparing two stimuli)
channel = 2;   
idx1 = strcmp({EEG.event.type}, {EEG.event(1).type});  % for non-CCEP data (this corresponds to stimulus type)
idx2 = strcmp({EEG.event.type}, {EEG.event(2).type});  % for non-CCEP data (this corresponds to stimulus type)
plot_ccep(EEG.data(:,:,idx1), EEG.times, {EEG.chanlocs.labels}, 'single', channel, []);
% plot_ccep(EEG.data(:,:,idx2), EEG.times, {EEG.chanlocs.labels}, 'single', channel, []);

% figure; pop_newtimef( EEG, 1, 1, [-203  789], [3         0.8] , 'topovec', 1, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', 'RA1', 'baseline',[0], 'alpha',0.001, 'mcorrect', 'fdr', 'padratio', 1, 'winsize', 26);

%%

% look through each pair?
% iPair = 1
pairs = unique({EEG.event.type});
avg_pair = {};
for iPair = 1:length(pairs)
    pair = nPairs{iPair};
    fprintf('---------- PAIR %g/%g: %s ----------\n', iPair, length(pairs), pair)
    trials = find(strcmpi({EEG.event.type}, pair));

    stimElec = extractBefore(pair, '-');
    respElec = extractAfter(pair, '-');
    stimElecIdx = strcmpi({EEG.chanlocs.labels}, stimElec);
    respElecIdx = strcmpi({EEG.chanlocs.labels}, respElec);

    if any(respElecIdx)
        % % plot all single trials for this pair
        figure; hold on
        for iTrial = 1:length(trials)
            % plot(EEG.times, squeeze(EEG.data(stimElecIdx, :, trials(iTrial))))
            plot(EEG.times, squeeze(EEG.data(respElecIdx, :, trials(iTrial))))
            % legend({'stim. elec' 'meas. elec'})
        end
    end
    title(sprintf("Pair %g/%g: %s", iPair, length(pairs), pair))
    pause(1); close(gcf)
end


%% CRP per eligible channel 

pairs = unique({EEG.event.type});

nPairs = length(pairs);
nEv = length(EEG.event);

% % Initialize
% crp_out = repmat(struct('data',[],'tt',[],'crp_parms',[],'crp_projs',[]), EEG.data, nPairs);
% average_ccep = NaN(EEG.nbchan, max(nPairs,1), length(EEG.times));

stim_pair_nr    = NaN(nEv,1);   
stim_pair_name  = cell(nEv,1); 

condition_type_counter = 0;
% get all stim + electrodes
stimEl1         = extractBefore({EEG.event.type},'-');
% get all stim - electrodes
stimEl2         = extractAfter({EEG.event.type},'-');

for iEv = 1:nEv

    % which electrodes are stimulated
    el1 = stimEl1{iEv};
    el2 = stimEl2{iEv};

    % is this trial a stimulation trial: does el1 have content?
    if ~isempty(el1) && ~isempty(el2)
        % if this trial type does not exist yet & is a stimulation trial
        if sum(strcmp({EEG.event.type},[el1 '-' el2]))==0 && ...  % does el1-el2 already exist?
                sum(strcmp({EEG.event.type},[el2 '-' el1]))==0    % group el2-el1 with el1-el2
            condition_type_counter = condition_type_counter+1;

            % find all trials with el1 & el2 | el2 & el1
            theseTrials = strcmp(stimEl1,el1) & strcmp(stimEl2,el2) | ...
                strcmp(stimEl2,el1) & strcmp(stimEl1,el2);
            trial_nrs = find(theseTrials==1);                   % number of trials of this type
            for ll = 1:sum(theseTrials)
                stim_pair_name{trial_nrs(ll),1} = [el1 '-' el2];
            end
            stim_pair_nr(theseTrials) = condition_type_counter;
        end
    end
end
clear el1 el2 trial_nrs theseTrials epoch_type_counter


% load data for each condition and put average in matrix
% average_ccep: matrix with average cceps (measured electrodes X stim pairs X time)
% average_ccep_names: names of stimulated channels

% set epoch parameters
% epoch_length            = 5; % in seconds, -1:3
% epoch_prestim_length    = 2; % in seconds
% tt = (1:epoch_length*fs)/fs - epoch_prestim_length;
tt = EEG.times ./ 1000;
epoch_length = EEG.xmax;
epoch_prestim_length = (EEG.times(EEG.times==0) - EEG.times(1)) ./ 1000;

% initialize output
average_ccep            = NaN(EEG.nbchan, max(stim_pair_nr), epoch_length * EEG.srate);
average_ccep_names      = cell(max(stim_pair_nr),1);
crp_out                 = [];

for iPair = 1:max(stim_pair_nr) % condition number
    disp(['loading data for condition ' int2str(iPair) ' out of ' int2str(max(stim_pair_nr))])
    
    % epochs of this condition
    these_epochs        = find(stim_pair_nr==iPair & epochs_include==1);
    
    if ~isempty(these_epochs)
    
        % save name of the current epoch
        average_ccep_names{iPair} = stim_pair_name{these_epochs(1)};

        % for this condition number (iEv), load each epoch    
        all_start           = round((events_table.onset(these_epochs)-epoch_prestim_length)*srate);
        all_end             = round((events_table.onset(these_epochs)+epoch_length-epoch_prestim_length)*srate);
        epoch_ranges        = [all_start all_end];
        [~,signaldata]      = readMef3(fileName,[], [], 'samples', epoch_ranges);  % read all channels, samples 0-1000 
        
        % exclude stimulated channels from the good channels
        stimEl1_nr          = find(ismember(channel_names,extractBefore(average_ccep_names{iPair},'-')));
        stimEl2_nr          = find(ismember(channel_names,extractAfter(average_ccep_names{iPair},'-')));

        good_channels_car   = setdiff(good_channels,[stimEl1_nr stimEl2_nr]);

        % adjusted Common Average Reference (CAR)
        if use_CAR == 1
            perc_channels   = 0.2;
            car_timeint     = [0.015 0.500];
            [signaldata]    = ccep_CAR64blocks_percent(signaldata,tt,good_channels_car,perc_channels,car_timeint);
        end

        signaldata          = permute(signaldata,[1 3 2]);

        % baseline subtract
        samples_base        = find(tt>baseline_t(1) & tt<baseline_t(2));
        data_epoch          = ieeg_baselinesubtract(signaldata,samples_base,'median');
        
        % run CRP for channels in the limbic network
        for ii = 1:size(channel_areas,1)
            if channel_areas(ii)>1  && ismember(ii,good_channels_car) % && size(data_epoch,2)>3% only limbic 
                crp_out(ii,iPair).data = squeeze(data_epoch(ii,:,:));
                crp_out(ii,iPair).tt   = tt;
                V                   = squeeze(data_epoch(ii,:,tt>t_win_cod(1) & tt<t_win_cod(2)));
                t_win               = tt(tt>t_win_cod(1) & tt<t_win_cod(2));
                [crp_parms, crp_projs] = CRP_method(V',t_win);
                crp_out(ii,iPair).crp_parms = crp_parms;
                crp_out(ii,iPair).crp_projs = crp_projs;
            else
                crp_out(ii,iPair).data = [];
                crp_out(ii,iPair).tt   = [];
                crp_out(ii,iPair).crp_parms = [];
                crp_out(ii,iPair).crp_projs = [];
            end
        end

        % put average in matrix
        average_ccep(:,iPair,:)        = squeeze(nanmean(data_epoch,2));

    else
        % save name of the current epoch
        average_ccep_names{iPair}      = stim_pair_name{find(stim_pair_nr==iPair,1)};
    end
    if return_single_trials         == 1 && max(stim_pair_nr) == 1
        single_trials = data_epoch;
    elseif return_single_trials     == 1 && max(stim_pair_nr) > 1
        disp('can not return single trials, only works for 1 stim pair')
        single_trials = [];
    end
    clear these_epochs_data data_epoch ll_start ll_end
end

% make stimulated electrodes a NaN
for iPair = 1:size(average_ccep,2) % epochs
    % stimulated electrode names
    el1     = extractBefore(average_ccep_names{iPair},'-');
    el2     = extractAfter(average_ccep_names{iPair},'-');
    
    % stimulated electrode index in the data
    el1_nr  = ismember(channel_names,el1);
    el2_nr  = ismember(channel_names,el2);
    
    % set to NaN
    average_ccep(el1_nr==1,iPair,:) = NaN;
    average_ccep(el2_nr==1,iPair,:) = NaN;
    
    clear el1 el2 el1_nr el2_nr% housekeeping
end
% For each stim pair, record the anatomical areas of the two electrodes.
average_ccep_areas = NaN(nPairs, 2);
for p = 1:nPairs
    if p > numel(pair_to_rows) || isempty(pair_to_rows{p}), continue; end
    r = pair_to_rows{p};   % [rowA, rowB]
    average_ccep_areas(p, :) = [master_table.area(r(1)), master_table.area(r(2))];
end

% % Export (uncomment to save)
outDir   = fullfile(data_path,'derivatives','stats',['sub-' bids_sub]);
if ~exist(outDir,'dir'), mkdir(outDir); end
statsFile = fullfile(outDir, ...
  ['sub-' bids_sub '_ses-' bids_ses '_task-' bids_task '_run-' bids_run '_crp.mat']);
save('average_ccep','average_ccep_names','average_ccep_areas', ...
             'timevec','fs','crp_out','channel_names','channel_areas','pair_to_rows','pair_to_names');
% disp(['Saved: ' statsFile])

disp("Done computing statistics on all pairs")




