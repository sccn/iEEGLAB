function [signaldata, out] = apply_ieeg_car(signaldata, ttt, good_channels, perc_channels, car_timeint, preserve_rank)
% APPLY_IEEG_CAR  Apply Common Average Reference (CAR) to iEEG data with block-wise variance selection.
%
%   Adapted from HAPwave / Mayo Clinic Multimodal Neuroimaging Lab.
%   Adds safety checks, dimension handling, optional rank check, and improved reporting.
%
%   [signaldata, out] = apply_ieeg_car(signaldata, ttt, good_channels, perc_channels, car_timeint, preserve_rank)
%
% INPUTS:
%   signaldata     - Channels × Time × Epochs numeric array
%   ttt            - Time vector in seconds (length must match Time dimension or Epoch dimension)
%   good_channels  - Indices of channels to include in CAR pool (exclude noisy or stim channels)
%   perc_channels  - Proportion of channels with lowest variance to include in CAR (0.1 = lowest 10%)
%   car_timeint    - Time window [start, end] in seconds to compute variance (default: [0.015 0.500])
%   preserve_rank  - Logical flag, if true, appends surrogate channel to reduce risk of rank loss.
%                    NOTE: For iEEG block-wise CAR with channel exclusions, rank preservation is not guaranteed.
%
% OUTPUTS:
%   signaldata     - Channels × Time × Epochs after CAR applied
%   out            - Struct array with CAR details for each block of 64 channels
%
% WARNINGS:
%   - Prints skipped block indices once at the end.
%   - Warns if rank changes after CAR.
%   - Warns if data appears to already have CAR applied.
%
% DH & HH Multimodal Neuroimaging Lab, Mayo Clinic, 2020
% Adapted & expanded by Cedric Cannard, 2025.
% Copyright (c) 2025 iEEGLAB Plugin Project.

if nargin < 4 || isempty(perc_channels)
    perc_channels = 0.30;
end
if nargin < 5 || isempty(car_timeint)
    car_timeint = [0.015 0.500];
end
if nargin < 6
    preserve_rank = false;
end

out = [];

% --- Dimension check: channels × time × epochs expected ---
if length(ttt) ~= size(signaldata, 2)
    if length(ttt) == size(signaldata, 3)
        warning('Permuting signaldata from channels×epochs×time to channels×time×epochs');
        signaldata = permute(signaldata, [1 3 2]);
    else
        error('Time vector length (%d) does not match any signaldata dimension (%d or %d).', ...
              length(ttt), size(signaldata, 2), size(signaldata, 3));
    end
end

% --- Detect possible pre-CAR data (high mean channel correlation) ---
mean_corr = mean(triu(corr(squeeze(mean(signaldata,3))'),1),'all','omitnan');
if mean_corr < 0.01
    warning('Data appears to have very low inter-channel correlation — CAR may have been applied already.');
end

% --- Rank before CAR ---
% r_before = rank(double(reshape(signaldata, size(signaldata,1), [])));
r_before = sum(eig(cov(double(signaldata(:,:)'))) > 1E-7);


% --- Create channel blocks BEFORE surrogate channel ---
nChans_orig = size(signaldata, 1);
set_nrs = 1:ceil(nChans_orig/64);
for ss = set_nrs
    set_inds = (ss*64-63) : ss*64; 
    set_inds = set_inds(set_inds <= nChans_orig);
    out(ss).channels_set = set_inds;
end

% --- Add surrogate channel if preserving rank ---
if preserve_rank
    signaldata(end+1, :, :) = 0; % append zero-filled surrogate channel
    surrogate_idx = size(signaldata, 1);
    % Add surrogate channel to last block
    out(end).channels_set(end+1) = surrogate_idx;
else
    surrogate_idx = [];
end

% --- Determine CAR channels for each block ---
skipped_blocks = [];
for ss = set_nrs
    these_channel_nrs = intersect(good_channels, out(ss).channels_set);
    if isempty(these_channel_nrs)
        skipped_blocks(end+1) = ss; %#ok<AGROW>
        out(ss).car_channels = [];
        continue;
    end
    
    % Get time mask and data in CAR window
    time_mask = ttt > car_timeint(1) & ttt < car_timeint(2);
    these_data = signaldata(these_channel_nrs, time_mask, :);
    
    % Concatenate trials (channels × all_time_samples)
    these_data_cat = reshape(these_data, size(these_data,1), []);
    
    % Variance-based selection
    chan_var = var(these_data_cat, [], 2);
    var_th = quantile(chan_var, perc_channels);
    chans_incl = setdiff(1:length(these_channel_nrs), find(chan_var > var_th));
    
    if isempty(chans_incl)
        warning('Block %d: No channels below variance threshold. Using all good channels.', ss);
        chans_incl = 1:length(these_channel_nrs);
    end
    out(ss).car_channels = these_channel_nrs(chans_incl);
end

% --- Print skipped blocks once ---
if ~isempty(skipped_blocks)
    warning('Skipped CAR for blocks: %s (no valid channels)', num2str(skipped_blocks));
end

% --- Apply CAR ---
for ss = set_nrs
    if isempty(out(ss).car_channels)
        continue;
    end
    car_data = mean(signaldata(out(ss).car_channels,:,:), 1);
    signaldata(out(ss).channels_set,:,:) = ...
        signaldata(out(ss).channels_set,:,:) - car_data;
end

% --- Remove surrogate channel if added ---
if preserve_rank && ~isempty(surrogate_idx)
    signaldata(surrogate_idx,:,:) = [];
end

% --- Rank after CAR ---
% r_after = rank(double(reshape(signaldata, size(signaldata,1), [])));
r_after = sum(eig(cov(double(signaldata(:,:)'))) > 1E-7);
if r_after < r_before
    warning('Effective data rank degraded during CAR: %d → %d', r_before, r_after);
end

end
