function plot_ccep(data, timevec, chan_names, view_type, chan_idx, trim_prop)
% plot_ccep - Plot Cortico-Cortical Evoked Potentials (CCEPs) from iEEG/EEG data
%
% Usage:
%   plot_ccep(signal, timevec, chan_names, view_type, chan_idx, trim_prop)
%
% Inputs:
%   data      - Data array [nChannels x nTimePoints x nTrials]
%   timevec     - Time vector [1 x nTimePoints] in seconds
%   chan_names  - (Optional) cell array of channel names for labeling
%   view_type   - 'all' (heatmap of all channels) or 'single' (plot one channel)
%   chan_idx    - Channel index (required if view_type = 'single')
%   trim_prop   - Proportion to trim for trimmed mean (e.g., 0.10 for 10%)
%
% Description:
%   This function visualizes CCEP responses across channels or for a single
%   channel using trimmed mean statistics to reduce outlier influence.
%   For 'all' view, a heatmap shows amplitude across channels and time.
%   For 'single' view, all trial waveforms are plotted in gray with the
%   trimmed mean overlaid in red.
%
% Cedric Cannard © iEEGLAB Plugin, 2025

% Argument positions: data=1 timevec=2 chan_names=3 view_type=4 chan_idx=5 trim_prop=6
% The nargin thresholds below must match those positions. They previously did
% not, so plot_ccep(d,t,n,'single') silently drew a heatmap and the valid
% 5-argument call plot_ccep(d,t,n,'single',3) raised a spurious error.
if nargin < 6 || isempty(trim_prop)
    trim_prop = 0.10; % default 10%
end
if nargin < 4 || isempty(view_type)
    view_type = 'all';
end
if nargin < 5, chan_idx = []; end
if strcmpi(view_type, 'single') && isempty(chan_idx)
    error('plot_ccep:noChanIdx', 'chan_idx must be provided when view_type = ''single''.');
end

switch lower(view_type)
    case 'all'
        % --- Collapse trials with trimmed mean ---
        data = squeeze(trimmean(data, trim_prop*100, 3)); % chan × time

        % --- Main plot ---
        figure('Color', 'w');
        imagesc(timevec, 1:size(data,1), data);
        % colormap(gca, jet);
        colormap("parula")  % bone hot winter summer spring turbo
        c = colorbar;
        ylabel(c, 'Amplitude (\muV)', 'FontWeight','bold', ...
            'FontSize',14, 'Rotation',-90);

        % --- Y-ticks (channels) ---
        if exist('chan_names','var') && ~isempty(chan_names)
            skip = 3;
            set(gca, 'YTick', 1:skip:numel(chan_names), ...
                'YTickLabel', chan_names(1:skip:end), ...
                'FontSize',10,'FontWeight','bold');
        else
            ylabel('Channels', 'FontSize', 14, 'FontWeight', 'bold');
        end

        % --- X-ticks (time) ---
        if exist('timevec','var') && ~isempty(timevec)
            skip = ceil(numel(timevec)/8);
            xt = timevec(1:skip:end);
            set(gca, 'XTick', xt, ...
                'XTickLabel', arrayfun(@(x) sprintf('%g', x), xt, 'UniformOutput', false), ...
                'FontSize', 10, 'FontWeight', 'bold');
        end

        % --- Labels ---
        xlabel('Time (s)', 'FontSize', 14, 'FontWeight','bold');
        ylabel('Channels', 'FontSize', 14, 'FontWeight','bold');

        % --- Symmetric color limits (robust to outliers) ---
        vals = data(:);
        lo = prctile(vals, 1);   % 1st percentile
        hi = prctile(vals, 99);  % 99th percentile
        clim_val = max(abs([lo hi]));
        clim([-clim_val clim_val]);

        % --- Styling ---
        set(gca, 'YDir', 'normal', 'LineWidth',1);
        xline(0, '--k', 'LineWidth', 1.5); % stim onset
        title(sprintf('CCEP Map (%.0f%% Trimmed Mean)', trim_prop*100), ...
            'FontSize', 18, 'FontWeight', 'bold');

    case 'single'

        if ndims(data)==3 && size(data,2)~=numel(timevec) && size(data,3)==numel(timevec)
            data = permute(data, [1 3 2]);  % enforce [chan x time x trials]
        end

        % Extract channel trials: time × trials
        chan_trials = squeeze(data(chan_idx,:,:))'; % trials × time
        chan_mean = trimmean(chan_trials, trim_prop*100, 1);

        figure('Color', 'w'); hold on;
        plot(timevec, chan_trials', 'Color', [0 0 0 0.3], 'LineWidth', 0.5);
        plot(timevec, chan_mean, 'Color', [0.8 0.1 0.1], 'LineWidth', 2);
        xline(0, '--k', 'LineWidth', 1);
        xlim([min(timevec) max(timevec)]);
        xlabel('Time (s)', 'FontSize', 14, 'FontWeight', 'bold');
        ylabel('Amplitude (\muV)', 'FontSize', 14, 'FontWeight', 'bold');

        if exist('chan_names','var') && ~isempty(chan_names)
            ch_label = chan_names{chan_idx};
        else
            ch_label = sprintf('Channel %d', chan_idx);
        end
        % no trim_prop in title here
        title(sprintf('%s - CCEP', ch_label), 'FontSize', 16, 'FontWeight', 'bold');
        box on;

    otherwise
        error('view_type must be ''all'' or ''single''.');
end
end
