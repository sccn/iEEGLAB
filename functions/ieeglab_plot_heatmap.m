function [com, info] = ieeglab_plot_heatmap(EEG, varargin)
% ieeglab_plot_heatmap() - Average over all epochs, every channel, as a
%                          channels x time heatmap. Works at any step once the
%                          data are epoched (before or after filtering,
%                          re-referencing, baseline removal, CCEP analysis).
%
% Usage:
%   ieeglab_plot_heatmap(EEG)
%   ieeglab_plot_heatmap(EEG, 'site', 'RA1-RA2', 'clim', [-200 200])
%   [com, info] = ieeglab_plot_heatmap(EEG, 'draw', false)   % values only
%
% Options:
%   'site'        average only the epochs of this stimulation site (CCEP).
%                 Default: all epochs.
%   'exclude_stim' CCEP: leave out each contact's own stimulation epochs, where
%                 it carries the artifact. Default true.
%   'clim'        colour limits in uV. Default: symmetric, 99th percentile of |values|.
%   'draw'        false returns the values in info without plotting.
%
% info.mean (channels x time, uV), info.times (ms), info.labels, info.n_trials
% (per channel, after exclusions).
%
% Cedric Cannard, iEEGLAB, 2026

p = inputParser;
p.addParameter('site', '');
p.addParameter('exclude_stim', true);
p.addParameter('clim', []);
p.addParameter('draw', true);
p.parse(varargin{:});
o = p.Results;

ieeglab_require_data(EEG, 'ieeglab_plot_heatmap');
if EEG.trials < 2
    error('ieeglab_plot_heatmap:notEpoched', ...
        'The heatmap averages epochs: epoch the data first (iEEGLAB > Preprocess iEEG data).');
end
labels = {EEG.chanlocs.labels};
tr = 1:EEG.trials;
isCCEP = strcmp(ieeglab_detect_mode(EEG), 'ccep');
what = 'all epochs';
if ~isempty(o.site)
    sites = ieeglab_epoch_sites(EEG);
    cs = ieeglab_canonical_site({char(o.site)}, EEG.chanlocs);
    tr = find(sites == cs(1));
    if isempty(tr)
        error('ieeglab_plot_heatmap:unknownSite', 'No epoch was stimulated at "%s". Sites: %s', ...
            char(o.site), strjoin(cellstr(unique(sites)), ', '));
    end
    what = sprintf('stimulation of %s', char(o.site));
end
D = EEG.data(:, :, tr);                             % in the data's class: a full recording is many GB
note = '';
if isCCEP && o.exclude_stim
    [~, stimIdx] = ieeglab_epoch_sites(EEG);
    for k = 1:numel(tr)
        sIdx = stimIdx{tr(k)};
        if ~isempty(sIdx), D(sIdx, :, k) = NaN; end
    end
    note = ', each contact''s own stimulation epochs excluded';
end
nTr = sum(any(isfinite(D), 2), 3);
avg = mean(D, 3, 'double', 'omitnan');

info = struct('mean', avg, 'times', double(EEG.times), 'labels', {labels}, 'n_trials', nTr(:)');
com = sprintf('ieeglab_plot_heatmap(EEG%s);', local_args(o));
if ~o.draw, return; end

cl = o.clim;
if isempty(cl)
    v = sort(abs(avg(isfinite(avg))));
    m = 1;
    if ~isempty(v), m = v(max(1, ceil(0.99 * numel(v)))); end   % 99th percentile, no toolbox
    if ~(m > 0), m = 1; end
    cl = [-m m];
end

ieeglab_figure('Name','iEEGLAB average heatmap');
h = imagesc(info.times, 1:numel(labels), avg);
set(h, 'AlphaData', isfinite(avg));                 % channels with no epoch left: grey
set(gca, 'Color', [0.8 0.8 0.8], 'YDir', 'reverse', 'TickLabelInterpreter', 'none', 'Layer', 'top');
colormap(gca, parula(256));                         % limits symmetric: 0 uV is the middle (green)
caxis(cl); %#ok<CAXIS> works before R2022a
cb = colorbar; ylabel(cb, 'Amplitude (\muV)');
step = max(1, ceil(numel(labels) / 60));            % at most ~60 labels
set(gca, 'YTick', 1:step:numel(labels), 'YTickLabel', labels(1:step:end), 'FontSize', 8);
if info.times(1) < 0 && info.times(end) > 0, xline(0, 'k-', 'LineWidth', 1); end
xlabel('Time (ms)'); ylabel('Channel');
title({sprintf('Average of %d epochs (%s%s)', numel(tr), what, note), ...
       'grey: no epoch left for this channel'}, 'Interpreter', 'none', 'FontWeight', 'normal');
end


function s = local_args(o)
s = '';
if ~isempty(o.site), s = [s sprintf(', ''site'', ''%s''', strrep(char(o.site), '''', ''''''))]; end
if ~o.exclude_stim, s = [s ', ''exclude_stim'', false']; end
if ~isempty(o.clim), s = [s sprintf(', ''clim'', %s', mat2str(o.clim))]; end
end
