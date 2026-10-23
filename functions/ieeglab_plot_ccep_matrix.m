function [img, info] = ieeglab_plot_ccep_matrix(M, metric, varargin)
% ieeglab_plot_ccep_matrix() - Draw the CCEP connectivity matrix.
%
% Usage:
%   ieeglab_plot_ccep_matrix(EEG.ieeglab.ccep_matrix)                 % responses
%   ieeglab_plot_ccep_matrix(M, 'amplitude_uv')                       % N1 amplitude (signed)
%   ieeglab_plot_ccep_matrix(M, 'latency_ms')                         % N1 latency
%   ieeglab_plot_ccep_matrix(M, 'explained_var')                      % CRP explained variance
%   ieeglab_plot_ccep_matrix(M, 'tR_ms')                              % CRP response duration
%   [img, info] = ieeglab_plot_ccep_matrix(M, metric, 'draw', false)  % compute only
%
% Rows are stimulation sites, columns recording contacts. Each cell is one
% site-contact pair: for 'response', teal = significant response (N1: the
% detection rule or permutation test of ieeglab_detect_n1; CRP: the permutation
% test of ieeglab_stats_subject), white = tested with no significant response,
% grey = not tested (the stimulated contacts, bad channels, too few trials).
% For continuous metrics, only significant responses are coloured unless
% 'sig_only' is false; white and grey keep their meaning.
%
% The 'draw', false form returns exactly what would be drawn, so the figure
% logic can be tested without a display.
%
% Cedric Cannard, iEEGLAB, 2026

if nargin < 2 || isempty(metric), metric = 'response'; end
p = inputParser;
p.addParameter('draw', true);
p.addParameter('sig_only', true);
p.addParameter('title', '');
p.parse(varargin{:});
o = p.Results;

if isfield(M,'ieeglab') && isfield(M.ieeglab,'ccep_matrix'), M = M.ieeglab.ccep_matrix; end
if ~isfield(M, metric)
    avail = intersect({'response','amplitude_uv','latency_ms','tR_ms','explained_var','p_adj'}, fieldnames(M));
    error('ieeglab_plot_ccep_matrix:badMetric', ...
        'Metric "%s" is not in this matrix. Available: %s', metric, strjoin(avail, ', '));
end

img = double(M.(metric));
isResp = strcmp(metric, 'response');
if ~isResp
    if o.sig_only, img(M.response ~= 1) = NaN; end
end

info = struct();
info.metric    = metric;
info.sites     = M.sites;
info.channels  = M.channels;
info.measured  = ~isnan(M.response);
info.n_sig     = nnz(M.response == 1);
info.n_tested  = nnz(~isnan(M.response));
labelsMetric = struct('response','Significant response', 'amplitude_uv','N1 amplitude (\muV, signed)', ...
    'latency_ms','N1 latency (ms after stimulus)', 'tR_ms','CRP \tau_R (ms after stimulus)', ...
    'explained_var','CRP explained variance (fraction, 0-1)', 'p_adj','Adjusted p');
info.colorlabel = labelsMetric.(metric);
if isempty(o.title)
    info.title = sprintf('CCEP connectivity (%s%s): %d of %d tested pairs respond', ...
        upper(M.source), local_method(M), info.n_sig, info.n_tested);
else
    info.title = o.title;
end
if ~o.draw, return; end

% ---------- draw ----------
ieeglab_figure('Name','CCEP connectivity matrix');
ax = axes; hold(ax,'on');
grey = [0.82 0.82 0.82];
set(ax, 'Color', grey);                         % shows through: not tested
% White layer under every tested pair; the metric is drawn on top, so white
% left visible means tested without a (significant) value.
image(ax, repmat(reshape([1 1 1], 1, 1, 3), size(M.response)), 'AlphaData', double(info.measured));
teal = [0.05 0.40 0.47];
if isResp
    shown = img; shown(M.response ~= 1) = NaN;
    imagesc(ax, shown, 'AlphaData', ~isnan(shown));
    colormap(ax, [teal; teal]);
    h = [patch(ax, NaN, NaN, teal), patch(ax, NaN, NaN, [1 1 1]), patch(ax, NaN, NaN, grey)];
    legend(ax, h, {'significant response', 'tested, no response', 'not tested'}, ...
        'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');
else
    imagesc(ax, img, 'AlphaData', ~isnan(img));
    v = img(isfinite(img));
    switch metric
        case 'amplitude_uv'                     % signed: diverging, white at 0
            m = max([abs(v); eps]);
            colormap(ax, local_diverging(256)); caxis(ax, [-m m]); %#ok<CAXIS>
        case 'explained_var'
            colormap(ax, parula); caxis(ax, [0 1]); %#ok<CAXIS>
        otherwise
            colormap(ax, parula);
            if ~isempty(v) && max(v) > min(v), caxis(ax, [min(v) max(v)]); end %#ok<CAXIS>
    end
    cb = colorbar(ax); ylabel(cb, info.colorlabel);
end
axis(ax, 'ij', 'tight');
xs = local_ticks(numel(M.channels), 80); ys = local_ticks(numel(M.sites), 60);
set(ax, 'XTick', xs, 'XTickLabel', M.channels(xs), 'XTickLabelRotation', 90, ...
        'YTick', ys, 'YTickLabel', M.sites(ys), 'TickLength', [0 0], ...
        'FontSize', 8, 'Layer', 'top', 'TickLabelInterpreter', 'none');
xlabel(ax, 'Recording contact', 'FontWeight', 'bold');
ylabel(ax, 'Stimulation site', 'FontWeight', 'bold');
title(ax, {info.title, info.colorlabel}, 'FontWeight', 'normal');
box(ax, 'on');
end

function k = local_ticks(n, maxLabels)
k = 1:max(1, ceil(n / maxLabels)):n;
end

function s = local_method(M)
% ", permutation, fdr, alpha 0.05" from the options the source table was made with
s = '';
if ~isfield(M,'source_opt') || ~isstruct(M.source_opt), return; end
o = M.source_opt;
if isfield(o,'method') && ~isempty(o.method), s = [s ', ' char(string(o.method))]; end
if isfield(o,'crp_test') && strcmpi(char(string(o.crp_test)),'ttest'), s = [s ', t-test at tau_R']; end
if isfield(o,'correct') && ~isempty(o.correct), s = [s ', ' char(string(o.correct))]; end
if isfield(o,'alpha') && ~isempty(o.alpha) && ~(isfield(o,'method') && strcmpi(char(string(o.method)),'sd'))
    s = [s sprintf(', alpha %g', o.alpha)];
end
end

function cmap = local_diverging(n)
% blue - white - red, white at zero
h = floor(n/2);
up = linspace(0, 1, h)';
dn = flipud(linspace(0, 1, n-h)');
cmap = [[up up ones(h,1)]; [ones(n-h,1) dn dn]];
end
