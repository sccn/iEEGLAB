function [img, info] = ieeglab_plot_ccep_matrix(M, metric, varargin)
% ieeglab_plot_ccep_matrix() - Draw the CCEP connectivity matrix.
%
% Usage:
%   ieeglab_plot_ccep_matrix(EEG.ieeglab.ccep_matrix)                 % responses
%   ieeglab_plot_ccep_matrix(M, 'amplitude_uv')                       % |N1| amplitude
%   ieeglab_plot_ccep_matrix(M, 'latency_ms')                         % N1 latency
%   [img, info] = ieeglab_plot_ccep_matrix(M, metric, 'draw', false)  % compute only
%
% Rows are stimulation sites, columns recording contacts. Grey cells were not
% measured (the stimulated contacts, or pairs with too few trials); white
% cells were measured with no significant response. For continuous metrics,
% only significant responses are coloured unless 'sig_only' is false.
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
    if strcmp(metric, 'amplitude_uv'), img = abs(img); end
    if o.sig_only, img(M.response ~= 1) = NaN; end
end

info = struct();
info.metric    = metric;
info.sites     = M.sites;
info.channels  = M.channels;
info.measured  = ~isnan(M.response);
info.n_sig     = nnz(M.response == 1);
info.n_tested  = nnz(~isnan(M.response));
labelsMetric = struct('response','Significant response', 'amplitude_uv','|N1| amplitude (\muV)', ...
    'latency_ms','N1 latency (ms)', 'tR_ms','CRP response duration \tau_R (ms)', ...
    'explained_var','CRP explained variance', 'p_adj','Adjusted p');
info.colorlabel = labelsMetric.(metric);
if isempty(o.title)
    info.title = sprintf('CCEP connectivity (%s): %d of %d tested pairs respond', ...
        upper(M.source), info.n_sig, info.n_tested);
else
    info.title = o.title;
end
if ~o.draw, return; end

% ---------- draw ----------
figure('Color','w', 'Name','CCEP connectivity matrix', 'NumberTitle','off');
ax = axes; hold(ax,'on');
grey = [0.82 0.82 0.82];
set(ax, 'Color', grey);               % shows through NaN (not measured)
if isResp
    shown = img; shown(isnan(M.response)) = NaN;
    imagesc(ax, shown, 'AlphaData', ~isnan(shown));
    colormap(ax, [1 1 1; 0.05 0.40 0.47]);
    clim(ax, [0 1]);
else
    imagesc(ax, img, 'AlphaData', ~isnan(img));
    % Not-measured cells grey, measured-but-not-significant cells white
    notSig = info.measured & isnan(img);
    if any(notSig(:))
        [r, c] = find(notSig);
        plot(ax, c, r, 's', 'MarkerSize', 4, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'none');
    end
    colormap(ax, parula);
    cb = colorbar(ax); ylabel(cb, info.colorlabel);
end
axis(ax, 'ij', 'tight');
set(ax, 'XTick', 1:numel(M.channels), 'XTickLabel', M.channels, 'XTickLabelRotation', 90, ...
        'YTick', 1:numel(M.sites), 'YTickLabel', M.sites, 'TickLength', [0 0], ...
        'FontSize', 8, 'Layer', 'top');
xlabel(ax, 'Recording contact', 'FontWeight', 'bold');
ylabel(ax, 'Stimulation site', 'FontWeight', 'bold');
title(ax, info.title, 'FontWeight', 'normal');
box(ax, 'on');
end
