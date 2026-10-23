function [com, info] = ieeglab_plot_crp(EEG, site, channel, varargin)
% ieeglab_plot_crp() - The CRP of one stimulation site and one recording contact,
%                      drawn as in the CRP authors' illustration script
%                      (CRP_illustration.m, Miller et al., 2023): the trials, how
%                      long the response lasts, the response shape and how well
%                      it fits each trial.
%
% Usage:
%   ieeglab_plot_crp(EEG)                          % dialog: choose the site, then the contact
%   ieeglab_plot_crp(EEG, 'RA3-RA4', 'RA7')
%   [com, info] = ieeglab_plot_crp(EEG, 'RA3-RA4', 'RA7', 'draw', false)   % values only
%
% Needs iEEGLAB > CCEP analysis with CRP. Panels:
%   A  Every trial (grey), their average over the whole window (black) and up to
%      the response duration tau_R (yellow), from before the pulse to the end of
%      the CRP window.
%   B  How similar the trials are for each candidate duration: every trial-pair
%      cross-projection (grey dots) and their mean (black). tau_R is the
%      duration where the mean peaks (red circle); the red ticks mark where the
%      mean is within 98% of that peak, i.e. how sharply tau_R is defined.
%   C  Every trial up to tau_R (black), the response shape scaled to that trial
%      (alpha x C, red) and what is left (residual, green), one row per trial.
%   D  For every trial, the fraction of its variance the shape explains (bars)
%      and the coefficient of determination (dots).
%
% info holds the values drawn: times_ms, trials, mean, tR_ms, tR_range_ms,
% shape_ms, shape, profile_ms, profile, S_all, alpha, explained_var, cod, and the
% row of the results table.
%
% Cedric Cannard, iEEGLAB, 2026

p = inputParser;
p.addParameter('draw', true);
p.addParameter('pre_ms', 100);       % how much of the pre-stimulus period panel A shows
p.parse(varargin{:});
o = p.Results;
com = ''; info = struct();

if ~isfield(EEG,'ieeglab') || ~isfield(EEG.ieeglab,'stats') || ~isfield(EEG.ieeglab.stats,'crp') ...
        || isempty(EEG.ieeglab.stats.crp)
    error('ieeglab_plot_crp:noCRP', 'No CRP results in this dataset: run iEEGLAB > CCEP analysis with CRP first.');
end
S = EEG.ieeglab.stats;
T = S.table;
allSites = string({S.crp.site});
allChans = string({S.crp.channel});

% ---- choose the pair
if nargin < 2 || isempty(site)
    us = unique(allSites, 'stable');
    nSig = arrayfun(@(s) nnz(T.significant & string(T.site) == s), us);
    [k, ok] = listdlg('ListString', compose('%s  (%d responding contacts)', us(:), nSig(:)), ...
        'SelectionMode', 'single', 'Name', 'iEEGLAB', 'PromptString', 'CRP: stimulation site', ...
        'ListSize', [320 300]);
    if ~ok || isempty(k), return; end
    site = us(k);
end
site = string(site);
if nargin < 3 || isempty(channel)
    rows = T(string(T.site) == site, :);
    if isempty(rows), error('ieeglab_plot_crp:unknownSite', 'No CRP result for site "%s".', site); end
    rows = sortrows(rows, {'significant','explained_var'}, {'descend','descend'});
    lab = compose('%s  %s  (tau_R %.0f ms, explained variance %.2f)', string(rows.channel), ...
        iff_str(rows.significant, 'responds', 'no response'), rows.tR_ms, rows.explained_var);
    [k, ok] = listdlg('ListString', lab, 'SelectionMode', 'single', 'Name', 'iEEGLAB', ...
        'PromptString', sprintf('CRP: contact (stimulation at %s)', site), 'ListSize', [420 360]);
    if ~ok || isempty(k), return; end
    channel = rows.channel{k};
end
channel = string(channel);
cs = ieeglab_canonical_site({char(site)}, {EEG.chanlocs.labels});
iRes = find(allSites == cs(1) & strcmpi(allChans, channel), 1);
if isempty(iRes)
    error('ieeglab_plot_crp:noPair', 'No CRP result for contact %s at site %s.', channel, site);
end
r = S.crp(iRes);
channel = string(r.channel);                        % as spelled in the dataset
row = T(string(T.site) == string(r.site) & string(T.channel) == channel, :);

% ---- the trials, as CRP saw them (the results must still describe these data)
sites = ieeglab_epoch_sites(EEG);
tr = find(sites == string(r.site));
ch = find(strcmp({EEG.chanlocs.labels}, char(channel)), 1);
if ~local_same_data(EEG, ch, tr, r.parms)
    error('ieeglab_plot_crp:stale', ['The CRP results no longer match the data (contact %s, %d trials ' ...
        'at %s now, %d when CRP ran): the dataset changed after the CCEP analysis. Run iEEGLAB > ' ...
        'CCEP analysis again.'], channel, numel(tr), r.site, size(r.parms.V_tR, 2));
end
win = S.opt.crp_window;
tshow = EEG.times >= -abs(o.pre_ms) & EEG.times <= min(win(2), EEG.times(end));
V = squeeze(double(EEG.data(ch, tshow, tr)));            % time x trials
if isvector(V), V = V(:); end
tR_ms = r.parms.tR * 1000;
prof_ms = r.projs.proj_tpts(:)' * 1000;
prof = r.projs.mean_proj_profile(:)';
[pk, iPk] = max(prof);
lo = iPk; while lo > 1 && prof(lo-1) >= 0.98 * pk, lo = lo - 1; end
hi = iPk; while hi < numel(prof) && prof(hi+1) >= 0.98 * pk, hi = hi + 1; end
info = struct('site', char(r.site), 'channel', char(channel), 'times_ms', double(EEG.times(tshow)), ...
    'trials', V, 'mean', mean(V, 2, 'omitnan'), 'tR_ms', tR_ms, 'tR_range_ms', prof_ms([lo hi]), ...
    'shape_ms', r.parms.parms_times(:)' * 1000, 'shape', r.parms.C(:)', ...
    'profile_ms', prof_ms, 'profile', prof, 'S_all', r.projs.S_all, 'alpha', r.parms.al(:)', ...
    'explained_var', r.parms.expl_var(:)', 'cod', r.parms.cod(:)', 'result', row);
com = sprintf('ieeglab_plot_crp(EEG, ''%s'', ''%s'');', char(r.site), char(channel));
if ~o.draw, return; end

% ---- figure
verdict = 'no significant response';
if ~isempty(row) && row.significant, verdict = 'significant response'; end
pTxt = '';
if ~isempty(row), pTxt = sprintf(', adjusted p = %.3g', row.p_adj); end
f = ieeglab_figure('Name', sprintf('CRP %s -> %s', r.site, channel), 'Position', [60 60 1150 720]);
yellow = [0.988 0.831 0];  red = [0.85 0.1 0.1];  green = [0.1 0.6 0.2];
inTR = info.times_ms >= win(1) & info.times_ms <= tR_ms;

% A: trials, average, average up to tau_R
ax = subplot(2, 3, [1 2]); hold(ax, 'on');
plot(ax, info.times_ms, V, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
yline(ax, 0, 'Color', [0.3 0.3 0.3]);
plot(ax, info.times_ms, info.mean, 'k', 'LineWidth', 2);
plot(ax, info.times_ms(inTR), info.mean(inTR), 'k', 'LineWidth', 4);
plot(ax, info.times_ms(inTR), info.mean(inTR), 'Color', yellow, 'LineWidth', 3);
xline(ax, 0, 'k:');
xline(ax, tR_ms, '--', sprintf('\\tau_R = %.0f ms', tR_ms), 'Color', red, 'LabelVerticalAlignment', 'bottom');
ylim(ax, local_ylim(V)); xlim(ax, info.times_ms([1 end]));
xlabel(ax, 'Time after the pulse (ms)'); ylabel(ax, 'Voltage (\muV)');
title(ax, sprintf('%s stimulated, recorded at %s: %d trials, %s%s', r.site, channel, numel(tr), verdict, pTxt), ...
    'FontWeight', 'normal', 'Interpreter', 'none');
box(ax, 'on');

% B: cross-projections and their mean against duration
ax = subplot(2, 3, 3); hold(ax, 'on');
Sd = info.S_all;
plot(ax, repmat(prof_ms, size(Sd, 1), 1), Sd, '.', 'Color', [0.3 0.3 0.3], 'MarkerSize', 3);
yline(ax, 0, 'Color', [0.6 0.6 0.6], 'LineWidth', 2);
plot(ax, prof_ms, prof, 'k', 'LineWidth', 2);
plot(ax, prof_ms(iPk), pk, 'o', 'Color', red, 'MarkerSize', 9, 'LineWidth', 1.5);
yl = ylim(ax); tick = 0.04 * diff(yl);
plot(ax, prof_ms([lo lo]), pk + [-tick tick], '-', 'Color', red, 'LineWidth', 1.5);
plot(ax, prof_ms([hi hi]), pk + [-tick tick], '-', 'Color', red, 'LineWidth', 1.5);
xlabel(ax, 'Candidate duration (ms)'); ylabel(ax, 'Projection magnitude (\muV s^{1/2})');
title(ax, sprintf('Trial-to-trial similarity; peak at %.0f ms', tR_ms), 'FontWeight', 'normal');
box(ax, 'on');

% C: each trial up to tau_R, alpha x C, residual
ax = subplot(2, 3, [4 5]); hold(ax, 'on');
Vt = r.parms.V_tR; fit = r.parms.C(:) * r.parms.al(:)'; res = r.parms.ep;
K = size(Vt, 2);
step = 1.2 * max(abs(Vt(:)));
if ~(step > 0), step = 1; end
tt = info.shape_ms;
for k = 1:K
    off = -(k - 1) * step;
    plot(ax, tt, Vt(:, k) + off, 'k', 'LineWidth', 1.5);
    plot(ax, tt, fit(:, k) + off, 'Color', red, 'LineWidth', 1);
    plot(ax, tt, res(:, k) + off, 'Color', green, 'LineWidth', 1);
end
set(ax, 'YTick', fliplr(-(0:K-1) * step), 'YTickLabel', fliplr(1:K));
xlim(ax, [tt(1) tt(end)]); ylim(ax, [-(K - 0.4) * step, 0.6 * step]);
xlabel(ax, 'Time after the pulse (ms)'); ylabel(ax, 'Trial');
legend(ax, {'trial', '\alpha \times shape', 'residual'}, 'Location', 'eastoutside', 'Box', 'off');
title(ax, sprintf('Each trial up to \\tau_R (row spacing %.0f \\muV)', step), 'FontWeight', 'normal');
box(ax, 'on');

% D: fit per trial
ax = subplot(2, 3, 6); hold(ax, 'on');
bar(ax, info.explained_var, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none');
plot(ax, info.cod, 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
yline(ax, mean(info.explained_var, 'omitnan'), '-', sprintf('mean %.2f', mean(info.explained_var, 'omitnan')));
ylim(ax, [min(0, min(info.cod)) 1]); xlim(ax, [0.4 K + 0.6]);
xlabel(ax, 'Trial'); ylabel(ax, 'Fraction of variance');
legend(ax, {'explained variance', 'coefficient of determination'}, 'Location', 'southoutside', 'Box', 'off');
title(ax, 'Fit of the shape to each trial', 'FontWeight', 'normal');
box(ax, 'on');
set(f, 'UserData', info);
end

function ok = local_same_data(EEG, ch, tr, parms)
% The trials of this contact and site, over the samples CRP kept, must be the
% ones CRP was fitted on (trials with NaN were left out by run_CRP).
ok = false;
if isempty(ch) || isempty(tr) || isempty(parms.V_tR) || isempty(parms.parms_times), return; end
i0 = find(abs(double(EEG.times) - parms.parms_times(1) * 1000) < 1e-6, 1);
n = size(parms.V_tR, 1);
if isempty(i0) || i0 + n - 1 > EEG.pnts, return; end
X = reshape(double(EEG.data(ch, i0:i0+n-1, tr)), n, numel(tr));
X = X(:, all(isfinite(X), 1));
if ~isequal(size(X), size(parms.V_tR)), return; end
ok = max(abs(X(:) - parms.V_tR(:))) <= 1e-6 * max(1, max(abs(parms.V_tR(:))));
end

function yl = local_ylim(V)
v = sort(abs(V(isfinite(V))));
m = 1;
if ~isempty(v), m = v(max(1, ceil(0.995 * numel(v)))); end
if ~(m > 0), m = 1; end
yl = [-m m] * 1.05;
end

function s = iff_str(tf, a, b)
s = repmat(string(b), size(tf));
s(logical(tf)) = string(a);
end
