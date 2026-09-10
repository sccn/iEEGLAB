function [EEG, com] = ieeglab_stats_subject(EEG, opt)
% ieeglab_stats_subject() - Within-subject CCEP statistics via Canonical
%                           Response Parameterization (CRP).
%
% For each stimulation site and each recording channel, fits the canonical
% response shape across trials and reports response duration, explained
% variance and significance. Results are stored in EEG.ieeglab.stats.
%
% Usage:
%   EEG = ieeglab_stats_subject(EEG)          % GUI
%   EEG = ieeglab_stats_subject(EEG, opt)     % headless
%
% Options (fields of opt):
%   .crp_window   [1x2] response window in MILLISECONDS. Default [15 500].
%                 The lower bound should sit after the stimulation artifact.
%   .alpha        significance level. Default 0.05.
%   .min_trials   minimum trials per stimulation site to attempt a fit. Default 5.
%   .exclude_stim logical, skip channels that were stimulated for that site.
%                 Default true - a stimulated contact has no meaningful CCEP.
%   .correct      multiple-comparison correction across channels within a site:
%                 'fdr' (Benjamini-Hochberg, default), 'bonferroni', or 'none'.
%   .plot         logical, draw a summary figure. Default true in the GUI path,
%                 false when opt is supplied, so scripted runs stay headless.
%   .verbose      logical. Default true.
%
% Output:
%   EEG.ieeglab.stats.table   - one row per (site, channel) with tR, p, and
%                               explained variance
%   EEG.ieeglab.stats.crp     - the full run_CRP output per (site, channel)
%   com                       - command string for the EEGLAB history
%
% Method: Miller KJ, Muller KR, Ojeda Valencia G, Huang H, Gregg NM,
% Worrell GA, Hermes D (2023). Canonical Response Parameterization:
% Quantifying the structure of responses to single-pulse intracranial
% electrical brain stimulation. PLoS Comput Biol, 19(5):e1011105.
%
% Cedric Cannard, iEEGLAB, 2026

com = '';
if nargin < 1 || isempty(EEG)
    error('ieeglab_stats_subject:noData', 'No dataset provided.');
end

interactive = (nargin < 2 || isempty(opt));
if interactive, opt = struct(); end

% Figures default ON in the GUI path and OFF when opt is supplied, so scripted
% and CI runs stay headless without having to remember to disable plotting.
def = struct('crp_window',[15 500], 'alpha',0.05, 'min_trials',5, ...
             'exclude_stim',true, 'correct','fdr', 'plot',interactive, ...
             'verbose',true);
fn = fieldnames(def);
for ii = 1:numel(fn)
    if ~isfield(opt, fn{ii}) || isempty(opt.(fn{ii}))
        opt.(fn{ii}) = def.(fn{ii});
    end
end

% ---------------- preconditions ----------------
if EEG.trials < 2
    error('ieeglab_stats_subject:notEpoched', ...
        ['CRP needs epoched data with at least 2 trials (this dataset has %d). ' ...
         'Run "Preprocess iEEG data" with epoching enabled first.'], EEG.trials);
end
if ~isfield(EEG,'times') || numel(EEG.times) ~= size(EEG.data,2)
    error('ieeglab_stats_subject:noTimes', 'EEG.times does not match the data.');
end

% ---------------- GUI ----------------
if interactive
    uigeom = {[1.4 0.6] [1.4 0.6] [1.4 0.6] [1.4 0.6] [1] [1]};
    uilist = { ...
        {'Style','text','string','Response window (ms, after stimulation artifact)'}, ...
        {'Style','edit','string',num2str(opt.crp_window),'tag','win'}, ...
        {'Style','text','string','Significance level (alpha)'}, ...
        {'Style','edit','string',num2str(opt.alpha),'tag','alpha'}, ...
        {'Style','text','string','Minimum trials per stimulation site'}, ...
        {'Style','edit','string',num2str(opt.min_trials),'tag','minTr'}, ...
        {'Style','text','string','Correction across channels'}, ...
        {'Style','popupmenu','string',{'FDR (Benjamini-Hochberg)','Bonferroni','None'},'value',1,'tag','corr'}, ...
        {'Style','checkbox','string','Skip channels that were stimulated for that site','value',opt.exclude_stim,'tag','exclStim'}, ...
        {'Style','checkbox','string','Plot summary figure','value',1,'tag','doPlot'} };
    [~, ~, ok, res] = inputgui('geometry', uigeom, 'uilist', uilist, ...
        'title', 'Within-subject CCEP statistics (CRP)', 'helpcom', 'pophelp(''ieeglab_stats_subject'');');
    if isempty(ok) || isempty(res), return; end
    opt.crp_window   = str2num(res.win);      %#ok<ST2NM>
    opt.alpha        = str2double(res.alpha);
    opt.min_trials   = str2double(res.minTr);
    opt.exclude_stim = logical(res.exclStim);
    opt.plot         = logical(res.doPlot);
    opt.correct      = lower(strtok(res.corr{1}));
    if strcmpi(opt.correct,'fdr(benjamini-hochberg)'), opt.correct = 'fdr'; end
    if numel(opt.crp_window) < 2 || opt.crp_window(2) <= opt.crp_window(1)
        error('ieeglab_stats_subject:badWindow', ...
            'Response window must be [start stop] in ms with stop > start; got %s.', mat2str(opt.crp_window));
    end
end

% ---------------- window ----------------
tt_s = double(EEG.times(:))' / 1000;
pick = tt_s >= opt.crp_window(1)/1000 & tt_s <= opt.crp_window(2)/1000;
if nnz(pick) < 10
    error('ieeglab_stats_subject:shortWindow', ...
        ['Response window [%g %g] ms selects only %d samples of the epoch ' ...
         '(epoch spans [%g %g] ms at %g Hz). CRP needs at least 10.'], ...
        opt.crp_window(1), opt.crp_window(2), nnz(pick), EEG.times(1), EEG.times(end), EEG.srate);
end
t_win = tt_s(pick);

% ---------------- stimulation sites ----------------
labels = string({EEG.chanlocs.labels});
sites  = local_site_per_epoch(EEG, labels);
[uSites, ~, grp] = unique(sites);
isCCEP = any(uSites ~= "");

if opt.verbose
    if isCCEP
        fprintf('[CRP] %d stimulation site(s), %d channel(s), window [%g %g] ms.\n', ...
            numel(uSites), EEG.nbchan, opt.crp_window(1), opt.crp_window(2));
    else
        fprintf('[CRP] No stimulation sites found; treating all %d trials as one condition.\n', EEG.trials);
    end
end

% ---------------- fit ----------------
rows = {};
crpAll = struct('site',{},'channel',{},'parms',{},'projs',{});
nSkipped = 0;

for g = 1:numel(uSites)
    tr = find(grp == g);
    if numel(tr) < opt.min_trials
        nSkipped = nSkipped + 1;
        if opt.verbose
            fprintf('[CRP] site %-14s skipped (%d trials < %d)\n', char(uSites(g)), numel(tr), opt.min_trials);
        end
        continue
    end

    stimIdx = [];
    if opt.exclude_stim && uSites(g) ~= ""
        stimIdx = find(ismember(upper(labels), upper(split(uSites(g), '-'))'));
    end

    pvals = nan(EEG.nbchan,1);
    tmp   = cell(EEG.nbchan,1);
    for ch = 1:EEG.nbchan
        if ismember(ch, stimIdx), continue; end
        V = squeeze(EEG.data(ch, pick, tr));
        if size(V,2) < 2 || all(~isfinite(V(:))), continue; end
        try
            [parms, projs] = run_CRP(double(V), t_win, struct('verbose',false));
        catch ME
            warning('ieeglab_stats_subject:crpFailed', ...
                'CRP failed for site %s channel %s: %s', char(uSites(g)), labels(ch), ME.message);
            continue
        end
        if isempty(parms) || ~isfield(parms,'tR') || isempty(parms.tR), continue; end
        pvals(ch) = projs.p_value_tR;
        tmp{ch}   = struct('parms',parms,'projs',projs);
    end

    % multiple comparisons across channels within this site
    padj = local_correct(pvals, opt.correct);

    for ch = 1:EEG.nbchan
        if isempty(tmp{ch}), continue; end
        p = tmp{ch};
        rows(end+1,:) = { char(uSites(g)), char(labels(ch)), ...
            p.parms.tR*1000, mean(p.parms.expl_var,'omitnan'), ...
            projs_snr(p.parms), pvals(ch), padj(ch), padj(ch) < opt.alpha }; %#ok<AGROW>
        crpAll(end+1) = struct('site',char(uSites(g)),'channel',char(labels(ch)), ...
            'parms',p.parms,'projs',p.projs); %#ok<AGROW>
    end
end

if isempty(rows)
    warning('ieeglab_stats_subject:noResults', ...
        ['No CRP fits succeeded. Most often this means too few trials per stimulation site ' ...
         '(minimum is %d; lower it in the dialog) or a response window outside the epoch.'], opt.min_trials);
    EEG.ieeglab.stats = struct('table', table(), 'crp', crpAll, 'opt', opt);
    return
end

T = cell2table(rows, 'VariableNames', ...
    {'site','channel','tR_ms','explained_var','snr','p','p_adj','significant'});
T = sortrows(T, {'site','p_adj'});

EEG.ieeglab.stats = struct('table', T, 'crp', crpAll, 'opt', opt);

if opt.verbose
    nSig = sum(T.significant);
    fprintf('\n[CRP] %d significant channel-site pairs of %d tested (%s-corrected, alpha=%g).\n', ...
        nSig, height(T), opt.correct, opt.alpha);
    if nSkipped > 0
        fprintf('[CRP] %d site(s) skipped for having fewer than %d trials.\n', nSkipped, opt.min_trials);
    end
    disp(head(T(T.significant,:), min(10, nSig)));
end

% ---------------- optional figure ----------------
if opt.plot
    local_plot_summary(T, uSites, labels, opt);
end

com = sprintf('EEG = ieeglab_stats_subject(EEG, %s);', local_struct2str(opt));

end

% ========================== local helpers ==========================

function s = projs_snr(parms)
if isfield(parms,'Vsnr') && ~isempty(parms.Vsnr)
    s = mean(parms.Vsnr, 'omitnan');
else
    s = NaN;
end
end

function sites = local_site_per_epoch(EEG, labels)
% Canonical, order-independent stimulation site for each epoch.
N = EEG.trials;
sites = strings(1, N);
if ~isfield(EEG,'epoch') || isempty(EEG.epoch) || numel(EEG.epoch) ~= N, return; end
for i = 1:N
    t = EEG.epoch(i).eventtype;
    l = [];
    if isfield(EEG.epoch,'eventlatency'), l = EEG.epoch(i).eventlatency; end
    if iscell(t)
        if iscell(l) && numel(l) == numel(t)
            [~, k] = min(cellfun(@(x) abs(double(x(1))), l));
        else
            k = 1;
        end
        t = t{k};
    end
    if isnumeric(t), t = num2str(t); end
    t = strtrim(string(t));
    if t == "", continue; end
    parts = regexp(char(t), '[,;+\-\/\|\s]+', 'split');
    parts = parts(~cellfun(@isempty, parts));
    hit = parts(ismember(upper(parts), upper(labels)));
    if numel(hit) >= 2
        sites(i) = string(strjoin(sort(hit), '-'));
    else
        sites(i) = t;
    end
end
end

function padj = local_correct(p, method)
padj = p;
ok = ~isnan(p);
n = sum(ok);
if n == 0, return; end
switch lower(method)
    case 'bonferroni'
        padj(ok) = min(1, p(ok) * n);
    case 'fdr'
        [ps, ord] = sort(p(ok));
        q = ps(:) .* n ./ (1:n)';
        q = min(1, flipud(cummin(flipud(q))));   % enforce monotonicity
        tmp = nan(n,1); tmp(ord) = q;
        padj(ok) = tmp;
    case 'none'
        % leave as is
    otherwise
        warning('ieeglab_stats_subject:badCorrection', ...
            'Unknown correction "%s"; reporting uncorrected p-values.', method);
end
end

function local_plot_summary(T, uSites, labels, opt)
sig = T(T.significant,:);
figure('Color','w','Name','CRP within-subject summary','NumberTitle','off');

subplot(1,2,1);
if isempty(sig)
    text(0.5,0.5,'No significant responses','HorizontalAlignment','center'); axis off
else
    histogram(sig.tR_ms, max(8, round(sqrt(height(sig)))));
    xlabel('Response duration \tau_R (ms)','FontWeight','bold');
    ylabel('Significant channel-site pairs','FontWeight','bold');
    title(sprintf('%d significant of %d tested', height(sig), height(T)));
    box on
end

subplot(1,2,2);
M = nan(numel(labels), numel(uSites));
for r = 1:height(T)
    ci = find(strcmp(cellstr(labels), T.channel{r}), 1);
    si = find(uSites == string(T.site{r}), 1);
    if ~isempty(ci) && ~isempty(si) && T.significant(r)
        M(ci, si) = T.explained_var(r);
    end
end
imagesc(M, 'AlphaData', ~isnan(M)); set(gca,'Color',[.94 .94 .94]);
colormap(parula); c = colorbar; ylabel(c,'Explained variance');
set(gca,'XTick',1:numel(uSites),'XTickLabel',cellstr(uSites),'XTickLabelRotation',45, ...
        'YTick',1:numel(labels),'YTickLabel',cellstr(labels),'FontSize',8);
xlabel('Stimulation site','FontWeight','bold');
ylabel('Recording channel','FontWeight','bold');
title(sprintf('Significant CCEPs (%s, \\alpha=%g)', opt.correct, opt.alpha));
end

function s = local_struct2str(opt)
% Compact, re-runnable literal for the EEGLAB history.
f = fieldnames(opt); parts = {};
for i = 1:numel(f)
    v = opt.(f{i});
    if ischar(v),        parts{end+1} = sprintf('''%s'',''%s''', f{i}, v);
    elseif islogical(v), parts{end+1} = sprintf('''%s'',%d', f{i}, v);
    elseif isnumeric(v), parts{end+1} = sprintf('''%s'',%s', f{i}, mat2str(v));
    end %#ok<AGROW>
end
s = ['struct(' strjoin(parts, ', ') ')'];
end
