function [EEG, com] = ieeglab_stats_subject(EEG, opt)
% ieeglab_stats_subject() - Within-subject CCEP analysis in one step: N1
%                           detection, Canonical Response Parameterization,
%                           the connectivity matrix, and export.
%
% Usage:
%   EEG = ieeglab_stats_subject(EEG)          % dialog
%   EEG = ieeglab_stats_subject(EEG, opt)     % headless
%
% Every stage can be switched on or off (dialog checkboxes, or these fields):
%   .run_n1        N1 amplitude and latency per site x contact. Default: on for CCEP data.
%   .run_crp       Canonical Response Parameterization. Default: on.
%   .run_matrix    connectivity matrix, sites x contacts. Default: on for CCEP data.
%   .matrix_source 'n1' (default) or 'crp' - which detector defines a response
%   .export_dir    write TSV/JSON/MAT results there. '' (default) = no export.
%   .plot          CRP summary and connectivity-matrix figures. Default: on in
%                  the dialog, off when opt is passed, so scripts stay headless.
%
% Shared:  .alpha (0.05)  .correct ('fdr'|'bonferroni'|'none')  .min_trials (5)
%          .exclude_stim (true)  .verbose (true)
% CRP:     .crp_window  [15 500] ms, after the stimulation artifact
% N1:      .n1_window   [10 100] ms   .n1_baseline [-500 -10] ms
%          .n1_method   'sd' | 'permutation'   .n1_threshold 3.4
%
% Results:
%   EEG.ieeglab.stats.table   CRP: tR, explained variance, p, per site x contact
%   EEG.ieeglab.stats.crp     full run_CRP output per site x contact
%   EEG.ieeglab.n1            N1 results (see ieeglab_detect_n1)
%   EEG.ieeglab.ccep_matrix   connectivity matrix (see ieeglab_ccep_matrix)
%
% Methods:
%   CRP - Miller KJ, et al. (2023). PLoS Comput Biol 19(5):e1011105.
%   N1  - van Blooijs D, et al. (2018). Hum Brain Mapp 39(11):4611-4622.
%
% Cedric Cannard, iEEGLAB, 2026

com = '';
if nargin < 1 || isempty(EEG)
    error('ieeglab_stats_subject:noData', 'No dataset provided.');
end
interactive = (nargin < 2 || isempty(opt));
if interactive, opt = struct(); end

isCCEP = strcmp(ieeglab_detect_mode(EEG), 'ccep');
def = struct('run_n1',isCCEP, 'run_crp',true, 'run_matrix',isCCEP, 'matrix_source','n1', ...
    'export_dir','', 'plot',interactive, 'alpha',0.05, 'correct','fdr', 'min_trials',5, ...
    'exclude_stim',true, 'verbose',true, 'crp_window',[15 500], ...
    'n1_window',[10 100], 'n1_baseline',[-500 -10], 'n1_method','sd', 'n1_threshold',3.4);
fn = fieldnames(def);
for ii = 1:numel(fn)
    if ~isfield(opt, fn{ii}) || (isempty(opt.(fn{ii})) && ~ischar(def.(fn{ii})))
        opt.(fn{ii}) = def.(fn{ii});
    end
end

% ---------------- preconditions ----------------
if EEG.trials < 2
    error('ieeglab_stats_subject:notEpoched', ...
        ['CCEP analysis needs epoched data with at least 2 trials (this dataset has %d). ' ...
         'Run "Preprocess iEEG data" with epoching enabled first.'], EEG.trials);
end
if ~isfield(EEG,'times') || numel(EEG.times) ~= size(EEG.data,2)
    error('ieeglab_stats_subject:noTimes', 'EEG.times does not match the data.');
end

% ---------------- dialog ----------------
if interactive
    opt = local_dialog(EEG, opt, isCCEP);
    if isempty(opt), return; end
end

% ---------------- CRP ----------------
if opt.run_crp
    EEG = local_crp(EEG, opt);
end

% ---------------- N1 ----------------
if opt.run_n1
    [EEG, ~] = ieeglab_detect_n1(EEG, struct('n1_window', opt.n1_window, 'baseline', opt.n1_baseline, ...
        'method', opt.n1_method, 'threshold', opt.n1_threshold, 'alpha', opt.alpha, ...
        'correct', opt.correct, 'min_trials', opt.min_trials, 'exclude_stim', opt.exclude_stim, ...
        'verbose', opt.verbose));
end

% ---------------- connectivity matrix ----------------
if opt.run_matrix
    src = lower(char(opt.matrix_source));
    if strcmp(src,'crp') && ~opt.run_crp
        warning('ieeglab_stats_subject:matrixNeedsCRP', 'Matrix source is CRP but CRP is off; using N1.');
        src = 'n1';
    end
    try
        EEG = ieeglab_ccep_matrix(EEG, struct('source', src, 'compute_missing', true, 'plot', false, ...
            'verbose', opt.verbose, 'n1_window', opt.n1_window, 'baseline', opt.n1_baseline, ...
            'alpha', opt.alpha, 'correct', opt.correct, 'min_trials', opt.min_trials));
    catch ME
        warning('ieeglab_stats_subject:matrixFailed', 'Connectivity matrix not built: %s', ME.message);
    end
end

% ---------------- figures ----------------
if opt.plot
    if opt.run_crp && isfield(EEG.ieeglab,'stats') && ~isempty(EEG.ieeglab.stats.table)
        local_plot_crp(EEG.ieeglab.stats.table, opt);
    end
    if opt.run_matrix && isfield(EEG.ieeglab,'ccep_matrix') && ~isempty(EEG.ieeglab.ccep_matrix)
        ieeglab_plot_ccep_matrix(EEG.ieeglab.ccep_matrix);
    end
end

% ---------------- export ----------------
if ~isempty(opt.export_dir)
    ieeglab_export(EEG, opt.export_dir, struct('verbose', opt.verbose));
end

com = sprintf('EEG = ieeglab_stats_subject(EEG, %s);', local_struct2str(rmfield(opt, 'plot')));
end

% ========================== stages ==========================

function EEG = local_crp(EEG, opt)
tt_s = double(EEG.times(:))' / 1000;
pick = tt_s >= opt.crp_window(1)/1000 & tt_s <= opt.crp_window(2)/1000;
if nnz(pick) < 10
    error('ieeglab_stats_subject:shortWindow', ...
        ['CRP window [%g %g] ms selects only %d samples of the epoch ' ...
         '(epoch spans [%g %g] ms at %g Hz). CRP needs at least 10.'], ...
        opt.crp_window(1), opt.crp_window(2), nnz(pick), EEG.times(1), EEG.times(end), EEG.srate);
end
t_win = tt_s(pick);

labels = string({EEG.chanlocs.labels});
sites  = ieeglab_epoch_sites(EEG);
[uSites, ~, grp] = unique(sites);
isBadCh = false(1, EEG.nbchan);
if isfield(EEG.chanlocs,'status')
    isBadCh = cellfun(@(x) ~isempty(x) && strcmpi(char(x),'bad'), {EEG.chanlocs.status});
end
if opt.verbose
    fprintf('[CRP] %d site(s)/condition(s), %d channel(s), window [%g %g] ms.\n', ...
        numel(uSites), EEG.nbchan, opt.crp_window(1), opt.crp_window(2));
end

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
        stimIdx = find(ismember(upper(labels), upper(ieeglab_site_tokens(uSites(g)))));
    end
    pvals = nan(EEG.nbchan,1);
    tmp   = cell(EEG.nbchan,1);
    for ch = 1:EEG.nbchan
        if ismember(ch, stimIdx) || isBadCh(ch), continue; end
        V = squeeze(EEG.data(ch, pick, tr));
        if size(V,2) < 2 || all(~isfinite(V(:))), continue; end
        try
            [parms, projs] = run_CRP(double(V), t_win, struct('verbose',false));
        catch ME
            warning('ieeglab_stats_subject:crpFailed', 'CRP failed for %s / %s: %s', ...
                char(uSites(g)), labels(ch), ME.message);
            continue
        end
        if isempty(parms) || ~isfield(parms,'tR') || isempty(parms.tR), continue; end
        pvals(ch) = projs.p_value_tR;
        tmp{ch}   = struct('parms',parms,'projs',projs);
    end
    padj = local_correct(pvals, opt.correct);
    for ch = 1:EEG.nbchan
        if isempty(tmp{ch}), continue; end
        p = tmp{ch};
        snr = NaN;
        if isfield(p.parms,'Vsnr') && ~isempty(p.parms.Vsnr), snr = mean(p.parms.Vsnr,'omitnan'); end
        rows(end+1,:) = { char(uSites(g)), char(labels(ch)), p.parms.tR*1000, ...
            mean(p.parms.expl_var,'omitnan'), snr, pvals(ch), padj(ch), padj(ch) < opt.alpha }; %#ok<AGROW>
        crpAll(end+1) = struct('site',char(uSites(g)),'channel',char(labels(ch)), ...
            'parms',p.parms,'projs',p.projs); %#ok<AGROW>
    end
end

if isempty(rows)
    warning('ieeglab_stats_subject:noResults', ...
        ['No CRP fits succeeded. Usually too few trials per stimulation site ' ...
         '(minimum %d) or a response window outside the epoch.'], opt.min_trials);
    T = table();
else
    T = cell2table(rows, 'VariableNames', ...
        {'site','channel','tR_ms','explained_var','snr','p','p_adj','significant'});
    T = sortrows(T, {'site','p_adj'});
end
EEG.ieeglab.stats = struct('table', T, 'crp', crpAll, 'opt', opt);
if opt.verbose && ~isempty(T)
    fprintf('[CRP] %d significant of %d site-channel pairs (%s, alpha=%g).', ...
        sum(T.significant), height(T), opt.correct, opt.alpha);
    if nSkipped, fprintf(' %d site(s) skipped for too few trials.', nSkipped); end
    fprintf('\n');
end
end

function opt = local_dialog(EEG, opt, isCCEP)
defDir = '';
if isfield(EEG,'filepath') && ~isempty(EEG.filepath)
    defDir = fullfile(EEG.filepath, 'derivatives', 'ieeglab');
end
    function browse(src, ~)
        d = uigetdir(pwd, 'Folder for the exported results');
        if isequal(d, 0), return; end
        set(findobj(ancestor(src,'figure'), 'tag', 'exportDir'), 'string', d);
    end
ccepOn = iff(isCCEP, 'on', 'off');
uigeom = {1 [1.4 0.6] [1.4 0.6] [1.4 0.6] [1.4 0.6] [1.4 0.6] 1 1 1 1 1 1 1 [0.2 0.65 0.15]};
% One flat row of controls, in geometry order (rows of uigeom hold 1, 2 or 3).
uilist = { ...
    {'style' 'text' 'string' 'Detection and statistics' 'fontweight' 'bold'}, ...
    {'style' 'text' 'string' 'CRP response window (ms, after the artifact)'}, ...
    {'style' 'edit' 'string' num2str(opt.crp_window) 'tag' 'crpWin'}, ...
    {'style' 'text' 'string' 'N1 search window (ms)' 'enable' ccepOn}, ...
    {'style' 'edit' 'string' num2str(opt.n1_window) 'tag' 'n1Win' 'enable' ccepOn}, ...
    {'style' 'text' 'string' 'Significance level (alpha)'}, ...
    {'style' 'edit' 'string' num2str(opt.alpha) 'tag' 'alpha'}, ...
    {'style' 'text' 'string' 'Minimum trials per stimulation site'}, ...
    {'style' 'edit' 'string' num2str(opt.min_trials) 'tag' 'minTr'}, ...
    {'style' 'text' 'string' 'Correction across channels'}, ...
    {'style' 'popupmenu' 'string' {'FDR (Benjamini-Hochberg)','Bonferroni','None'} 'tag' 'corr'}, ...
    {'style' 'checkbox' 'string' 'Skip the stimulated contacts for their own site' 'value' opt.exclude_stim 'tag' 'exclStim'}, ...
    {'style' 'checkbox' 'string' 'Detect N1 responses (amplitude, latency)' 'value' opt.run_n1 'tag' 'runN1' 'enable' ccepOn}, ...
    {'style' 'checkbox' 'string' 'Canonical Response Parameterization (CRP)' 'value' opt.run_crp 'tag' 'runCRP'}, ...
    {'style' 'checkbox' 'string' 'Build the connectivity matrix (sites x contacts)' 'value' opt.run_matrix 'tag' 'runMat' 'enable' ccepOn}, ...
    {'style' 'checkbox' 'string' 'Plot figures' 'value' 1 'tag' 'doPlot'}, ...
    {'style' 'checkbox' 'string' 'Export results (TSV / JSON / MAT) to the folder below' 'value' 0 'tag' 'doExport'}, ...
    {'style' 'text' 'string' ''}, ...
    {'style' 'text' 'string' 'Folder'}, ...
    {'style' 'edit' 'string' defDir 'tag' 'exportDir' 'horizontalalignment' 'left'}, ...
    {'style' 'pushbutton' 'string' 'Browse...' 'callback' @browse} };
assert(numel(uilist) == sum(cellfun(@numel, uigeom)), 'ieeglab_stats_subject: dialog geometry and controls disagree');
[res, ~, ~, out] = inputgui('geometry', uigeom, 'uilist', uilist, ...
    'title', 'iEEGLAB - CCEP analysis', 'helpcom', 'pophelp(''ieeglab_stats_subject'');');
if isempty(res), opt = []; return; end
opt.crp_window   = str2num(out.crpWin); %#ok<ST2NM>
if isCCEP, opt.n1_window = str2num(out.n1Win); end %#ok<ST2NM>
opt.alpha        = str2double(out.alpha);
opt.min_trials   = str2double(out.minTr);
corrs = {'fdr','bonferroni','none'}; opt.correct = corrs{out.corr};
opt.exclude_stim = logical(out.exclStim);
opt.run_n1       = isCCEP && logical(out.runN1);
opt.run_crp      = logical(out.runCRP);
opt.run_matrix   = isCCEP && logical(out.runMat);
opt.plot         = logical(out.doPlot);
opt.export_dir   = '';
if logical(out.doExport), opt.export_dir = strtrim(out.exportDir); end
if numel(opt.crp_window) < 2 || opt.crp_window(2) <= opt.crp_window(1)
    error('ieeglab_stats_subject:badWindow', 'CRP window must be [start stop] ms with stop > start; got %s.', mat2str(opt.crp_window));
end
if ~isfinite(opt.alpha) || opt.alpha <= 0 || opt.alpha >= 1
    error('ieeglab_stats_subject:badAlpha', 'alpha must be between 0 and 1; got %s.', out.alpha);
end
if logical(out.doExport) && isempty(opt.export_dir)
    error('ieeglab_stats_subject:noExportDir', 'Export is ticked but no folder was given.');
end
end

function L = local_flatten(C)
L = {};
for i = 1:size(C,1)
    for j = 1:size(C,2)
        if ~isempty(C{i,j}), L{end+1} = C{i,j}; end %#ok<AGROW>
    end
end
end

function local_plot_crp(T, opt)
sig = T(T.significant,:);
figure('Color','w','Name','CRP summary','NumberTitle','off');
if isempty(sig)
    text(0.5,0.5,'No significant responses','HorizontalAlignment','center'); axis off
else
    histogram(sig.tR_ms, max(8, round(sqrt(height(sig)))));
    xlabel('Response duration \tau_R (ms)','FontWeight','bold');
    ylabel('Significant site-contact pairs','FontWeight','bold');
    title(sprintf('CRP: %d significant of %d tested (%s, \\alpha=%g)', height(sig), height(T), opt.correct, opt.alpha));
    box on
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
        q = min(1, flipud(cummin(flipud(q))));
        tmp = nan(n,1); tmp(ord) = q;
        padj(ok) = tmp;
    case 'none'
    otherwise
        warning('ieeglab_stats_subject:badCorrection', 'Unknown correction "%s"; reporting uncorrected p.', method);
end
end

function s = local_struct2str(opt)
f = fieldnames(opt); parts = {};
for i = 1:numel(f)
    v = opt.(f{i});
    if ischar(v),        parts{end+1} = sprintf('''%s'',''%s''', f{i}, v); %#ok<AGROW>
    elseif islogical(v), parts{end+1} = sprintf('''%s'',%d', f{i}, v); %#ok<AGROW>
    elseif isnumeric(v), parts{end+1} = sprintf('''%s'',%s', f{i}, mat2str(v)); %#ok<AGROW>
    end
end
s = ['struct(' strjoin(parts, ', ') ')'];
end

function out = iff(c, a, b)
if c, out = a; else, out = b; end
end
