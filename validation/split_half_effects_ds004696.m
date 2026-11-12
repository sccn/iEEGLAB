function S = split_half_effects_ds004696(dsDir, outDir, schemes)
% split_half_effects_ds004696() - Effect sizes of the re-referencing schemes
%       of compare_references_ds004696 on ds004696 sub-02, without selection
%       bias: pairs are selected on one half of the trials and measured on the
%       other.
%
% Usage:
%   S = split_half_effects_ds004696('C:\data\ds004696', 'C:\data\ref_compare_sub02')
%   S = split_half_effects_ds004696(dsDir, outDir, {'hapwave', 'carla'})
%
% dsDir    the ds004696 BIDS folder (for the recording's channel list)
% outDir  the folder of compare_references_ds004696, with preprocessed.mat and
%          ica_decomposition.mat (run that function first). Writes
%          split_half_pairs.tsv (every pair and fold) and split_half_summary.json.
% schemes  which schemes to measure (default: all six; 'hapwave' is always run,
%          as it selects the pairs)
%
% Comparing schemes on the pairs that respond in the paper favours the paper's
% reference: those pairs were chosen with its decisions, so trials where its
% noise happened to be low are over-represented among them, and its effect size
% on the pairs without a response is biased low. Here the trials of each
% stimulation site are split into odd and even (recording order). On one half,
% the pairs are selected with the published reference (the lowest-variance
% subset, 'varsubset') and the paper's test: CRP's t-test at the response
% duration, Benjamini-Yekutieli over all pairs, p < 0.05 (Ojeda Valencia et al.,
% 2023). On the other half, every scheme's effect size d is measured as in
% compare_references_ds004696 (each trial, median baseline -500 to -50 ms
% removed, projected onto the unit-norm mean of the other trials from 15 ms to
% 1 s; d = mean / SD). The halves are then swapped. Sites with fewer than 8
% trials are left out, so each half has at least 4. The schemes re-reference all
% trials of a site at once, as in the analysis itself: CARLA and the published
% reference choose their reference contacts from both halves.
%
% Summary: mean d on the selected responding and non-responding pairs, and each
% scheme's difference from the published reference, pooled over both folds,
% with 95% cluster-robust t-intervals over the stimulation sites (CR1; Cameron &
% Miller, 2015).
%
% References:
%   Cameron, A. C., & Miller, D. L. (2015). A practitioner's guide to
%     cluster-robust inference. Journal of Human Resources, 50(2), 317-372.
%   Ojeda Valencia, G., et al. (2023). Signatures of electrical stimulation
%     driven network interactions in the human limbic system. Journal of
%     Neuroscience, 43(39), 6697-6711.
%
% Cedric Cannard, iEEGLAB, 2026

t00 = tic;
allSchemes = {'none', 'car', 'hapwave', 'carla', 'ica', 'ica_fieldtrip'};
if nargin < 3 || isempty(schemes), schemes = allSchemes; end
schemes = allSchemes(ismember(allSchemes, union(cellstr(schemes), {'hapwave'})));
effWin = [15 1000]; effBase = [-500 -50]; minTrials = 8;
selOpt = struct('run_n1', false, 'run_crp', true, 'run_matrix', false, 'crp_test', 'ttest', ...
    'correct', 'by', 'correct_scope', 'subject', 'min_trials', 4, 'plot', false, 'verbose', false);

L = load(fullfile(outDir, 'preprocessed.mat'), 'EEG'); EEG0 = L.EEG; clear L
fprintf('Preprocessed epochs loaded: %d contacts x %d epochs\n', EEG0.nbchan, EEG0.trials);
if ~isfield(EEG0.etc.ieeglab_import, 'recording_channels')    % saved before the import recorded it
    rec = fullfile(dsDir, 'sub-02', 'ses-ieeg01', 'ieeg', 'sub-02_ses-ieeg01_task-ccep_run-01_ieeg.mefd');
    I = ieeglab_import(rec, 'info', true);
    EEG0.etc.ieeglab_import.recording_channels = I.labels;   % the published reference's blocks
end
[sites, stimIdx] = ieeglab_epoch_sites(EEG0);
labels = {EEG0.chanlocs.labels};
uSites = unique(sites, 'stable');
nTr = arrayfun(@(s) nnz(sites == s), uSites);
uSites = uSites(nTr >= minTrials);
fold = zeros(1, EEG0.trials);                               % 1 odd, 2 even, 0 left out
for g = 1:numel(uSites)
    tr = find(sites == uSites(g));
    fold(tr(1:2:end)) = 1; fold(tr(2:2:end)) = 2;
end
fprintf('%d of %d sites with at least %d trials; %d and %d epochs in the two halves.\n', numel(uSites), ...
    numel(nTr), minTrials, nnz(fold == 1), nnz(fold == 2));

% ---------------- each scheme: d on each half ----------------
D = struct();
fDec = fullfile(outDir, 'ica_decomposition.mat');
for s = 1:numel(schemes)
    sc = schemes{s}; t0 = tic;
    switch sc
        case 'none',    EEG = EEG0;
        case 'car',     [~, EEG] = evalc('pop_ieeglab_reref(EEG0, ''method'', ''car'')');
        case 'hapwave', [~, EEG] = evalc('pop_ieeglab_reref(EEG0, ''method'', ''varsubset'', ''fraction'', 0.2)');
        case 'carla',   [~, EEG] = evalc('pop_ieeglab_reref(EEG0, ''method'', ''carla'')');
        case {'ica', 'ica_fieldtrip'}
            Ld = load(fDec, 'info');
            o = struct('weights', local_if(strcmp(sc, 'ica'), 'normalized', 'fieldtrip'), 'decomposition', Ld.info);
            [~, EEG] = evalc('ieeglab_icaref(EEG0, o)');
    end
    if ~isequal({EEG.chanlocs.labels}, labels), error('split_half:channels', '%s changed the channels.', sc); end
    for f = 1:2
        D.(sc){f} = local_d(EEG, sites, stimIdx, uSites, fold == f, effWin, effBase);
    end
    if strcmp(sc, 'hapwave')                                % selection on each half
        for f = 1:2
            [~, Ef] = evalc('pop_select(EEG, ''trial'', find(fold == f))');
            [~, Ef] = evalc('ieeglab_stats_subject(Ef, selOpt)');
            Sel{f} = Ef.ieeglab.stats.table; clear Ef %#ok<AGROW>
        end
    end
    clear EEG
    fprintf('%s: re-referencing and effect sizes %.0f s\n', sc, toc(t0));
end

% ---------------- pairs: selected on one half, measured on the other ----------------
rows = cell(2, 1);
for f = 1:2
    T = Sel{f};
    T = T(:, {'site', 'channel', 'p', 'p_adj', 'significant', 'n_trials'});
    T.fold = repmat(f, height(T), 1);
    key = string(T.site) + "|" + string(T.channel);
    m = 3 - f;                                              % measured on the other half
    for s = 1:numel(schemes)
        Dm = D.(schemes{s}){m};
        [in, i] = ismember(key, Dm.key);
        v = nan(height(T), 1); v(in) = Dm.d(i(in));
        T.(['d_' schemes{s}]) = v;
    end
    rows{f} = T;
end
P = [rows{1}; rows{2}];
dcols = strcat('d_', schemes);
ok = ~isnan(P.significant) & all(isfinite(P{:, dcols}), 2);
P = P(ok, :);
writetable(P, fullfile(outDir, 'split_half_pairs.tsv'), 'FileType', 'text', 'Delimiter', '\t');

% ---------------- summary ----------------
resp = P.significant == 1; g = string(P.site);
S = struct('scheme', schemes, 'resp', NaN, 'resp_ci', [], 'non', NaN, 'non_ci', [], ...
    'diff', NaN, 'diff_ci', [], 'diff_non', NaN, 'diff_non_ci', []);
for s = 1:numel(schemes)
    x = P.(dcols{s}); x0 = P.d_hapwave;
    [S(s).resp, S(s).resp_ci] = local_crt(x(resp), g(resp));
    [S(s).non, S(s).non_ci] = local_crt(x(~resp), g(~resp));
    [S(s).diff, S(s).diff_ci] = local_crt(x(resp) - x0(resp), g(resp));
    [S(s).diff_non, S(s).diff_non_ci] = local_crt(x(~resp) - x0(~resp), g(~resp));
end
info = struct('n_sites', numel(uSites), 'min_trials', minTrials, 'n_pairs_fold', [nnz(P.fold == 1), nnz(P.fold == 2)], ...
    'n_resp_fold', [nnz(resp & P.fold == 1), nnz(resp & P.fold == 2)], 'n_sites_resp', numel(unique(g(resp))), ...
    'n_sites_non', numel(unique(g(~resp))), 'schemes', S);
fid = fopen(fullfile(outDir, 'split_half_summary.json'), 'w');
fprintf(fid, '%s', jsonencode(info, 'PrettyPrint', true)); fclose(fid);
fprintf('Pairs per half %d and %d; responding %d and %d (published reference, paper''s test).\n', ...
    info.n_pairs_fold, info.n_resp_fold);
for s = 1:numel(S)
    fprintf(['%-14s responding %5.2f [%5.2f, %5.2f]  non %5.2f [%5.2f, %5.2f]  minus published: %+5.2f ' ...
        '[%+5.2f, %+5.2f]  non %+5.2f [%+5.2f, %+5.2f]\n'], S(s).scheme, S(s).resp, S(s).resp_ci, S(s).non, ...
        S(s).non_ci, S(s).diff, S(s).diff_ci, S(s).diff_non, S(s).diff_non_ci);
end
fprintf('Done in %.0f min.\n', toc(t00) / 60);
end

function Dm = local_d(EEG, sites, stimIdx, uSites, use, win, base)
% d of every site-contact pair on the epochs in use (see compare_references_ds004696)
tw = EEG.times >= win(1) & EEG.times <= win(2);
tb = EEG.times >= base(1) & EEG.times <= base(2);
labels = string({EEG.chanlocs.labels});
key = strings(0, 1); d = [];
for g = 1:numel(uSites)
    tr = find(sites == uSites(g) & use);
    X = double(EEG.data(:, :, tr));
    X = X(:, tw, :) - median(X(:, tb, :), 2, 'omitnan');           % contacts x time x trials
    K = numel(tr);
    M = (sum(X, 3) - X) / (K - 1);                                 % mean of the other trials
    x = reshape(sum(X .* M, 2) ./ sqrt(sum(M .^ 2, 2)), size(X, 1), K);
    dg = mean(x, 2) ./ std(x, 0, 2);
    keep = ~ismember((1:EEG.nbchan)', unique(vertcat(stimIdx{tr})));
    key = [key; uSites(g) + "|" + labels(keep)']; %#ok<AGROW>
    d = [d; dg(keep)]; %#ok<AGROW>
end
Dm = struct('key', key, 'd', d);
end

function [m, ci] = local_crt(x, g)
% mean with its 95% cluster-robust (CR1) t-interval, clusters g, G - 1 df
[~, ~, gi] = unique(g); G = max(gi); m = mean(x);
se = sqrt(G / (G - 1) * sum(accumarray(gi(:), x(:) - m) .^ 2)) / numel(x);
q = sqrt((G - 1) * (1 / betaincinv(0.05, (G - 1) / 2, 0.5) - 1));   % t(0.975, G - 1)
ci = [m - q * se, m + q * se];
end

function s = local_if(c, a, b)
if c, s = a; else, s = b; end
end
