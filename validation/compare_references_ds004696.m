function S = compare_references_ds004696(dsDir, outDir, schemes)
% compare_references_ds004696() - The CCEP analysis of ds004696 sub-02 (Ojeda
%       Valencia et al., 2023) with each re-referencing scheme of iEEGLAB, and
%       each scheme's result against the published CRP of that participant.
%
% Usage:
%   S = compare_references_ds004696('C:\data\ds004696', 'C:\data\ref_compare')
%   S = compare_references_ds004696(dsDir, outDir, {'car', 'hapwave'})
%
% dsDir    the ds004696 BIDS folder with sub-02 and derivatives/stats/sub-02 (the
%          published CRP, sub-02_ses-ieeg01_task-ccep_run-01_crp.mat)
% outDir   where the results go (about 1 GB): the preprocessed epochs
%          (preprocessed.mat, reused when present), per scheme a table of every
%          site-contact pair (pairs_<scheme>.tsv), the effect size of every pair
%          (effects_<scheme>.tsv) and the trial averages (averages_<scheme>.mat),
%          the ICA decomposition (ica_decomposition.mat, shared by both ICA
%          schemes), and summary.json. A scheme whose files exist is not run
%          again, so the script resumes after an interruption.
% schemes  which schemes to run (default: all six below)
%
% Preprocessing as Part A of the step-by-step page (the window's defaults:
% stimulation blanking, line noise at 60, 120 and 180 Hz, epochs -1 to 1 s,
% outlier trials and bad channels removed, baseline -500 to -50 ms). Then, on
% the same epochs:
%   none           the recording's own reference
%   car            common average of all good contacts
%   hapwave        the paper's reference, an adjusted common average: in each
%                  block of 64 channels of the recording (channels.tsv order),
%                  the 20% of good contacts with the lowest variance from 15 to
%                  500 ms ('varsubset', Ojeda Valencia et al., 2023)
%   carla          CARLA, iEEGLAB's default for CCEP data (Huang et al., 2024)
%   ica            ICA re-referencing, iEEGLAB's default scaling
%                  ('normalized'; Michelmann et al., 2018)
%   ica_fieldtrip  the same decomposition, scaled as the authors' code ('fieldtrip')
% and the CCEP analysis of Part A (CRP, permutation test, FDR over all pairs,
% connectivity matrix). Each scheme's decisions and response durations (tau_R)
% are compared with the published CRP, whose decision is the paper's: p at
% tau_R with the Benjamini-Yekutieli correction over all pairs of the
% participant, below 0.05.
%
% Effect size of each pair: each trial (median baseline -500 to -50 ms removed)
% is projected onto the unit-norm average of the other trials of the site over
% the paper's CRP window, 15 ms to 1 s; d is the mean of these projections over
% their standard deviation (the response in units of its trial-to-trial
% variability). Its 95% CI is from the jackknife over trials, which is
% conservative: the projections share trials, so a CI treating them as
% independent would be too narrow.
%
% References:
%   Huang, H., et al. (2024). CARLA: Adjusted common average referencing for
%     cortico-cortical evoked potential data. Journal of Neuroscience Methods,
%     407, 110153.
%   Michelmann, S., et al. (2018). Data-driven re-referencing of intracranial
%     EEG based on independent component analysis (ICA). Journal of
%     Neuroscience Methods, 307, 125-137.
%   Ojeda Valencia, G., et al. (2023). Signatures of electrical stimulation
%     driven network interactions in the human limbic system. Journal of
%     Neuroscience, 43(39), 6697-6711.
%
% Cedric Cannard, iEEGLAB, 2026

t00 = tic;
if ~isfolder(outDir), mkdir(outDir); end
allSchemes = {'none', 'car', 'hapwave', 'carla', 'ica', 'ica_fieldtrip'};
if nargin < 3 || isempty(schemes), schemes = allSchemes; end
schemes = cellstr(schemes);
if ~all(ismember(schemes, allSchemes))
    error('compare_references:scheme', 'Schemes are %s.', strjoin(allSchemes, ', '));
end
statsOpt = struct('run_n1', false, 'run_crp', true, 'run_matrix', true, 'matrix_source', 'crp', ...
    'plot', false, 'verbose', false);
tKeep = [-200 800];                                         % ms kept in the averages
effWin = [15 1000]; effBase = [-500 -50];                   % ms, effect sizes

% ---------------- preprocessed epochs (once) ----------------
fPre = fullfile(outDir, 'preprocessed.mat');
rec = fullfile(dsDir, 'sub-02', 'ses-ieeg01', 'ieeg', 'sub-02_ses-ieeg01_task-ccep_run-01_ieeg.mefd');
if isfile(fPre)
    L = load(fPre, 'EEG'); EEG0 = L.EEG; clear L
    fprintf('Preprocessed epochs loaded: %d contacts x %d epochs\n', EEG0.nbchan, EEG0.trials);
    if ~isfield(EEG0.etc.ieeglab_import, 'recording_channels')    % saved before the import recorded it
        I = ieeglab_import(rec, 'info', true);
        EEG0.etc.ieeglab_import.recording_channels = I.labels;
    end
else
    t0 = tic;
    [~, EEG0] = evalc('ieeglab_import(rec)');
    [~, EEG0] = evalc('ieeglab_load(EEG0, struct(''event_field'', ''electrical_stimulation_site''))');
    [~, EEG0] = evalc(['ieeglab_preprocess(EEG0, struct(''apply_blank'', true, ''apply_notch'', true, ' ...
        '''notch'', [60 120 180], ''apply_epoch'', true, ''epoch_window'', [-1000 1000], ''reject_trials'', true, ' ...
        '''remove_rare_cond'', true, ''min_trials'', 5, ''remove_no_coords'', true, ''remove_bad_channels'', true, ' ...
        '''apply_baseline'', true, ''baseline_period'', [-500 -50], ''plot'', false))']);
    EEG = EEG0; save(fPre, 'EEG', '-v7.3'); clear EEG
    fprintf('Preprocessed: %d contacts x %d epochs (%.0f s)\n', EEG0.nbchan, EEG0.trials, toc(t0));
end
[sites, stimIdx] = ieeglab_epoch_sites(EEG0);
uSites = unique(sites, 'stable');
labels = {EEG0.chanlocs.labels};
tk = EEG0.times >= tKeep(1) & EEG0.times <= tKeep(2);
P = load(fullfile(dsDir, 'derivatives', 'stats', 'sub-02', 'sub-02_ses-ieeg01_task-ccep_run-01_crp.mat'));
pub = local_published(P);
fprintf('Published CRP: %d pairs, %d respond (Benjamini-Yekutieli over all pairs)\n', numel(pub.key), nnz(pub.significant));

% ---------------- each scheme ----------------
fDec = fullfile(outDir, 'ica_decomposition.mat');
for s = 1:numel(schemes)
    sc = schemes{s};
    fT = fullfile(outDir, ['pairs_' sc '.tsv']); fA = fullfile(outDir, ['averages_' sc '.mat']);
    fE = fullfile(outDir, ['effects_' sc '.tsv']);
    havePairs = isfile(fT) && isfile(fA);
    if havePairs && isfile(fE)
        T = local_attach_published(readtable(fT, 'FileType', 'text', 'Delimiter', '\t'), pub);
        writetable(T, fT, 'FileType', 'text', 'Delimiter', '\t');
        fprintf('%s: done before; published decisions updated\n', sc); continue
    end
    t0 = tic; nBroad = NaN;
    switch sc
        case 'none',    EEG = EEG0;
        case 'car',     [~, EEG] = evalc('pop_ieeglab_reref(EEG0, ''method'', ''car'')');
        case 'hapwave', [~, EEG] = evalc('pop_ieeglab_reref(EEG0, ''method'', ''varsubset'', ''fraction'', 0.2)');
        case 'carla',   [~, EEG] = evalc('pop_ieeglab_reref(EEG0, ''method'', ''carla'')');
        case {'ica', 'ica_fieldtrip'}
            o = struct('weights', local_if(strcmp(sc, 'ica'), 'normalized', 'fieldtrip'));
            if isfile(fDec), L = load(fDec, 'info'); o.decomposition = L.info; clear L; end
            [~, EEG, info] = evalc('ieeglab_icaref(EEG0, o)');
            if ~isfile(fDec), save(fDec, 'info', '-v7'); end
            nBroad = nnz(info.broad);
    end
    tRef = toc(t0);
    [tf, loc] = ismember({EEG.chanlocs.labels}, labels);
    if ~all(tf), error('compare_references:channels', '%s changed the channels.', sc); end

    % effect size of every pair
    t1 = tic;
    Ef = local_effects(EEG, loc, sites, uSites, stimIdx, labels, effWin, effBase);
    Ef.scheme = repmat(string(sc), height(Ef), 1);
    writetable(Ef, fE, 'FileType', 'text', 'Delimiter', '\t');
    tEff = toc(t1);

    tC = NaN;
    if ~havePairs
        % trial averages per site (sites x contacts x time), the stimulated contacts NaN;
        % and what the scheme removed: the trial average of before minus after
        A = nan(numel(uSites), numel(labels), nnz(tk), 'single'); R = A;
        for g = 1:numel(uSites)
            tr = find(sites == uSites(g));
            a = mean(double(EEG.data(:, tk, tr)), 3);
            r = mean(double(EEG0.data(loc, tk, tr)), 3) - a;
            stim = unique(vertcat(stimIdx{tr}));
            a(ismember(loc, stim), :) = NaN; r(ismember(loc, stim), :) = NaN;
            A(g, loc, :) = single(a); R(g, loc, :) = single(r);
        end
        times = EEG0.times(tk); site_names = cellstr(uSites); %#ok<NASGU>
        save(fA, 'A', 'R', 'times', 'site_names', 'labels', '-v7');

        % CCEP analysis
        t2 = tic;
        [~, E] = evalc('ieeglab_stats_subject(EEG, statsOpt)');
        T = E.ieeglab.stats.table; clear E
        T = T(:, intersect({'site', 'channel', 'tR_ms', 'explained_var', 'p', 'p_adj', 'significant', 'n_trials'}, ...
            T.Properties.VariableNames, 'stable'));
        T.scheme = repmat(string(sc), height(T), 1);
        tC = toc(t2);
    else
        T = readtable(fT, 'FileType', 'text', 'Delimiter', '\t');
    end
    clear EEG
    T = local_attach_published(T, pub);
    writetable(T, fT, 'FileType', 'text', 'Delimiter', '\t');
    m = T.pub_tested == 1 & ~isnan(T.significant);
    fprintf(['%s: re-referencing %.0f s, effect sizes %.0f s%s; %d of %d pairs respond; published pairs %d ' ...
        '(%d respond): same decision %.1f%%; median d %.2f%s\n'], sc, tRef, tEff, ...
        local_if(~isnan(tC), sprintf(', CCEP analysis %.0f s', tC), ''), nnz(T.significant == 1), ...
        nnz(~isnan(T.significant)), nnz(m), nnz(T.pub_significant(m) == 1), ...
        100 * mean((T.significant(m) == 1) == (T.pub_significant(m) == 1)), median(Ef.d, 'omitnan'), ...
        local_if(~isnan(nBroad), sprintf('; %d broad components removed', nBroad), ''));
end

% ---------------- summary of the schemes run so far ----------------
done = allSchemes(cellfun(@(x) isfile(fullfile(outDir, ['pairs_' x '.tsv'])), allSchemes));
S = struct('scheme', done, 'pairs_tested', NaN, 'pairs_significant', NaN, 'published_pairs', NaN, ...
    'same_decision_pct', NaN, 'published_significant', NaN, 'significant_both', NaN, 'significant_only_ieeglab', NaN, ...
    'significant_only_published', NaN, 'tR_median_abs_diff_ms', NaN, 'tR_within_50ms_pct', NaN, ...
    'removed_rms_uv_10_300ms', NaN, 'response_rms_uv_10_300ms', NaN);
for s = 1:numel(done)
    T = readtable(fullfile(outDir, ['pairs_' done{s} '.tsv']), 'FileType', 'text', 'Delimiter', '\t');
    m = T.pub_tested == 1 & ~isnan(T.significant);
    is = T.significant(m) == 1; ps = T.pub_significant(m) == 1;
    b = m & T.significant == 1 & T.pub_significant == 1;
    d = abs(T.tR_ms(b) - T.pub_tR_ms(b));
    S(s).pairs_tested = nnz(~isnan(T.significant)); S(s).pairs_significant = nnz(T.significant == 1);
    S(s).published_pairs = nnz(m); S(s).same_decision_pct = 100 * mean(is == ps);
    S(s).published_significant = nnz(ps); S(s).significant_both = nnz(is & ps);
    S(s).significant_only_ieeglab = nnz(is & ~ps); S(s).significant_only_published = nnz(~is & ps);
    S(s).tR_median_abs_diff_ms = median(d); S(s).tR_within_50ms_pct = 100 * mean(d <= 50);
    L = load(fullfile(outDir, ['averages_' done{s} '.mat']), 'A', 'R', 'times');
    w = L.times >= 10 & L.times <= 300;
    S(s).removed_rms_uv_10_300ms = median(sqrt(mean(double(L.R(:, :, w)).^2, 3)), 'all', 'omitnan');
    S(s).response_rms_uv_10_300ms = median(sqrt(mean(double(L.A(:, :, w)).^2, 3)), 'all', 'omitnan');
end
fid = fopen(fullfile(outDir, 'summary.json'), 'w'); fprintf(fid, '%s', jsonencode(S, 'PrettyPrint', true)); fclose(fid);
disp(struct2table(S));
fprintf('Done in %.0f min.\n', toc(t00) / 60);
end

function Ef = local_effects(EEG, loc, sites, uSites, stimIdx, labels, win, base)
% Per site-contact pair: d of the leave-one-out projections over win (ms), after a
% per-trial median baseline over base (ms); jackknife 95% CI over the trials.
tw = EEG.times >= win(1) & EEG.times <= win(2);
tb = EEG.times >= base(1) & EEG.times <= base(2);
nC = EEG.nbchan; nS = numel(uSites); n = nS * nC;
[site, chan] = deal(strings(n, 1));
[nTr, d, se, lo, hi, pm, ps] = deal(nan(n, 1));
for g = 1:nS
    tr = find(sites == uSites(g)); K = numel(tr);
    X = double(EEG.data(:, :, tr));
    X = X(:, tw, :) - median(X(:, tb, :), 2, 'omitnan');          % contacts x time x trials
    [dg, xm, xs] = local_loo_d(X);
    dj = zeros(nC, K);
    for i = 1:K, dj(:, i) = local_loo_d(X(:, :, [1:i-1, i+1:K])); end
    sg = sqrt((K - 1) / K * sum((dj - mean(dj, 2)).^2, 2));
    q = sqrt((K - 1) * (1 / betaincinv(0.05, (K - 1) / 2, 0.5) - 1));   % t(0.975, K - 1)
    stim = ismember(loc(:), unique(vertcat(stimIdx{tr})));
    dg(stim) = NaN;
    k = (g - 1) * nC + (1:nC)';
    site(k) = uSites(g); chan(k) = string(labels(loc)); nTr(k) = K;
    d(k) = dg; se(k) = sg; lo(k) = dg - q * sg; hi(k) = dg + q * sg; pm(k) = xm; ps(k) = xs;
end
Ef = table(site, chan, nTr, d, se, lo, hi, pm, ps, 'VariableNames', ...
    {'site', 'channel', 'n_trials', 'd', 'd_se', 'd_ci_low', 'd_ci_high', 'proj_uv', 'proj_sd_uv'});
Ef = Ef(~isnan(Ef.d), :);                                          % not the stimulated contacts
end

function [d, xm, xs] = local_loo_d(X)
% X: contacts x time x trials. Each trial projected onto the unit-norm mean of the others.
K = size(X, 3);
M = (sum(X, 3) - X) / (K - 1);
x = reshape(sum(X .* M, 2) ./ sqrt(sum(M .^ 2, 2)), size(X, 1), K);
xm = mean(x, 2); xs = std(x, 0, 2); d = xm ./ xs;
end

function pub = local_published(P)
% Every pair of the published CRP: key "site|contact", p at tau_R, its
% Benjamini-Yekutieli adjustment over all pairs of the participant, the paper's
% decision (adjusted p < 0.05), tau_R in ms
key = strings(0, 1); p = []; tR = [];
for s = 1:numel(P.average_ccep_names)
    cs = ieeglab_canonical_site(P.average_ccep_names(s), P.channel_names);
    for ch = 1:size(P.crp_out, 1)
        Q = P.crp_out(ch, s);
        if isempty(Q.crp_parms) || ~isfield(Q.crp_projs, 'p_value_tR'), continue; end
        key(end + 1, 1) = cs(1) + "|" + string(P.channel_names{ch}); %#ok<AGROW>
        p(end + 1, 1) = Q.crp_projs.p_value_tR; %#ok<AGROW>
        tR(end + 1, 1) = Q.crp_parms.tR * 1000; %#ok<AGROW>
    end
end
m = numel(p); [ps, o] = sort(p);
q = ps * m * sum(1 ./ (1:m)) ./ (1:m)';
q = min(1, flipud(cummin(flipud(q))));
padj = zeros(m, 1); padj(o) = q;
pub = struct('key', key, 'p', p, 'p_adj', padj, 'significant', double(padj < 0.05), 'tR_ms', tR);
end

function T = local_attach_published(T, pub)
key = string(T.site) + "|" + string(T.channel);
[in, i] = ismember(key, pub.key);
T.pub_tested = in;
for f = {'p', 'p_adj', 'significant', 'tR_ms'}
    v = nan(height(T), 1); v(in) = pub.(f{1})(i(in));
    T.(['pub_' f{1}]) = v;
end
end

function s = local_if(c, a, b)
if c, s = a; else, s = b; end
end
