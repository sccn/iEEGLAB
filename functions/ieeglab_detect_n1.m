function [EEG, T, com] = ieeglab_detect_n1(EEG, opt)
% ieeglab_detect_n1() - Detect N1 evoked responses in CCEP data.
%
% For every stimulation site x recording contact, tests whether the trial-
% averaged response contains an early deflection larger than chance, and
% reports its amplitude and latency. The most widely reported CCEP measure;
% complements CRP (ieeglab_stats_subject), which describes the response shape.
%
% Usage:
%   EEG         = ieeglab_detect_n1(EEG)
%   [EEG, T]    = ieeglab_detect_n1(EEG, opt)
%
% Options:
%   .n1_window     [1x2] ms, where to look for the peak. Default [10 100].
%                  Must start after the stimulation artifact.
%   .baseline      [1x2] ms. Defines both the reference level (its median is
%                  subtracted) and the noise level (its SD). Default [-500 -10].
%   .method        'permutation' (default) - sign-flip test on the maximum of
%                                 the whole window: flipping the sign of random
%                                 trials cancels a consistent evoked response
%                                 but not noise, and taking the same maximum in
%                                 each permutation accounts for searching the
%                                 window. Exact (all sign patterns) for up to
%                                 9 trials.
%                  'sd'           - the van Blooijs et al. (2018) detection RULE:
%                                 a peak above threshold x baseline SD. A rule,
%                                 not a test: p is NaN and no correction applies.
%   .n_perm        permutations. Default 1000.
%   .threshold     SD multiple for method 'sd'. Default 3.4.
%   .alpha         significance level. Default 0.05.
%   .correct       'fdr' (default), 'bonferroni' or 'none', across contacts
%                  within each stimulation site (permutation method).
%   .min_trials    minimum trials per site (at least 2). Default 5.
%   .exclude_stim  skip each site's stimulated contacts. Default true.
%   .require_peak  the maximum must be a local extremum inside the window, not
%                  its first or last sample (a slope running out of the window
%                  is not an N1). Default true.
%   .polarity      'abs' (default) counts the largest deflection of either sign;
%                  'negative' counts only negative-going peaks, the convention
%                  used by erdetect; 'positive' only positive-going ones. On
%                  HAPwave sub-02, 37% of responsive pairs have their largest
%                  early deflection positive, so the two conventions disagree
%                  for about a third of the data.
%   .min_baseline_sd  floor on the baseline SD in uV, so that a quiet contact
%                  cannot produce a large z. Default 50, which with the default
%                  threshold of 3.4 gives erdetect's effective 170 uV criterion.
%                  Set 0 to disable.
%
% erdetect-equivalent settings:
%   opt = struct('method','sd', 'threshold',3.4, 'min_baseline_sd',50, ...
%                'polarity','negative', 'n1_window',[9 90], 'baseline',[-1000 -100]);
%   A baseline that does not fit the epoch is clamped to the available
%   pre-stimulus period with a warning, rather than silently returning no
%   detections.
%   .verbose       default true
%
% Output T (also EEG.ieeglab.n1.table): site, channel, n1_amplitude_uv,
% n1_latency_ms, baseline_sd, z, p, p_adj, significant, n_trials.
%
% Amplitude is measured from the baseline median. Polarity depends on the
% reference, so the largest ABSOLUTE deflection is taken and its sign kept.
%
% Recomputing N1 invalidates a connectivity matrix built from N1.
%
% References:
%   van Blooijs D, et al. (2018). Evoked directional network characteristics of
%   epileptogenic tissue derived from single pulse electrical stimulation.
%   Human Brain Mapping 39(11):4611-4622.
%   ER-detect (MultimodalNeuroimagingLab), the reference implementation of the
%   SD rule used here: https://github.com/MultimodalNeuroimagingLab/erdetect
%   Maris E, Oostenveld R (2007). Nonparametric statistical testing of EEG- and
%   MEG-data. J Neurosci Methods 164(1):177-190.
%
% Cedric Cannard, iEEGLAB, 2026

com = '';
if nargin < 2 || isempty(opt), opt = struct(); end
def = struct('n1_window',[10 100], 'baseline',[-500 -10], 'method','permutation', ...
             'n_perm',1000, 'threshold',3.4, 'alpha',0.05, 'correct','fdr', ...
             'min_trials',5, 'exclude_stim',true, 'require_peak',true, ...
             'polarity','abs', 'min_baseline_sd',50, 'verbose',true);
fn = fieldnames(def);
for i = 1:numel(fn)
    if ~isfield(opt,fn{i}) || isempty(opt.(fn{i})), opt.(fn{i}) = def.(fn{i}); end
end
opt.method  = lower(char(opt.method));
opt.correct = lower(char(opt.correct));
if ~any(strcmp(opt.method, {'permutation','sd'}))
    error('ieeglab_detect_n1:badMethod', 'method must be ''permutation'' or ''sd''; got ''%s''.', opt.method);
end
opt.polarity = lower(char(opt.polarity));
if ~any(strcmp(opt.polarity, {'abs','negative','positive'}))
    error('ieeglab_detect_n1:badPolarity', ...
        'polarity must be ''abs'', ''negative'' or ''positive''; got ''%s''.', opt.polarity);
end
if ~isnumeric(opt.min_baseline_sd) || ~isscalar(opt.min_baseline_sd) || opt.min_baseline_sd < 0
    error('ieeglab_detect_n1:badMinSd', 'min_baseline_sd must be a non-negative scalar (uV).');
end
if opt.min_trials < 2
    warning('ieeglab_detect_n1:minTrials', 'min_trials raised to 2: a site needs at least two trials.');
    opt.min_trials = 2;
end

if EEG.trials < 2
    error('ieeglab_detect_n1:notEpoched', ...
        'N1 detection needs epoched data with at least 2 trials (this dataset has %d).', EEG.trials);
end
dmode = ieeglab_detect_mode(EEG);
if ~strcmp(dmode,'ccep')
    warning('ieeglab_detect_n1:notCCEP', ...
        ['N1 detection is defined for CCEP data; this dataset looks like "%s". ' ...
         'Results are grouped by event type but the "stimulation site" reading does not apply.'], dmode);
end

t  = double(EEG.times(:))';
iN = t >= opt.n1_window(1) & t <= opt.n1_window(2);
iB = t >= opt.baseline(1)  & t <= opt.baseline(2);
if nnz(iN) < 3
    error('ieeglab_detect_n1:badWindow', ...
        'Search window [%g %g] ms selects %d samples of an epoch spanning [%g %g] ms.', ...
        opt.n1_window(1), opt.n1_window(2), nnz(iN), t(1), t(end));
end
if nnz(iB) < 10
    % The requested baseline does not fit this epoch. Clamp it to the pre-stimulus
    % period that exists rather than returning nothing, which is the failure mode
    % erdetect has with its default baseline of -1000 to -100 ms on short epochs.
    lo = max(opt.baseline(1), t(1));
    hi = min(opt.baseline(2), -1000/EEG.srate);          % one sample before zero
    iB = t >= lo & t <= hi;
    if nnz(iB) < 10
        error('ieeglab_detect_n1:badBaseline', ...
            ['Baseline [%g %g] ms selects only %d samples of an epoch spanning [%g %g] ms; ' ...
             'need at least 10 for a noise estimate.'], ...
            opt.baseline(1), opt.baseline(2), nnz(iB), t(1), t(end));
    end
    warning('ieeglab_detect_n1:baselineClamped', ...
        'Baseline [%g %g] ms does not fit this epoch; using [%g %g] ms (%d samples).', ...
        opt.baseline(1), opt.baseline(2), lo, hi, nnz(iB));
    opt.baseline = [lo hi];
end
use = iB | iN;                  % only these samples are ever needed
iBu = iB(use); iNu = iN(use); tN = t(iN);

labels = string({EEG.chanlocs.labels});
[sites, stimIdxEp] = ieeglab_epoch_sites(EEG);
[uSites, ~, grp] = unique(sites);
isBadCh = false(1, EEG.nbchan);
if isfield(EEG.chanlocs,'status')
    isBadCh = cellfun(@(x) ~isempty(x) && strcmpi(char(x),'bad'), {EEG.chanlocs.status});
end

rows = {};
for g = 1:numel(uSites)
    tr = find(grp == g);
    if numel(tr) < opt.min_trials
        if opt.verbose
            fprintf('[N1] site %-14s skipped (%d trials < %d)\n', char(uSites(g)), numel(tr), opt.min_trials);
        end
        continue
    end
    stimIdx = [];
    if opt.exclude_stim, stimIdx = unique(vertcat(stimIdxEp{tr})); end

    pvals = nan(EEG.nbchan,1);
    tmp   = cell(EEG.nbchan,1);
    for ch = 1:EEG.nbchan
        if ismember(ch, stimIdx) || isBadCh(ch), continue; end
        X = reshape(double(EEG.data(ch, use, tr)), nnz(use), numel(tr));   % samples x trials
        X = X(:, all(isfinite(X), 1));                                      % trials with NaN left out
        K = size(X, 2);
        if K < opt.min_trials, continue; end

        if strcmp(opt.method, 'permutation')
            [S, exact] = local_signs(K, opt.n_perm);
        else
            S = ones(K, 1); exact = false;
        end
        [stat, amp, lat, sd] = local_stat(X * S / K, iBu, iNu, tN, opt.require_peak, ...
                                          opt.polarity, opt.min_baseline_sd);

        switch opt.method
            case 'permutation'
                if exact
                    p = mean(stat >= stat(1));                     % all patterns, observed included
                else
                    p = (1 + sum(stat(2:end) >= stat(1))) / numel(stat);
                end
            otherwise
                p = NaN;                                           % a rule, not a test
        end
        pvals(ch) = p;
        tmp{ch} = struct('amp',amp, 'lat',lat, 'sd',sd, 'z',stat(1), 'K',K);
    end

    if strcmp(opt.method, 'permutation')
        padj = local_correct(pvals, opt.correct);
    else
        padj = nan(size(pvals));
    end
    for ch = 1:EEG.nbchan
        if isempty(tmp{ch}), continue; end
        r = tmp{ch};
        if strcmp(opt.method, 'permutation')
            sig = padj(ch) < opt.alpha && isfinite(r.amp);
        else
            sig = r.z >= opt.threshold && isfinite(r.amp);
        end
        rows(end+1,:) = { char(uSites(g)), char(labels(ch)), r.amp, r.lat, ...
            r.sd, r.z, pvals(ch), padj(ch), sig, r.K }; %#ok<AGROW>
    end
end

% Results derived from an older N1 table are no longer valid
if isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'ccep_matrix') && ~isempty(EEG.ieeglab.ccep_matrix) ...
        && isfield(EEG.ieeglab.ccep_matrix,'source') && strcmp(EEG.ieeglab.ccep_matrix.source, 'n1')
    EEG.ieeglab = rmfield(EEG.ieeglab, 'ccep_matrix');
    if opt.verbose, fprintf('[N1] The connectivity matrix built from the previous N1 results was removed.\n'); end
end

com = sprintf('EEG = ieeglab_detect_n1(EEG, %s);', ieeglab_literal(opt));
if isempty(rows)
    warning('ieeglab_detect_n1:noResults', ...
        'No N1 fits succeeded. Most often too few trials per site (minimum %d).', opt.min_trials);
    T = table();
    EEG.ieeglab.n1 = struct('table', T, 'opt', opt);
    return
end

T = cell2table(rows, 'VariableNames', ...
    {'site','channel','n1_amplitude_uv','n1_latency_ms','baseline_sd','z','p','p_adj','significant','n_trials'});
T = sortrows(T, {'site','channel'});
EEG.ieeglab.n1 = struct('table', T, 'opt', opt);

if opt.verbose
    if strcmp(opt.method,'permutation')
        fprintf('\n[N1] %d responses of %d site-contact pairs (sign-flip max test, %d permutations, %s-corrected, alpha=%g).\n', ...
            sum(T.significant), height(T), opt.n_perm, opt.correct, opt.alpha);
    else
        fprintf('\n[N1] %d responses of %d site-contact pairs (rule: peak > %g x baseline SD; not a significance test).\n', ...
            sum(T.significant), height(T), opt.threshold);
    end
    if any(T.significant)
        s = T(T.significant,:);
        fprintf('[N1] Latency %.0f-%.0f ms (median %.0f), |amplitude| %.0f-%.0f uV.\n', ...
            min(s.n1_latency_ms), max(s.n1_latency_ms), median(s.n1_latency_ms), ...
            min(abs(s.n1_amplitude_uv)), max(abs(s.n1_amplitude_uv)));
    end
end

% Plausibility: a CCEP N1 is tens to a few hundred microvolts. Thousands mean
% residual stimulation artifact rather than a neural response.
BIG = 1000;
nBig = sum(T.significant & abs(T.n1_amplitude_uv) > BIG);
if nBig > 0
    warning('ieeglab_detect_n1:implausibleAmplitude', ...
        ['%d detected response(s) exceed %g uV, far above a physiological N1. They are very ' ...
         'likely residual stimulation artifact. Check that the artifact was blanked before ' ...
         'filtering and downsampling, and that n1_window starts late enough (now %g ms).'], ...
        nBig, BIG, opt.n1_window(1));
end
if EEG.srate < 500
    warning('ieeglab_detect_n1:lowSampleRate', ...
        ['Sampling rate is %g Hz (%.1f ms between samples). N1 latency is only resolved to that ' ...
         'precision, and residual artifact is hard to separate from an early response.'], ...
        EEG.srate, 1000/EEG.srate);
end
end

% ===================== local helpers =====================

function [stat, amp, lat, sd] = local_stat(A, iB, iN, tN, requirePeak, polarity, minSd)
% A: samples x patterns (column 1 is the observed average). Statistic per
% column: largest deflection in the window, from the baseline median, in units
% of baseline SD. polarity selects which deflections count ('abs', 'negative'
% or 'positive'); minSd is a floor on the baseline SD, so that threshold x SD
% cannot fall below threshold x minSd on a quiet contact (erdetect's rule).
if nargin < 6 || isempty(polarity), polarity = 'abs'; end
if nargin < 7 || isempty(minSd),    minSd    = 0;     end
A = A - median(A(iB,:), 1, 'omitnan');
sdAll = std(A(iB,:), 0, 1, 'omitnan');
sdAll = max(sdAll, minSd);
seg = A(iN,:);
switch polarity
    case 'negative', a = -seg; a(a < 0) = 0;
    case 'positive', a =  seg; a(a < 0) = 0;
    otherwise,       a = abs(seg);
end
n = size(seg, 1);
if requirePeak && n >= 3
    ext = false(size(a));
    ext(2:n-1,:) = a(2:n-1,:) >= a(1:n-2,:) & a(2:n-1,:) >= a(3:n,:);
    a(~ext) = 0;
end
[mx, k] = max(a, [], 1);
stat = mx ./ sdAll;
stat(~isfinite(stat)) = 0;
sd = sdAll(1);
if mx(1) > 0
    amp = seg(k(1), 1);
    lat = tN(k(1));
else
    amp = NaN; lat = NaN;                      % no peak inside the window
end
end

function [S, exact] = local_signs(K, nPerm)
% Sign patterns, K x P. Column 1 is the identity (the observed data). All 2^K
% patterns when that is no more than nPerm + 1, otherwise random ones.
exact = 2^K <= nPerm + 1;
if exact
    B = dec2bin(0:2^K-1, K) == '1';
    S = double(~B') * 2 - 1;          % first column: all +1
else
    S = [ones(K,1), (randi(2, K, nPerm) * 2 - 3)];
end
end

function padj = local_correct(p, method)
padj = p;
ok = ~isnan(p);
n = sum(ok);
if n == 0, return; end
switch lower(method)
    case 'bonferroni'
        padj(ok) = min(1, p(ok)*n);
    case 'fdr'
        [ps, ord] = sort(p(ok));
        q = ps(:) .* n ./ (1:n)';
        q = min(1, flipud(cummin(flipud(q))));
        tmp = nan(n,1); tmp(ord) = q;
        padj(ok) = tmp;
    case 'none'
    otherwise
        warning('ieeglab_detect_n1:badCorrection', 'Unknown correction "%s"; reporting uncorrected p-values.', method);
end
end
