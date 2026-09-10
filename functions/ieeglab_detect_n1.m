function [EEG, T] = ieeglab_detect_n1(EEG, opt)
% ieeglab_detect_n1() - Detect N1 evoked responses in CCEP data.
%
% For every stimulation site x recording channel, tests whether the trial-
% averaged response contains a significant early deflection, and reports its
% amplitude and latency. This is the most widely reported CCEP measure and
% complements the CRP fit in ieeglab_stats_subject: N1 asks "is there a
% response and when does it peak", CRP asks "what shape is it and how long
% does it last".
%
% Usage:
%   EEG      = ieeglab_detect_n1(EEG)
%   [EEG, T] = ieeglab_detect_n1(EEG, opt)
%
% Options:
%   .n1_window     [1x2] ms, search window for the peak. Default [10 100].
%                  Must start after the blanked stimulation artifact.
%   .baseline      [1x2] ms, window for the noise estimate. Default [-500 -10].
%   .threshold     multiples of the baseline SD a peak must exceed. Default 3.4,
%                  following van Blooijs et al. Ignored when method='permutation'.
%   .method        'sd'          - peak must exceed threshold x baseline SD (default)
%                  'permutation' - non-parametric null built by sign-flipping
%                                  trials; slower but makes no distributional
%                                  assumption
%   .n_perm        permutations when method='permutation'. Default 500.
%   .alpha         significance level. Default 0.05.
%   .correct       'fdr' (default), 'bonferroni' or 'none', applied across
%                  channels within each stimulation site.
%   .min_trials    minimum trials per site. Default 5.
%   .exclude_stim  skip the stimulated contacts for their own site. Default true.
%   .verbose       logical. Default true.
%
% Output:
%   T - table with site, channel, n1_amplitude_uv, n1_latency_ms, baseline_sd,
%       z, p, p_adj, significant
%   Also stored in EEG.ieeglab.n1
%
% Sign convention: the N1 is a negative deflection by convention, but polarity
% depends on the reference, so the largest ABSOLUTE deflection in the window is
% taken and its sign reported in n1_amplitude_uv.
%
% References:
%   van Blooijs D, et al. (2018). Evoked directional network characteristics of
%   epileptogenic tissue derived from single pulse electrical stimulation.
%   Human Brain Mapping 39(11):4611-4622.
%   Ojeda Valencia G, et al. (2023). Signatures of electrical stimulation driven
%   network interactions in the human limbic system. J Neurosci 43(39):6697-6711.
%
% Cedric Cannard, iEEGLAB, 2026

if nargin < 2, opt = struct(); end
def = struct('n1_window',[10 100], 'baseline',[-500 -10], 'threshold',3.4, ...
             'method','sd', 'n_perm',500, 'alpha',0.05, 'correct','fdr', ...
             'min_trials',5, 'exclude_stim',true, 'verbose',true);
fn = fieldnames(def);
for i = 1:numel(fn)
    if ~isfield(opt,fn{i}) || isempty(opt.(fn{i})), opt.(fn{i}) = def.(fn{i}); end
end

if EEG.trials < 2
    error('ieeglab_detect_n1:notEpoched', ...
        'N1 detection needs epoched data with at least 2 trials (this dataset has %d).', EEG.trials);
end

dmode = ieeglab_detect_mode(EEG);
if ~strcmp(dmode,'ccep')
    warning('ieeglab_detect_n1:notCCEP', ...
        ['N1 detection is defined for CCEP data; this dataset looks like "%s". ' ...
         'Results will be grouped by event type but the "stimulation site" ' ...
         'interpretation does not apply.'], dmode);
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
    error('ieeglab_detect_n1:badBaseline', ...
        'Baseline [%g %g] ms selects only %d samples; need at least 10 for a noise estimate.', ...
        opt.baseline(1), opt.baseline(2), nnz(iB));
end

labels = string({EEG.chanlocs.labels});
sites  = local_sites(EEG, labels);
[uSites, ~, grp] = unique(sites);

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
    if opt.exclude_stim && uSites(g) ~= ""
        stimIdx = find(ismember(upper(labels), upper(split(uSites(g),'-'))'));
    end

    pvals = nan(EEG.nbchan,1);
    tmp   = cell(EEG.nbchan,1);
    for ch = 1:EEG.nbchan
        if ismember(ch, stimIdx), continue; end
        trials = squeeze(EEG.data(ch,:,tr));           % [T x nTrials]
        if size(trials,2) < 2, continue; end
        avg = mean(trials, 2, 'omitnan');

        bslSD = std(avg(iB), 0, 'omitnan');
        if ~isfinite(bslSD) || bslSD == 0, continue; end

        seg = avg(iN);
        [~, k] = max(abs(seg));
        amp = seg(k);
        tN  = t(iN);
        lat = tN(k);
        z   = abs(amp) / bslSD;

        switch lower(opt.method)
            case 'permutation'
                % Null: flip the sign of a random subset of trials, so any
                % consistent evoked deflection cancels while noise does not.
                null = zeros(opt.n_perm,1);
                for pI = 1:opt.n_perm
                    sgn = (randi(2, 1, size(trials,2))*2 - 3);     % +1 / -1
                    a   = mean(trials .* sgn, 2, 'omitnan');
                    s   = std(a(iB), 0, 'omitnan');
                    if ~isfinite(s) || s == 0, null(pI) = 0; continue; end
                    null(pI) = max(abs(a(iN))) / s;
                end
                p = (1 + sum(null >= z)) / (opt.n_perm + 1);
            otherwise
                % Threshold on the baseline SD, converted to a normal-tail p so
                % that the same correction machinery applies.
                p = 2 * (1 - normcdf(z));
                if z < opt.threshold, p = max(p, opt.alpha*1.001); end
        end

        pvals(ch) = p;
        tmp{ch} = struct('amp',amp, 'lat',lat, 'sd',bslSD, 'z',z);
    end

    padj = local_correct(pvals, opt.correct);
    for ch = 1:EEG.nbchan
        if isempty(tmp{ch}), continue; end
        r = tmp{ch};
        rows(end+1,:) = { char(uSites(g)), char(labels(ch)), r.amp, r.lat, ...
            r.sd, r.z, pvals(ch), padj(ch), padj(ch) < opt.alpha }; %#ok<AGROW>
    end
end

if isempty(rows)
    warning('ieeglab_detect_n1:noResults', ...
        'No N1 fits succeeded. Most often too few trials per site (minimum %d).', opt.min_trials);
    T = table();
    EEG.ieeglab.n1 = struct('table', T, 'opt', opt);
    return
end

T = cell2table(rows, 'VariableNames', ...
    {'site','channel','n1_amplitude_uv','n1_latency_ms','baseline_sd','z','p','p_adj','significant'});
T = sortrows(T, {'site','p_adj'});
EEG.ieeglab.n1 = struct('table', T, 'opt', opt);

if opt.verbose
    fprintf('\n[N1] %d significant of %d site-channel pairs (%s-corrected, alpha=%g, method=%s).\n', ...
        sum(T.significant), height(T), opt.correct, opt.alpha, opt.method);
    if any(T.significant)
        s = T(T.significant,:);
        fprintf('[N1] Latency %.0f-%.0f ms (median %.0f), |amplitude| %.0f-%.0f uV.\n', ...
            min(s.n1_latency_ms), max(s.n1_latency_ms), median(s.n1_latency_ms), ...
            min(abs(s.n1_amplitude_uv)), max(abs(s.n1_amplitude_uv)));
    end
end

% Plausibility check. A genuine CCEP N1 is on the order of tens to a few
% hundred microvolts. Amplitudes in the thousands mean the detector has locked
% onto residual stimulation artifact rather than a neural response - usually
% because the search window opens too early, the artifact was not blanked, or
% a long zero-phase filter rang forward into the window.
BIG = 1000;   % uV
nBig = sum(abs(T.n1_amplitude_uv) > BIG);
if nBig > 0
    warning('ieeglab_detect_n1:implausibleAmplitude', ...
        ['%d of %d detected peaks exceed %g uV, which is far above a physiological ' ...
         'CCEP N1 (tens to a few hundred uV). These are very likely residual ' ...
         'stimulation artifact, not neural responses.\n' ...
         'Check, in order: (1) was the artifact blanked BEFORE filtering and before ' ...
         'downsampling; (2) does n1_window start late enough (currently %g ms); ' ...
         '(3) is a long zero-phase high-pass ringing into the window.'], ...
         nBig, height(T), BIG, opt.n1_window(1));
end
if EEG.srate < 500
    warning('ieeglab_detect_n1:lowSampleRate', ...
        ['Sampling rate is %g Hz, giving %.1f ms between samples. N1 latency can only ' ...
         'be resolved to that precision, and residual stimulation artifact is hard to ' ...
         'separate from an early response at this rate.'], EEG.srate, 1000/EEG.srate);
end
end

% ===================== local helpers =====================

function sites = local_sites(EEG, labels)
N = EEG.trials;
sites = strings(1,N);
if ~isfield(EEG,'epoch') || isempty(EEG.epoch) || numel(EEG.epoch) ~= N, return; end
for i = 1:N
    ty = EEG.epoch(i).eventtype;
    l  = [];
    if isfield(EEG.epoch,'eventlatency'), l = EEG.epoch(i).eventlatency; end
    if iscell(ty)
        if iscell(l) && numel(l) == numel(ty)
            [~,k] = min(cellfun(@(x) abs(double(x(1))), l));
        else
            k = 1;
        end
        ty = ty{k};
    end
    if isnumeric(ty), ty = num2str(ty); end
    ty = strtrim(string(ty));
    if ty == "", continue; end
    parts = regexp(char(ty), '[-+/|]', 'split');
    parts = parts(~cellfun(@isempty, parts));
    hit = parts(ismember(upper(parts), upper(labels)));
    if numel(hit) >= 2
        sites(i) = string(strjoin(sort(hit), '-'));   % order-independent
    else
        sites(i) = ty;
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
        padj(ok) = min(1, p(ok)*n);
    case 'fdr'
        [ps, ord] = sort(p(ok));
        q = ps(:) .* n ./ (1:n)';
        q = min(1, flipud(cummin(flipud(q))));
        tmp = nan(n,1); tmp(ord) = q;
        padj(ok) = tmp;
    case 'none'
    otherwise
        warning('ieeglab_detect_n1:badCorrection', ...
            'Unknown correction "%s"; reporting uncorrected p-values.', method);
end
end
