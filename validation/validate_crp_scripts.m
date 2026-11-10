function T = validate_crp_scripts(crpDir)
% validate_crp_scripts() - iEEGLAB's canonical response parameterization
%       (run_CRP) against the authors' code (CRP_method.m; Miller et al., 2023)
%       on the authors' own sample data, every file of it.
%
% Usage:
%   T = validate_crp_scripts('C:\data\crp_scripts')
%
% crpDir  the authors' code and data (osf.io/tx3yq, crp_scripts_Mar31_2023.zip;
%         also github.com/kaijmiller/crp_scripts): CRP_method.m and
%         sampledata/*.mat, each holding data (samples x trials, uV) and t (s).
%
% As the authors' CRP_illustration.m, each file is cut to 15 ms < t <= 1 s
% (or its end, when shorter) and the same matrix goes to both functions. The
% illustration's optional rejection of outlier trials is not applied: what is
% compared is the method. Reported per file: response duration tau_R, the
% mean projection profile and its variance, the projections, the canonical
% shape C (whose sign is arbitrary: aligned first), the trial weights alpha,
% the residuals, the explained variance, and the t-tests at tau_R and over the
% whole window (t and p; p at tau_R is the published significance), as the
% largest difference relative to the authors' values.
% run_CRP computes the projections in its own way (running Gram matrices
% instead of the authors' loop over durations) and its p-values with betainc;
% the parameterization at tau_R (eigenvectors of V'V, C, alpha, residuals)
% follows the authors' lines, so agreement there shows the same input and
% duration, not an independent implementation. iEEGLAB's pipeline tests CRP by
% permutation by default (ieeglab_stats_subject), which is not compared here,
% and cuts its window as 15 <= t <= 1000 ms (one more sample when 15 ms falls
% on a sample).
%
% Reference:
%   Miller, K. J., Muller, K.-R., Ojeda Valencia, G., Huang, H., Gregg, N. M.,
%     Worrell, G. A., & Hermes, D. (2023). Canonical response parameterization:
%     Quantifying the structure of responses to single-pulse intracranial
%     electrical brain stimulation. PLOS Computational Biology, 19(5), e1011105.
%
% Cedric Cannard, iEEGLAB, 2026

if exist(fullfile(crpDir, 'CRP_method.m'), 'file') ~= 2
    error('validate_crp_scripts:noCode', 'No CRP_method.m in %s.', crpDir);
end
addpath(crpDir);
cleanup = onCleanup(@() rmpath(crpDir));
files = dir(fullfile(crpDir, 'sampledata', '*.mat'));
files = files(~startsWith({files.name}, '._'));
rows = cell(numel(files), 1);
rel = @(a, b) max(abs(a(:) - b(:))) / max(max(abs(b(:))), eps);
for i = 1:numel(files)
    s = load(fullfile(files(i).folder, files(i).name), 'data', 't');
    t = s.t(:)';
    keep = t > 0.015 & t <= 1;
    V = s.data(keep, :); tw = t(keep);
    [pA, jA] = CRP_method(V, tw);                          % authors
    [pB, jB] = run_CRP(V, tw, struct('verbose', false));   % iEEGLAB
    sameTR = jA.tR_index == jB.tR_index;
    r = struct();
    r.file = string(erase(files(i).name, '.mat'));
    r.n_trials = size(V, 2);
    r.tR_ms_authors = 1000 * pA.tR;
    r.tR_ms_ieeglab = 1000 * pB.tR;
    r.same_tR = sameTR;
    r.profile_rel = rel(jB.mean_proj_profile, jA.mean_proj_profile);
    r.profile_var_rel = rel(jB.var_proj_profile, jA.var_proj_profile);
    r.projections_rel = rel(jB.S_all, jA.S_all);
    r.C_rel = NaN; r.alpha_rel = NaN; r.residual_rel = NaN; r.expl_var_max_diff = NaN;
    if sameTR                                              % (different lengths otherwise)
        sg = sign(pA.C' * pB.C);
        r.C_rel = rel(sg * pB.C, pA.C);
        r.alpha_rel = rel(sg * pB.al, pA.al);
        r.residual_rel = rel(pB.ep, pA.ep);
        r.expl_var_max_diff = max(abs(pB.expl_var(:) - pA.expl_var(:)));
    end
    r.t_tR_authors = jA.t_value_tR;
    r.t_rel = abs(jB.t_value_tR - jA.t_value_tR) / abs(jA.t_value_tR);
    r.p_tR_authors = jA.p_value_tR;
    r.p_tR_ieeglab = jB.p_value_tR;
    r.t_full_rel = abs(jB.t_value_full - jA.t_value_full) / abs(jA.t_value_full);
    r.p_full_diff = abs(jB.p_value_full - jA.p_value_full);
    rows{i} = struct2table(r);
end
T = vertcat(rows{:});
fprintf(['\n%d sample files: tau_R the same in %d. Largest relative difference over all files: profile %.1e ' ...
    '(variance %.1e), projections %.1e, t at tau_R %.1e, t full %.1e; p at tau_R %.1e, p full %.1e. Same ' ...
    'computation once tau_R agrees: C %.1e, alpha %.1e, residuals %.1e, explained variance %.1e.\n'], ...
    height(T), nnz(T.same_tR), max(T.profile_rel), max(T.profile_var_rel), max(T.projections_rel), max(T.t_rel), ...
    max(T.t_full_rel), max(abs(T.p_tR_ieeglab - T.p_tR_authors)), max(T.p_full_diff), max(T.C_rel), ...
    max(T.alpha_rel), max(T.residual_rel), max(T.expl_var_max_diff));
disp(T(:, {'file', 'n_trials', 'tR_ms_authors', 'tR_ms_ieeglab', 'p_tR_authors', 'p_tR_ieeglab'}));
end
