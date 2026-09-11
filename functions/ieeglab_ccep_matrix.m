function [EEG, M, com] = ieeglab_ccep_matrix(EEG, opt)
% ieeglab_ccep_matrix() - The CCEP connectivity matrix: stimulation sites x
%                         recording contacts.
%
% Usage:
%   EEG         = ieeglab_ccep_matrix(EEG)
%   [EEG, M]    = ieeglab_ccep_matrix(EEG, opt)
%
% The standard summary of a single-pulse stimulation experiment: one row per
% stimulated pair, one column per recording contact, and in each cell whether
% stimulating that pair evoked a significant response at that contact. Cells
% for the stimulated contacts themselves, marked-bad contacts and pairs not
% tested are NaN - "not measured" stays distinct from "no response", including
% in the degrees.
%
% Options:
%   .source           'n1' (default) or 'crp' - which detector defines a response
%   .compute_missing  run the detector when its results are absent, or were
%                     computed with options other than those passed here.
%                     Default true. Any N1 option (n1_window, baseline, method,
%                     n_perm, threshold, alpha, correct, min_trials,
%                     exclude_stim, require_peak) or CRP option (crp_window, ...)
%                     given here is forwarded.
%   .plot             draw the matrix. Default false.
%   .verbose          default true
%
% Output M (also EEG.ieeglab.ccep_matrix):
%   .sites        stimulation sites, in montage order of their contacts
%   .channels     recording contacts, montage order
%   .response     1 significant, 0 tested and not significant, NaN not tested
%   .p_adj        corrected p-values (NaN for the 'sd' rule)
%   .amplitude_uv, .latency_ms       when source = 'n1' (all tested pairs)
%   .tR_ms, .explained_var           when source = 'crp'
%   .out_degree   contacts each site drives (NaN if the site tested nothing)
%   .in_degree    sites that evoke a response at each contact (NaN if never tested)
%   .n_tested, .n_significant, .density, .source
%   .source_opt, .source_height      what the matrix was built from; export and
%                                    plots use them to detect a stale matrix
%
% Cedric Cannard, iEEGLAB, 2026

com = '';
if nargin < 2 || isempty(opt), opt = struct(); end
def = struct('source','n1', 'compute_missing',true, 'plot',false, 'verbose',true);
f = fieldnames(def);
for i = 1:numel(f)
    if ~isfield(opt,f{i}) || isempty(opt.(f{i})), opt.(f{i}) = def.(f{i}); end
end
src = lower(char(string(opt.source)));
n1Names  = {'n1_window','baseline','method','n_perm','threshold','alpha','correct','min_trials','exclude_stim','require_peak'};
crpNames = {'crp_window','alpha','correct','min_trials','exclude_stim','n_perm'};

switch src
    case 'n1'
        have = isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'n1') && ...
               isfield(EEG.ieeglab.n1,'table') && ~isempty(EEG.ieeglab.n1.table);
        requested = local_pass(opt, n1Names);
        if have && local_differs(requested, EEG.ieeglab.n1.opt), have = false; end
        if ~have
            if ~opt.compute_missing
                error('ieeglab_ccep_matrix:noN1', 'No current N1 results on the dataset. Run ieeglab_detect_n1 first.');
            end
            requested.verbose = opt.verbose;
            EEG = ieeglab_detect_n1(EEG, requested);
        end
        T = EEG.ieeglab.n1.table; srcOpt = EEG.ieeglab.n1.opt;
        metrics = {'amplitude_uv','n1_amplitude_uv'; 'latency_ms','n1_latency_ms'};
    case 'crp'
        have = isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'stats') && ...
               isfield(EEG.ieeglab.stats,'table') && ~isempty(EEG.ieeglab.stats.table);
        requested = local_pass(opt, crpNames);
        if have && local_differs(requested, EEG.ieeglab.stats.opt), have = false; end
        if ~have
            if ~opt.compute_missing
                error('ieeglab_ccep_matrix:noCRP', 'No current CRP results on the dataset. Run ieeglab_stats_subject first.');
            end
            requested.plot = false; requested.run_n1 = false; requested.run_matrix = false;
            requested.run_crp = true; requested.export_dir = ''; requested.verbose = opt.verbose;
            EEG = ieeglab_stats_subject(EEG, requested);
        end
        T = EEG.ieeglab.stats.table; srcOpt = EEG.ieeglab.stats.opt;
        metrics = {'tR_ms','tR_ms'; 'explained_var','explained_var'};
    otherwise
        error('ieeglab_ccep_matrix:badSource', 'source must be ''n1'' or ''crp''; got ''%s''.', src);
end
if isempty(T)
    error('ieeglab_ccep_matrix:empty', ...
        'The %s results table is empty - usually too few trials per stimulation site.', upper(src));
end

% ---------- axes ----------
channels = {EEG.chanlocs.labels}';
sites = unique(cellstr(string(T.site)), 'stable');
sites = local_order_sites(sites, channels);
nS = numel(sites); nC = numel(channels);

M = struct();
M.sites = sites;
M.channels = channels;
M.response = nan(nS, nC);
M.p_adj    = nan(nS, nC);
for k = 1:size(metrics,1), M.(metrics{k,1}) = nan(nS, nC); end

[~, si] = ismember(cellstr(string(T.site)), sites);
[~, ci] = ismember(upper(strtrim(cellstr(string(T.channel)))), upper(strtrim(channels)));
ok = si > 0 & ci > 0;
lin = sub2ind([nS nC], si(ok), ci(ok));
M.response(lin) = double(T.significant(ok));
M.p_adj(lin)    = T.p_adj(ok);
for k = 1:size(metrics,1)
    M.(metrics{k,1})(lin) = T.(metrics{k,2})(ok);
end

testedRow = any(~isnan(M.response), 2);
testedCol = any(~isnan(M.response), 1)';
M.out_degree = sum(M.response == 1, 2);   M.out_degree(~testedRow) = NaN;
M.in_degree  = sum(M.response == 1, 1)';  M.in_degree(~testedCol)  = NaN;
M.n_tested      = nnz(~isnan(M.response));
M.n_significant = nnz(M.response == 1);
M.density       = M.n_significant / max(M.n_tested, 1);
M.source        = src;
M.source_opt    = srcOpt;
M.source_height = height(T);

EEG.ieeglab.ccep_matrix = M;

if opt.verbose
    [~, top] = max(M.out_degree);
    fprintf(['[CCEP matrix] %d sites x %d contacts from %s: %d of %d tested pairs respond ' ...
             '(density %.2f). Most connected site: %s (%d contacts).\n'], ...
        nS, nC, upper(src), M.n_significant, M.n_tested, M.density, sites{top}, M.out_degree(top));
end
if opt.plot
    ieeglab_plot_ccep_matrix(M);
end
com = sprintf('EEG = ieeglab_ccep_matrix(EEG, %s);', ieeglab_literal(opt, {'plot'}));
end

% ======================= helpers =======================

function s = local_order_sites(s, channels)
% Order sites by the montage position of their first, then second, contact
% (tokens resolved as everywhere else: trimmed, case- and punctuation-tolerant),
% so the matrix shows lead-by-lead block structure.
key = inf(numel(s), 2);
[~, idx] = ieeglab_canonical_site(s, channels);
for i = 1:numel(s)
    p = sort(idx{i});
    if ~isempty(p), key(i,1) = p(1); end
    if numel(p) > 1, key(i,2) = p(2); end
end
[~, ord] = sortrows(key);
s = s(ord);
end

function o = local_pass(opt, names)
o = struct();
for i = 1:numel(names)
    if isfield(opt, names{i}) && ~isempty(opt.(names{i})), o.(names{i}) = opt.(names{i}); end
end
end

function tf = local_differs(requested, stored)
% True when any option explicitly requested differs from the one the stored
% results were computed with.
tf = false;
f = fieldnames(requested);
for i = 1:numel(f)
    if ~isfield(stored, f{i}), tf = true; return; end
    a = requested.(f{i}); b = stored.(f{i});
    if (ischar(a) || isstring(a)) || (ischar(b) || isstring(b))
        if ~strcmpi(char(string(a)), char(string(b))), tf = true; return; end
    elseif ~isequal(double(a), double(b))
        tf = true; return
    end
end
end
