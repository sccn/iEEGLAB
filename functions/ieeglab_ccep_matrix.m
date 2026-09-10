function [EEG, M, com] = ieeglab_ccep_matrix(EEG, opt)
% ieeglab_ccep_matrix() - The CCEP connectivity matrix: stimulation sites x
%                         recording channels.
%
% Usage:
%   EEG         = ieeglab_ccep_matrix(EEG)
%   [EEG, M]    = ieeglab_ccep_matrix(EEG, opt)
%
% This is the standard summary of a single-pulse stimulation experiment: one
% row per stimulated pair, one column per recording contact, and in each cell
% whether stimulating that pair evoked a significant response at that contact.
% Cells for the stimulated contacts themselves, and for pairs not tested, are
% NaN - "not measured" is kept distinct from "no response".
%
% Options:
%   .source           'n1' (default) or 'crp' - which detector defines a response
%   .compute_missing  run the detector first if its results are absent. Default true.
%   .plot             draw the matrix. Default false.
%   .verbose          default true
%
% Output M (also stored in EEG.ieeglab.ccep_matrix):
%   .sites        nSites x 1 cellstr, ordered by the montage position of their contacts
%   .channels     nChans x 1 cellstr, montage order
%   .response     nSites x nChans: 1 significant, 0 tested and not significant, NaN not tested
%   .p_adj        corrected p-values
%   .amplitude_uv, .latency_ms       when source = 'n1'
%   .tR_ms, .explained_var           when source = 'crp'
%   .out_degree   nSites x 1: how many contacts each site drives
%   .in_degree    nChans x 1: how many sites evoke a response at each contact
%   .n_tested, .n_significant, .source
%
% Cedric Cannard, iEEGLAB, 2026

com = '';
if nargin < 2 || isempty(opt), opt = struct(); end
def = struct('source','n1', 'compute_missing',true, 'plot',false, 'verbose',true);
f = fieldnames(def);
for i = 1:numel(f)
    if ~isfield(opt,f{i}) || isempty(opt.(f{i})), opt.(f{i}) = def.(f{i}); end
end
src = lower(char(opt.source));

% ---------- results table ----------
switch src
    case 'n1'
        have = isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'n1') && ...
               isfield(EEG.ieeglab.n1,'table') && ~isempty(EEG.ieeglab.n1.table);
        if ~have
            if ~opt.compute_missing
                error('ieeglab_ccep_matrix:noN1', 'No N1 results on the dataset. Run ieeglab_detect_n1 first.');
            end
            [EEG, ~] = ieeglab_detect_n1(EEG, local_pass(opt, {'n1_window','baseline','threshold','method', ...
                'n_perm','alpha','correct','min_trials','exclude_stim','verbose'}));
        end
        T = EEG.ieeglab.n1.table;
        metrics = {'amplitude_uv','n1_amplitude_uv'; 'latency_ms','n1_latency_ms'};
    case 'crp'
        have = isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'stats') && ...
               isfield(EEG.ieeglab.stats,'table') && ~isempty(EEG.ieeglab.stats.table);
        if ~have
            if ~opt.compute_missing
                error('ieeglab_ccep_matrix:noCRP', 'No CRP results on the dataset. Run ieeglab_stats_subject first.');
            end
            o = local_pass(opt, {'crp_window','alpha','min_trials','correct','exclude_stim','verbose'});
            o.plot = false; o.run_n1 = false; o.run_matrix = false; o.export_dir = '';
            EEG = ieeglab_stats_subject(EEG, o);
        end
        T = EEG.ieeglab.stats.table;
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
[~, ci] = ismember(upper(cellstr(string(T.channel))), upper(channels));
ok = si > 0 & ci > 0;
lin = sub2ind([nS nC], si(ok), ci(ok));
M.response(lin) = double(T.significant(ok));
M.p_adj(lin)    = T.p_adj(ok);
for k = 1:size(metrics,1)
    M.(metrics{k,1})(lin) = T.(metrics{k,2})(ok);
end

M.out_degree    = sum(M.response == 1, 2);
M.in_degree     = sum(M.response == 1, 1)';
M.n_tested      = nnz(~isnan(M.response));
M.n_significant = nnz(M.response == 1);
M.density       = M.n_significant / max(M.n_tested, 1);
M.source        = src;

EEG.ieeglab.ccep_matrix = M;

if opt.verbose
    [~, top] = max(M.out_degree);
    fprintf(['[CCEP matrix] %d sites x %d channels from %s: %d of %d tested pairs respond ' ...
             '(density %.2f). Most connected site: %s (%d contacts).\n'], ...
        nS, nC, upper(src), M.n_significant, M.n_tested, M.density, sites{top}, M.out_degree(top));
end
if opt.plot
    ieeglab_plot_ccep_matrix(M);
end
com = sprintf('EEG = ieeglab_ccep_matrix(EEG, struct(''source'',''%s''));', src);
end

% ======================= helpers =======================

function s = local_order_sites(s, channels)
% Order sites by the montage position of their first, then second, contact,
% so the matrix shows lead-by-lead block structure.
U = upper(channels);
key = inf(numel(s), 2);
for i = 1:numel(s)
    [tf, loc] = ismember(upper(ieeglab_site_tokens(s{i})), U);
    p = sort(loc(tf));
    if ~isempty(p), key(i,1) = p(1); end
    if numel(p) > 1, key(i,2) = p(2); end
end
[~, ord] = sortrows(key);
s = s(ord);
end

function o = local_pass(opt, names)
o = struct();
for i = 1:numel(names)
    if isfield(opt, names{i}), o.(names{i}) = opt.(names{i}); end
end
end
