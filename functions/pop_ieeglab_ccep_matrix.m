function com = pop_ieeglab_ccep_matrix(EEG, metric)
% pop_ieeglab_ccep_matrix() - Choose what the CCEP connectivity matrix shows,
%                             then draw it (ieeglab_plot_ccep_matrix).
%
% Usage:
%   pop_ieeglab_ccep_matrix(EEG)                    % dialog
%   pop_ieeglab_ccep_matrix(EEG, 'explained_var')
%
% Offers what the dataset holds: significant responses, N1 amplitude and
% latency (N1 matrix), CRP response duration and explained variance (CRP
% results, also when the matrix itself was built from N1).
%
% Cedric Cannard, iEEGLAB, 2026

com = '';
[ok, why] = ieeglab_matrix_current(EEG);
if ~ok
    error('pop_ieeglab_ccep_matrix:noMatrix', ...
        'Cannot plot the connectivity matrix: %s. Run iEEGLAB > CCEP analysis.', why);
end
M = EEG.ieeglab.ccep_matrix;
names = {'response','amplitude_uv','latency_ms','tR_ms','explained_var'};
descs = {'Significant responses (which pairs respond)', 'N1 amplitude (signed, uV)', ...
         'N1 latency (ms)', 'CRP response duration tau_R (ms)', 'CRP explained variance (0-1)'};
fromCRP = false(size(names));
for k = 4:5
    if ~isfield(M, names{k}) && isfield(EEG.ieeglab,'stats') && isfield(EEG.ieeglab.stats,'table') ...
            && ~isempty(EEG.ieeglab.stats.table)
        fromCRP(k) = true;                            % N1 matrix, CRP results available
    end
end
avail = cellfun(@(n) isfield(M, n), names) | fromCRP;

if nargin < 2 || isempty(metric)
    idx = find(avail);
    [sel, okSel] = listdlg('ListString', descs(idx), 'SelectionMode', 'single', ...
        'Name', 'iEEGLAB', 'PromptString', 'Connectivity matrix: what to show', 'ListSize', [320 120]);
    if ~okSel || isempty(sel), return; end
    metric = names{idx(sel)};
end
k = find(strcmp(names, metric), 1);
if ~isempty(k) && fromCRP(k)
    [~, M] = ieeglab_ccep_matrix(EEG, struct('source','crp', 'compute_missing',false, 'verbose',false));
end
ieeglab_plot_ccep_matrix(M, metric);
com = sprintf('pop_ieeglab_ccep_matrix(EEG, ''%s'');', metric);
end
