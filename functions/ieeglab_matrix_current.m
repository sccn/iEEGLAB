function [tf, why] = ieeglab_matrix_current(EEG)
% ieeglab_matrix_current() - True when EEG.ieeglab.ccep_matrix still matches the
%                            N1 or CRP results it was built from.
%
% Usage:
%   [tf, why] = ieeglab_matrix_current(EEG)
%
% The matrix records the options and size of its source table. Re-running N1
% or CRP with other settings leaves a matrix that no longer describes the
% results next to it; export and the electrode plots check this first, so
% they never mix two analyses. why says what is wrong when tf is false.
%
% Cedric Cannard, iEEGLAB, 2026

tf = false;
why = 'there is no connectivity matrix on the dataset';
if ~isfield(EEG,'ieeglab') || ~isfield(EEG.ieeglab,'ccep_matrix') || isempty(EEG.ieeglab.ccep_matrix)
    return
end
M = EEG.ieeglab.ccep_matrix;
switch M.source
    case 'n1',  fld = 'n1';
    case 'crp', fld = 'stats';
    otherwise
        why = sprintf('the matrix names an unknown source "%s"', M.source);
        return
end
if ~isfield(EEG.ieeglab, fld) || ~isfield(EEG.ieeglab.(fld), 'table') || ~isfield(EEG.ieeglab.(fld), 'opt')
    why = sprintf('the %s results it was built from are no longer on the dataset', upper(M.source));
    return
end
S = EEG.ieeglab.(fld);
if ~isfield(M,'source_height') || ~isfield(M,'source_opt') || height(S.table) ~= M.source_height ...
        || ~isequaln(S.opt, M.source_opt)
    why = sprintf(['the %s results were recomputed after the matrix was built; ' ...
        'rebuild it with ieeglab_ccep_matrix'], upper(M.source));
    return
end
tf = true;
why = '';
end
