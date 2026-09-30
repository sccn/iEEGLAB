function ieeglab_require_data(EEG, caller)
% ieeglab_require_data() - Error clearly when no dataset is loaded.
%
% Usage:
%   ieeglab_require_data(EEG, 'ieeglab_load')
%
% EEGLAB passes EEG = [] to menu callbacks when nothing is loaded, which
% otherwise fails later with "Dot indexing is not supported for variables of
% type double".
%
% Cedric Cannard, iEEGLAB, 2026

if nargin < 2, caller = 'iEEGLAB'; end
if isempty(EEG) || ~isstruct(EEG) || ~isfield(EEG, 'data') || isempty(EEG.data)
    error([caller ':noData'], ['%s: no dataset is loaded. Import or load the iEEG data first ' ...
        '(File > Import data > From BIDS folder structure, or File > Load existing dataset).'], caller);
end
end
