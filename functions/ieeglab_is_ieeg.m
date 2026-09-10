function tf = ieeglab_is_ieeg(EEG)
% ieeglab_is_ieeg() - True when a dataset is intracranial EEG.
%
% Usage:
%   tf = ieeglab_is_ieeg(EEG)
%
% True when the dataset has been through iEEGLAB (EEG.ieeglab exists) or when
% at least half of its channels are typed SEEG, ECOG or DBS in EEG.chanlocs.
%
% Used to route EEGLAB-style topography requests to ieeglab_topoplot: a 2D
% scalp interpolation is not meaningful for electrodes that sit inside the
% brain, so iEEG data should be drawn as electrodes on the brain instead.
%
% Cedric Cannard, iEEGLAB, 2026

tf = false;
if ~isstruct(EEG), return; end
if isfield(EEG,'ieeglab') && ~isempty(EEG.ieeglab)
    tf = true;
    return
end
if isfield(EEG,'chanlocs') && ~isempty(EEG.chanlocs) && isfield(EEG.chanlocs,'type')
    ty = cellfun(@local_upper, {EEG.chanlocs.type}, 'UniformOutput', false);
    tf = mean(ismember(ty, {'SEEG','ECOG','DBS','IEEG'})) >= 0.5;
end
end

function s = local_upper(x)
if isempty(x), s = ''; else, s = upper(strtrim(char(string(x)))); end
end
