function [files, where] = ieeglab_fsaverage_surfaces()
% ieeglab_fsaverage_surfaces() - fsaverage's pial surfaces, if installed.
%
% Usage:
%   [files, where] = ieeglab_fsaverage_surfaces()   % {} and '' when not found
%
% fsaverage is FreeSurfer's average brain, the surface on which coordinates in
% fsaverage (MNI305) space are defined. It is looked for in
% $SUBJECTS_DIR/fsaverage, $FREESURFER_HOME/subjects/fsaverage and MNE's copy
% of it (~/mne_data/MNE-fsaverage-data/fsaverage, from
% mne.datasets.fetch_fsaverage). files are surf/lh.pial and surf/rh.pial
% (FreeSurfer format, read by ieeglab_read_gifti). Which brain suits a
% recording is decided by ieeglab_brain_surfaces.
%
% Cedric Cannard, iEEGLAB, 2026

files = {}; where = '';
home = getenv('USERPROFILE'); if isempty(home), home = getenv('HOME'); end
cands = {};
if ~isempty(getenv('SUBJECTS_DIR')),    cands{end+1} = fullfile(getenv('SUBJECTS_DIR'), 'fsaverage'); end
if ~isempty(getenv('FREESURFER_HOME')), cands{end+1} = fullfile(getenv('FREESURFER_HOME'), 'subjects', 'fsaverage'); end
if ~isempty(home), cands{end+1} = fullfile(home, 'mne_data', 'MNE-fsaverage-data', 'fsaverage'); end
for k = 1:numel(cands)
    f = fullfile(cands{k}, 'surf', {'lh.pial', 'rh.pial'});
    if all(cellfun(@(x) exist(x, 'file') == 2, f))
        files = f; where = cands{k}; return
    end
end
end
