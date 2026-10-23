function files = ieeglab_remembered_surfaces(EEG, files)
% ieeglab_remembered_surfaces() - Brain surface files remembered for a dataset.
%
% Usage:
%   files = ieeglab_remembered_surfaces(EEG)      % recall
%   ieeglab_remembered_surfaces(EEG, files)       % remember
%
% Recall returns the surfaces stored on the dataset (EEG.ieeglab.opt.surf_files),
% else the ones last used for the same subject of the same dataset (BIDS root
% folder) in any MATLAB session (MATLAB preference iEEGLAB/surf_files): sub-02 of
% another dataset is another person. Files that no longer exist are dropped.
%
% Cedric Cannard, iEEGLAB, 2026

if nargin > 1
    try
        setpref('iEEGLAB', 'surf_files', struct('subject', local_subject(EEG), 'root', local_root(EEG), ...
            'files', {cellstr(files)}));
    catch
    end
    return
end

files = {};
if isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'opt') && isstruct(EEG.ieeglab.opt) && ...
   isfield(EEG.ieeglab.opt,'surf_files') && ~isempty(EEG.ieeglab.opt.surf_files)
    files = cellstr(EEG.ieeglab.opt.surf_files);
else
    try
        last = getpref('iEEGLAB', 'surf_files');
        sameRoot = isfield(last, 'root') && strcmpi(char(last.root), local_root(EEG));   % both '' when no folder is known
        if isequal(last.subject, local_subject(EEG)) && sameRoot, files = cellstr(last.files); end
    catch
    end
end
files = files(cellfun(@(f) exist(f,'file') == 2, files));
files = files(:)';
end

function s = local_subject(EEG)
s = '';
if isfield(EEG,'subject') && ~isempty(EEG.subject), s = char(string(EEG.subject)); end
if isempty(s) && isfield(EEG,'filename'), s = char(EEG.filename); end
tok = regexp(s, '^sub-([A-Za-z0-9]+)', 'tokens', 'once');   % also from a BIDS file name
if ~isempty(tok), s = tok{1}; end
end

function r = local_root(EEG)
% The dataset's root folder: above sub-<label> (and ses-/ieeg folders), or the data folder
r = '';
if isfield(EEG,'filepath') && ~isempty(EEG.filepath), r = char(EEG.filepath); end
if isempty(r) && isfield(EEG,'BIDS') && isstruct(EEG.BIDS) && isfield(EEG.BIDS,'sourcefile') && ~isempty(EEG.BIDS.sourcefile)
    r = fileparts(char(EEG.BIDS.sourcefile));
end
p = r;
for k = 1:4
    [up, name] = fileparts(p);
    if isempty(name) || strcmp(up, p), break; end
    if startsWith(name, 'sub-'), r = up; break; end
    p = up;
end
end
