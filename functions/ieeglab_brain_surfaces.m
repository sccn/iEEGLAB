function [files, info] = ieeglab_brain_surfaces(EEG, varargin)
% ieeglab_brain_surfaces() - The brain surfaces that match the space of the
%                            electrode coordinates, and why.
%
% Usage:
%   [files, info] = ieeglab_brain_surfaces(EEG)
%   [files, info] = ieeglab_brain_surfaces(EEG, 'check_fit', false)
%
% The space comes from the BIDS *_coordsystem.json next to the recording
% (iEEGCoordinateSystem). Then:
%   participant space (ACPC, ScanRAS, Other, individual, T1w, ... or not given)
%       the participant's own pial surfaces: GIfTI *pial*.gii next to the data,
%       or in derivatives/freesurfer/sub-<label> (GIfTI, or FreeSurfer
%       surf/lh.pial and rh.pial). A template
%       would not line up: none is chosen when they are missing.
%   fsaverage (fsaverage, fsaverageSym, fsaverage3-6, MNI305)
%       fsaverage's pial surfaces (FreeSurfer or MNE installation; see
%       ieeglab_fsaverage_surfaces). Coordinates in this space are fsaverage
%       surface coordinates: they line up exactly.
%   MNI152 and other MNI templates
%       fsaverage's pial surfaces too, as an approximation: MNI305 and MNI152
%       differ by a few millimetres.
% When the space is not given and the participant has no surfaces, fsaverage is
% chosen only if the contacts fit it.
%
% The chosen surfaces are checked against the contacts (ieeglab_check_coords):
% when most contacts lie outside them, info.fit_ok is false and info.message
% says so.
%
% info fields: space (as in coordsystem.json, '' if not given), kind
% ('participant', 'fsaverage', 'fsaverage_approx' or 'none'), how (short
% description for a dialog), message (one or two sentences for the user: what
% was chosen and why, or what to do), folder, template_ok (the coordinates are
% in a template space, so a template brain is a fair approximation), fit_ok, fit
% (see ieeglab_check_coords).
%
% Cedric Cannard, iEEGLAB, 2026

p = inputParser;
p.addParameter('check_fit', true);
p.parse(varargin{:});
o = p.Results;

info = struct('space', local_coordsystem(EEG), 'kind', 'none', 'how', '', 'message', '', ...
    'folder', '', 'template_ok', false, 'fit_ok', true, 'fit', []);
sp = info.space;
isFs  = any(strcmpi(sp, {'fsaverage','fsaverageSym','fsaverage3','fsaverage4','fsaverage5','fsaverage6','MNI305'}));
isMni = ~isFs && (startsWith(sp, 'MNI', 'IgnoreCase', true) || startsWith(sp, 'ICBM', 'IgnoreCase', true) || ...
                  strcmpi(sp, 'IXI549Space') || strcmpi(sp, 'fsLR'));
% other BIDS template spaces: no surfaces for them here, the MNI template is an approximation
isOther = ~isFs && ~isMni && (any(strcmpi(sp, {'Talairach','IXI549Space'})) || startsWith(sp, 'NIHPD', 'IgnoreCase', true) || ...
          startsWith(sp, 'OASIS30', 'IgnoreCase', true) || startsWith(sp, 'UNCInfant', 'IgnoreCase', true));
info.template_ok = isFs || isMni || isOther;   % a template brain is a fair approximation
spaceTxt = sp; if isempty(spaceTxt), spaceTxt = 'not given in coordsystem.json'; end

[own, ownDir] = local_participant_surfaces(EEG);
[fs, fsDir] = ieeglab_fsaverage_surfaces();
getFs = ['fsaverage comes with FreeSurfer ($FREESURFER_HOME/subjects/fsaverage), or from Python: ' ...
         'import mne; mne.datasets.fetch_fsaverage().'];
files = {};
if isFs
    if ~isempty(fs)
        files = fs; info.kind = 'fsaverage'; info.folder = fsDir;
        info.how = 'the fsaverage pial surfaces (the coordinates are in fsaverage space)';
        info.message = sprintf('Coordinates in %s space: drawn on fsaverage''s pial surfaces (%s).', sp, fsDir);
    else
        info.message = sprintf(['Coordinates in %s space, but fsaverage was not found. %s ' ...
            'Until then, the template brain is an approximation.'], sp, getFs);
    end
elseif isMni
    if ~isempty(fs)
        files = fs; info.kind = 'fsaverage_approx'; info.folder = fsDir;
        info.how = sprintf('the fsaverage pial surfaces (the coordinates are in %s: close, not exact)', sp);
        info.message = sprintf(['Coordinates in %s space: drawn on fsaverage''s pial surfaces (%s), which are in ' ...
            'MNI305; the two templates differ by a few millimetres.'], sp, fsDir);
    else
        info.message = sprintf('Coordinates in %s space: drawn on the template brain (MNI). %s', sp, getFs);
    end
elseif isOther
    info.message = sprintf(['Coordinates in the %s template space: no surfaces for it here; drawn on the ' ...
        'template brain (MNI), an approximation.'], sp);
elseif ~isempty(own)
    files = own; info.kind = 'participant'; info.folder = ownDir;
    info.how = 'the participant''s pial surfaces';
    info.message = sprintf('Coordinates in the participant''s space (%s): drawn on the participant''s pial surfaces (%s).', ...
        spaceTxt, ownDir);
elseif isempty(sp) && ~isempty(fs) && o.check_fit && local_fits(EEG, fs)
    files = fs; info.kind = 'fsaverage'; info.folder = fsDir;
    info.how = 'the fsaverage pial surfaces (no coordinate system given; the contacts fit them)';
    info.message = sprintf(['No coordinate system given and no surfaces for this participant; the contacts fit ' ...
        'fsaverage, so it is drawn (%s).'], fsDir);
else
    info.message = sprintf(['Coordinates in the participant''s own space (%s), and no pial surfaces were found for ' ...
        'this participant (derivatives/freesurfer/sub-<label>). Select them, or build them from the T1 with ' ...
        'FreeSurfer; a template brain will not line up with these coordinates.'], spaceTxt);
end

if o.check_fit && ~isempty(files)
    [info.fit_ok, info.fit] = local_check(EEG, files);
    if ~info.fit_ok
        info.message = sprintf(['%s But %.0f%% of the contacts lie outside these surfaces: the coordinates and ' ...
            'the surfaces are probably not in the same space.'], info.message, 100 * info.fit.frac_outside);
    end
end
end

% ======================= helpers =======================

function [files, folder] = local_participant_surfaces(EEG)
% derivatives/freesurfer/sub-<label> (or freesurfer/sub-<label>) above the data folder
files = {}; folder = '';
p = local_data_folder(EEG);
if ~isempty(p)                                % next to the data (the shipped tutorial)
    g = dir(fullfile(p, '*pial*.gii'));
    if ~isempty(g), files = fullfile({g.folder}, {g.name}); folder = p; return; end
end
sub = local_subject(EEG);
if isempty(sub), return; end
for up = 1:6
    if isempty(p), return; end
    for cand = {fullfile(p, 'derivatives', 'freesurfer', ['sub-' sub]), fullfile(p, 'freesurfer', ['sub-' sub])}
        d = cand{1};
        if ~isfolder(d), continue; end
        g = [dir(fullfile(d, '*pial*.gii')); dir(fullfile(d, 'surf', '*pial*.gii'))];
        if ~isempty(g)
            files = fullfile({g.folder}, {g.name}); folder = d; return
        end
        f = fullfile(d, 'surf', {'lh.pial', 'rh.pial'});
        if all(cellfun(@(x) exist(x, 'file') == 2, f))
            files = f; folder = d; return
        end
    end
    parent = fileparts(p);
    if strcmp(parent, p), return; end
    p = parent;
end
end

function p = local_data_folder(EEG)
p = '';
if isfield(EEG,'filepath') && ~isempty(EEG.filepath), p = EEG.filepath; end
if isempty(p) && isfield(EEG,'BIDS') && isstruct(EEG.BIDS) && isfield(EEG.BIDS,'sourcefile') && ~isempty(EEG.BIDS.sourcefile)
    p = fileparts(char(EEG.BIDS.sourcefile));
end
end

function s = local_subject(EEG)
s = '';
if isfield(EEG,'subject') && ~isempty(EEG.subject), s = char(string(EEG.subject)); end
if isempty(s) && isfield(EEG,'filename') && ~isempty(EEG.filename), s = char(EEG.filename); end
if isempty(s) && isfield(EEG,'BIDS') && isstruct(EEG.BIDS) && isfield(EEG.BIDS,'sourcefile')
    [~, s] = fileparts(char(EEG.BIDS.sourcefile));
end
tok = regexp(s, '^sub-([A-Za-z0-9]+)', 'tokens', 'once');
if ~isempty(tok), s = tok{1}; end
end

function cs = local_coordsystem(EEG)
% The coordsystem.json of the electrodes file that was loaded (same name, BIDS
% pairs them, space- entity included); else one of this participant's next to
% the data, in the ieeg, session or participant folder.
cs = '';
if isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'opt') && isfield(EEG.ieeglab.opt,'elec_tsv') ...
        && ischar(EEG.ieeglab.opt.elec_tsv) && ~isempty(EEG.ieeglab.opt.elec_tsv)
    cs = local_read_cs(regexprep(EEG.ieeglab.opt.elec_tsv, '_electrodes\.tsv$', '_coordsystem.json'));
    if ~isempty(cs), return; end
end
sub = local_subject(EEG);
folders = {local_data_folder(EEG)};
if isfield(EEG,'BIDS') && isstruct(EEG.BIDS) && isfield(EEG.BIDS,'sourcefile') && ~isempty(EEG.BIDS.sourcefile)
    folders{end+1} = fileparts(char(EEG.BIDS.sourcefile));
end
for k = 1:numel(folders)
    p = folders{k};
    for up = 0:2                                  % ieeg folder, session folder, participant folder
        if isempty(p), break; end
        d = dir(fullfile(p, '*_coordsystem.json'));
        if ~isempty(sub)                          % this participant's only, when several are there
            mine = startsWith({d.name}, ['sub-' sub '_']);
            if any(mine), d = d(mine); end
        end
        for i = 1:numel(d)
            cs = local_read_cs(fullfile(d(i).folder, d(i).name));
            if ~isempty(cs), return; end
        end
        p = fileparts(p);
    end
end
end

function cs = local_read_cs(f)
cs = '';
if exist(f, 'file') ~= 2, return; end
try
    j = jsondecode(fileread(f));
    if isfield(j, 'iEEGCoordinateSystem'), cs = strtrim(char(j.iEEGCoordinateSystem)); end
catch
end
end

function xyz = local_xyz(EEG)
xyz = nan(numel(EEG.chanlocs), 3);
if ~isfield(EEG, 'chanlocs') || ~isfield(EEG.chanlocs, 'X'), xyz = zeros(0, 3); return; end
for c = 1:numel(EEG.chanlocs)
    v = [EEG.chanlocs(c).X EEG.chanlocs(c).Y EEG.chanlocs(c).Z];
    if numel(v) == 3, xyz(c, :) = double(v); end
end
xyz = xyz(all(isfinite(xyz), 2), :);
end

function [ok, fit] = local_check(EEG, files)
ok = true; fit = [];
xyz = local_xyz(EEG);
if isempty(xyz), return; end
try
    meshes = cellfun(@ieeglab_read_gifti, files, 'UniformOutput', false);
    [ok, fit] = ieeglab_check_coords(xyz, meshes, 'verbose', false);
catch
end
end

function tf = local_fits(EEG, files)
[tf, fit] = local_check(EEG, files);
tf = tf && ~isempty(fit) && fit.median_dist_mm < 10;
end
