function EEG = ieeglab_vis_elec(EEG, opt)
% ieeglab_vis_elec() - Draw the electrodes in 3D, on the brain surfaces if given.
%
% Usage:
%   EEG = ieeglab_vis_elec(EEG)                                    % asks for the surfaces once
%   EEG = ieeglab_vis_elec(EEG, struct('reselect', true))          % choose other surfaces
%   EEG = ieeglab_vis_elec(EEG, struct('surf_files', {{'pial.L.surf.gii','pial.R.surf.gii'}}))
%
% - Surfaces (.gii) given in opt.surf_files or picked in the dialog are drawn
%   under the electrodes; they must be in the same space as the coordinates.
%   Once chosen they are remembered and used by every brain figure (this one
%   and ieeglab_topoplot) without asking again.
% - Without surfaces, a template cortex ('cortex.mat', variable cortex_highres)
%   is used, with the electrodes aligned to it (orientation, rigid PCA, uniform
%   scale); without it, the dipfit standard BEM head.


% -------- resolve surface files: explicit opt, then remembered files, then file picker --------
% The surfaces chosen once are remembered (on the dataset, and for the same
% subject across datasets and sessions). The next call asks whether to reuse
% them, so they can still be changed (issue #2: a silent reuse left no way to).
% Pass opt.surf_files to skip all dialogs.
ieeglab_require_data(EEG, 'ieeglab_vis_elec');
if nargin < 2, opt = struct(); end

surf_files = {};   % full paths
if isfield(opt,'surf_files')
    sf = opt.surf_files;
    if ischar(sf) || isstring(sf), sf = cellstr(sf); end
    if ~isempty(sf), surf_files = local_resolve_paths(sf, EEG.filepath); end

elseif ~(isfield(opt,'reselect') && isequal(opt.reselect, true)) && ~isempty(ieeglab_remembered_surfaces(EEG))
    % Asked once: every brain figure then uses the same surfaces. They are
    % changed with opt.reselect (menu: iEEGLAB > Choose brain surfaces).
    surf_files = ieeglab_remembered_surfaces(EEG);
    [~, nm, ext] = cellfun(@fileparts, surf_files, 'UniformOutput', false);
    fprintf('[vis_elec] Using the brain surfaces chosen earlier: %s\n           To choose others: iEEGLAB > Choose brain surfaces.\n', ...
        strjoin(strcat(nm, ext), ', '));

else
    % Instructions go to the command window first (file-window titles are not
    % shown on macOS), then into the dialogs.
    prev = ieeglab_remembered_surfaces(EEG);
    startDir = local_start_dir(EEG, prev);
    offered = prev; how = 'the brain surfaces used last time';
    if isempty(offered)                      % none remembered: look for the subject's pial files
        d = dir(fullfile(startDir, '*pial*.gii'));
        offered = cellfun(@(n) fullfile(startDir, n), {d.name}, 'UniformOutput', false);
        how = 'pial surface files found for this subject';
    end
    bar = repmat('=', 1, 78);
    fprintf('\n%s\n iEEGLAB > Visualize electrodes: which brain to draw the electrodes on\n%s\n', bar, bar);
    fprintf([' WHAT TO SELECT: the subject''s PIAL SURFACE files (GIfTI, .gii), one per hemisphere,\n' ...
             '   for example  pial.L.surf.gii  and  pial.R.surf.gii\n' ...
             '   (select both at once: Ctrl+click on Windows/Linux, Cmd+click on a Mac).\n' ...
             ' WHERE: in a BIDS dataset, derivatives/freesurfer/sub-<label>. The file window opens in:\n' ...
             '   %s\n' ...
             ' The surfaces must be in the same space as the electrode coordinates.\n' ...
             ' No surfaces for this subject? Choose "Template brain" (or Cancel in the file window).\n'], startDir);
    if ~isempty(offered)
        [~, nm, ext] = cellfun(@fileparts, offered, 'UniformOutput', false);
        fprintf(' FOUND (%s):\n%s', how, sprintf('   %s\n', offered{:}));
        fprintf('%s\n\n', bar);
        choice = questdlg(sprintf('Draw the electrodes on %s?\n\n%s\n\nin %s', how, ...
            strjoin(strcat(nm, ext), newline), fileparts(offered{1})), 'iEEGLAB: brain surfaces', ...
            'Use these', 'Select other files', 'Template brain', 'Use these');
    else
        fprintf('%s\n\n', bar);
        choice = questdlg(sprintf(['WHAT TO SELECT in the next window: the subject''s pial surface files ' ...
            '(GIfTI, .gii), one per hemisphere, e.g. pial.L.surf.gii and pial.R.surf.gii. Select both at ' ...
            'once with Ctrl+click (Cmd+click on a Mac).\n\nWHERE: in a BIDS dataset, ' ...
            'derivatives/freesurfer/sub-<label>.\n\nThe surfaces must be in the same space as the ' ...
            'electrode coordinates. Without subject surfaces, choose Template brain.']), ...
            'iEEGLAB: brain surfaces', 'Select pial files', 'Template brain', 'Select pial files');
    end
    if isempty(choice), return; end   % dialog closed: nothing drawn, dataset unchanged

    switch choice
        case 'Use these'
            surf_files = offered;
        case 'Template brain'
            surf_files = {};
        otherwise
            fprintf('[vis_elec] File window open: select the PIAL files (.gii), both hemispheres (Ctrl/Cmd+click).\n');
            [fn, fp] = uigetfile( ...
                {'*.gii;*.stl;*.obj;*.ply','Surface files (*.gii, *.stl, *.obj, *.ply)'; ...
                 '*.*','All files (*.*)'}, ...
                'SELECT PIAL SURFACES (.gii): pial.L and pial.R, Ctrl/Cmd+click for both', ...
                startDir, 'MultiSelect','on');
            if ~isequal(fn,0)
                if ischar(fn) || isstring(fn), fn = cellstr(fn); end
                % Keep FULL paths: meshes outside EEG.filepath must still load.
                surf_files = cellfun(@(x) fullfile(fp, x), fn, 'UniformOutput', false);
            end
    end
    if ~isempty(surf_files), ieeglab_remembered_surfaces(EEG, surf_files); end
end
if ~isempty(surf_files), EEG.ieeglab.opt.surf_files = surf_files; end

ieeglab_figure('Name','iEEGLAB electrodes'); hold on
try icadefs; set(gcf, 'color', BACKCOLOR); catch, end  % eeglab color

% Subject pial-surface plotting (GIfTI read by ieeglab_read_gifti)
if ~isempty(surf_files)        
        
    try icadefs; set(gcf, 'color', BACKCOLOR); catch; end  % eeglab color

    % Electrode coordinates, keeping the mapping back to EEG.chanlocs intact.
    % [EEG.chanlocs.X] silently drops channels with an empty X, which used to
    % make the hemisphere mask shorter than chanlocs and mis-index the scatter.
    hasXYZ = arrayfun(@(c) ~isempty(c.X) && ~isempty(c.Y) && ~isempty(c.Z) && ...
                           all(isfinite([c.X c.Y c.Z])), EEG.chanlocs);
    XYZ = nan(numel(EEG.chanlocs), 3);
    XYZ(hasXYZ,:) = cell2mat(arrayfun(@(c) [c.X c.Y c.Z], EEG.chanlocs(hasXYZ), 'UniformOutput', false)');
    if any(~hasXYZ)
        fprintf('[vis_elec] %d/%d channel(s) have no 3D coordinates and are not plotted.\n', ...
            nnz(~hasXYZ), numel(hasXYZ));
    end

    % Issue #3: render EVERY selected mesh rather than picking one per
    % hemisphere. surf_files{logicalMask} with more than one match assigned a
    % comma-separated list to a scalar, which was a hard error.
    nRendered = 0;
    meshes = {};
    for iSurf = 1:numel(surf_files)
        f = surf_files{iSurf};
        if exist(f,'file') ~= 2
            warning('ieeglab_vis_elec:missingSurface', 'Surface file not found, skipping: %s', f);
            continue
        end
        try
            g  = ieeglab_read_gifti(f);
            tH = ieeg_RenderGifti(g, ieeglab_surface_sulc(g));   % sulci darker than gyri
            tH.FaceAlpha = 0.15;
            nRendered = nRendered + 1;
            meshes{end+1} = g; %#ok<AGROW>
        catch ME
            warning('ieeglab_vis_elec:renderFailed', ...
                'Could not render %s: %s', f, ME.message);
        end
    end
    if nRendered == 0
        error('ieeglab_vis_elec:noSurfaceRendered', ...
            'None of the %d selected surface file(s) could be rendered. See the warnings above.', numel(surf_files));
    end
    % Warn when the contacts do not sit inside the surfaces (different spaces)
    ieeglab_check_coords(XYZ(hasXYZ,:), meshes);

    % Plot all electrodes that have coordinates, in one call. Hemisphere is no
    % longer inferred from sign(X) - a depth lead crossing the midline broke
    % that assumption, and with every mesh drawn there is nothing to gate on.
    show = hasXYZ(:)';
    s = scatter3(XYZ(show,1), XYZ(show,2), XYZ(show,3), 'o', 'Filled');
    s.SizeData = 20;
    s.MarkerFaceColor = [.9 .5 .5];
    s.MarkerEdgeColor = [0 0 0];
    ieeg_viewLight(90,0)

    axis equal off
    [~, shown] = cellfun(@(p) fileparts(p), surf_files, 'UniformOutput', false);
    title(sprintf('%d electrode(s) on %d surface(s): %s', ...
        nnz(show), nRendered, strjoin(shown, ', ')), 'Interpreter','none');

else
    % No subject surfaces: dipfit standard BEM brain (template, MNI space)
    E = nan(numel(EEG.chanlocs), 3);   % channels without coordinates stay NaN (not drawn)
    for c = 1:numel(EEG.chanlocs)
        v = [EEG.chanlocs(c).X EEG.chanlocs(c).Y EEG.chanlocs(c).Z];
        if numel(v) == 3, E(c,:) = double(v); end
    end
    E = E * ieeglab_draw_template_brain(E);
    s = scatter3(E(:,1), E(:,2), E(:,3), 'o', 'filled');
    s.SizeData = 10; s.MarkerFaceColor = [.9 .5 .5]; s.MarkerEdgeColor = [0 0 0];
    try icadefs; set(gcf, 'color', BACKCOLOR); catch; end  % eeglab color
    axis equal off
    view([-135 20]);
    title('Visualization using standard BEM template (fallback when no Freesurfer pial surface files is detected)','Interpreter','none');

end
ieeglab_rotate3d(gcf);

end

function s = local_subject(EEG)
    s = '';
    if isfield(EEG,'subject') && ~isempty(EEG.subject), s = char(string(EEG.subject)); end
    if isempty(s), s = char(EEG.filename); end   % BIDS file name: sub-<label>_...
    tok = regexp(s, '^sub-([A-Za-z0-9]+)', 'tokens', 'once');
    if ~isempty(tok), s = tok{1}; end
end

% --- helper: where the file picker opens: last folder used, else the subject's
% derivatives/freesurfer folder of the BIDS dataset, else the dataset folder ---
function d = local_start_dir(EEG, prev)
    d = '';
    if ~isempty(prev), d = fileparts(prev{1}); end
    p = EEG.filepath;
    if isempty(p) && isfield(EEG,'BIDS') && isstruct(EEG.BIDS) && isfield(EEG.BIDS,'sourcefile') && ~isempty(EEG.BIDS.sourcefile)
        p = fileparts(char(EEG.BIDS.sourcefile));     % dataset read from a raw BIDS file, not saved yet
    end
    home = p;
    sub = ['sub-' local_subject(EEG)];
    for up = 1:6
        if ~isempty(d) || isempty(p), break; end
        for cand = {fullfile(p, 'derivatives', 'freesurfer', sub), fullfile(p, 'freesurfer', sub)}
            if isfolder(cand{1}), d = cand{1}; break; end
        end
        parent = fileparts(p);
        if strcmp(parent, p), break; end
        p = parent;
    end
    if isempty(d) || ~isfolder(d), d = home; end
    if isempty(d) || ~isfolder(d), d = pwd; end
end

% --- helper: return just the filename from a path or name ---
function nm = get_name_only(p)
    [~,nm,ext] = fileparts(char(p));
    nm = [nm ext];
end

% --- helper: accept bare filenames or full paths, resolve against the dataset ---
function out = local_resolve_paths(sf, filepath)
    out = {};
    for k = 1:numel(sf)
        p = char(sf{k});
        if exist(p,'file') == 2
            out{end+1} = p; %#ok<AGROW>
            continue
        end
        cand = fullfile(filepath, get_name_only(p));
        if exist(cand,'file') == 2
            out{end+1} = cand; %#ok<AGROW>
        else
            warning('ieeglab_vis_elec:surfaceNotFound', ...
                'Surface file not found as "%s" nor in the dataset folder as "%s".', p, cand);
        end
    end
end
