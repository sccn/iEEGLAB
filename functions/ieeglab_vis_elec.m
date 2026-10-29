function EEG = ieeglab_vis_elec(EEG, opt)
% ieeglab_vis_elec() - Draw the electrodes in 3D, on the brain surfaces if given.
%
% Usage:
%   EEG = ieeglab_vis_elec(EEG)                                    % the brain chosen once, never asks
%   EEG = ieeglab_vis_elec(EEG, struct('surf_files', {{'pial.L.surf.gii','pial.R.surf.gii'}}))
%   EEG = ieeglab_vis_elec(EEG, struct('reselect', true))          % choose the brain first
%
% - The brain is chosen once with iEEGLAB > Choose brain surfaces
%   (pop_ieeglab_surfaces) and used by every brain figure (this one and
%   ieeglab_topoplot) without asking. Without a choice, the brain that matches
%   the space of the coordinates (ieeglab_brain_surfaces): the participant's
%   pial surfaces, or fsaverage's for template coordinates when FreeSurfer or
%   MNE's fsaverage is installed.
% - opt.surf_files (GIfTI .gii, or FreeSurfer files such as lh.pial) are drawn
%   instead; they must be in the same space as the coordinates.
% - Without surfaces, a template cortex ('cortex.mat', variable cortex_highres)
%   is used, with the electrodes aligned to it (orientation, rigid PCA, uniform
%   scale); without it, the dipfit standard BEM head.


% -------- which brain: opt.surf_files, else the one chosen, else the one matching the coordinates --------
ieeglab_require_data(EEG, 'ieeglab_vis_elec');
if nargin < 2, opt = struct(); end
if isfield(opt,'reselect') && isequal(opt.reselect, true)
    [EEG, com] = pop_ieeglab_surfaces(EEG);
    if isempty(com), return; end   % dialog closed: nothing drawn, dataset unchanged
end

surf_files = {};   % full paths
if isfield(opt,'surf_files')
    sf = opt.surf_files;
    if ischar(sf) || isstring(sf), sf = cellstr(sf); end
    if ~isempty(sf), surf_files = local_resolve_paths(sf, EEG.filepath); end
    if ~isempty(surf_files)
        EEG.ieeglab.opt.surf_files = surf_files;
        if isfield(EEG.ieeglab.opt, 'brain'), EEG.ieeglab.opt = rmfield(EEG.ieeglab.opt, 'brain'); end
    end

elseif isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'opt') && isfield(EEG.ieeglab.opt,'brain') && ...
        strcmp(EEG.ieeglab.opt.brain, 'template')
    fprintf('[vis_elec] Template brain, as chosen in iEEGLAB > Choose brain surfaces.\n');

elseif ~isempty(ieeglab_remembered_surfaces(EEG))
    surf_files = ieeglab_remembered_surfaces(EEG);
    [~, nm, ext] = cellfun(@fileparts, surf_files, 'UniformOutput', false);
    fprintf('[vis_elec] Brain surfaces chosen earlier: %s\n           To choose others: iEEGLAB > Choose brain surfaces.\n', ...
        strjoin(strcat(nm, ext), ', '));

else
    [surf_files, bi] = ieeglab_brain_surfaces(EEG);
    fprintf('[vis_elec] %s\n           To choose others: iEEGLAB > Choose brain surfaces.\n', bi.message);
end

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
    [az, el] = ieeglab_brain_view(EEG, XYZ(show,:));   % sEEG: front, slightly above; ECoG: face-on
    ieeg_viewLight(az, el)

    axis equal off
    shown = cellfun(@get_name_only, surf_files, 'UniformOutput', false);
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
    [az, el] = ieeglab_brain_view(EEG, E); view([az el]);   % sEEG: front, slightly above; ECoG: face-on
    [~, bi] = ieeglab_brain_surfaces(EEG, 'check_fit', false);
    if bi.template_ok
        title({'Template brain (MNI): an approximation for these template coordinates', ...
               'fsaverage''s pial surfaces are drawn instead when FreeSurfer or MNE''s fsaverage is installed'}, ...
              'Interpreter', 'none', 'FontWeight', 'normal');
    else
        title({'Template brain: the coordinates are in the participant''s space and will not line up with it', ...
               'Give the participant''s pial surfaces: iEEGLAB > Choose brain surfaces'}, ...
              'Interpreter', 'none', 'FontWeight', 'normal');
    end

end
ieeglab_rotate3d(gcf);

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
