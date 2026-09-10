function EEG = ieeglab_vis_elec(EEG, opt)
% ieeglab_vis_elec
% - If subject pial *.gii surfaces exist in EEG.filepath, uses your original code (unchanged).
% - Else: uses high-res cortex from 'cortex.mat' (variable cortex_highres),
%         auto-aligns electrodes (orientation + rigid PCA + uniform scale),
%         and plots a single rotatable 3D.
% - If cortex.mat is missing, falls back to dipfit/standard_BEM (smooth).


% -------- resolve surface files: explicit opt, then file picker, then folder scan --------
% Issue #2: a selection cached on EEG.ieeglab.opt used to short-circuit the file
% picker, so calling the menu a second time silently reused the first choice with
% no way to change it. The picker now always opens in the interactive path; the
% cached list is only a starting directory hint. Pass opt.surf_files to skip it.
if nargin < 2, opt = struct(); end

surf_files = {};   % full paths
if isfield(opt,'surf_files') && ~isempty(opt.surf_files)
    sf = opt.surf_files;
    if ischar(sf) || isstring(sf), sf = cellstr(sf); end
    surf_files = local_resolve_paths(sf, EEG.filepath);

else
    start_dir = EEG.filepath;
    if isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'opt') && ...
       isfield(EEG.ieeglab.opt,'surf_files') && ~isempty(EEG.ieeglab.opt.surf_files)
        prev = cellstr(EEG.ieeglab.opt.surf_files);
        p = fileparts(char(prev{1}));
        if ~isempty(p) && isfolder(p), start_dir = p; end
    end
    if isempty(start_dir) || ~isfolder(start_dir), start_dir = pwd; end

    [fn, fp] = uigetfile( ...
        {'*.gii;*.stl;*.obj;*.ply','Surface files (*.gii, *.stl, *.obj, *.ply)'; ...
         '*.*','All files (*.*)'}, ...
        'Select anatomical surface file(s) - you may select several (e.g. LH and RH)', ...
        start_dir, 'MultiSelect','on');

    if isequal(fn,0)
        % Cancelled -> scan the dataset folder for pial or white surfaces
        d = dir(fullfile(EEG.filepath, '*.gii'));
        names = {d.name}';
        hit = names(contains(names,'pial') | contains(names,'white'));
        surf_files = cellfun(@(x) fullfile(EEG.filepath, x), hit, 'UniformOutput', false);
        if ~isempty(surf_files)
            fprintf('[vis_elec] No file selected; using %d surface(s) found in %s.\n', ...
                numel(surf_files), EEG.filepath);
        end
    else
        if ischar(fn) || isstring(fn), fn = cellstr(fn); end
        % Keep FULL paths: the directory the user browsed to was previously
        % discarded, so any mesh outside EEG.filepath failed to load.
        surf_files = cellfun(@(x) fullfile(fp, x), fn, 'UniformOutput', false);
        EEG.ieeglab.opt.surf_files = surf_files;
    end
end

figure('color','w'); hold on
try icadefs; set(gcf, 'color', BACKCOLOR); catch, end  % eeglab color

% Subject pial-surface plotting from vistasoft
if ~isempty(surf_files)        
        
    try icadefs; set(gcf, 'color', BACKCOLOR); catch; end  % eeglab color

    assert(exist('gifti','file') ~= 0, ...
        ['The gifti toolbox is required to read .gii surfaces but is not on the MATLAB path.\n' ...
         'Install vistasoft (https://github.com/vistalab/vistasoft) and add it with addpath(genpath(...)), ' ...
         'or run ieeglab_check_install for details.']);

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
    for iSurf = 1:numel(surf_files)
        f = surf_files{iSurf};
        if exist(f,'file') ~= 2
            warning('ieeglab_vis_elec:missingSurface', 'Surface file not found, skipping: %s', f);
            continue
        end
        try
            g  = gifti(f);
            tH = ieeg_RenderGifti(g);
            tH.FaceAlpha = 0.1;
            nRendered = nRendered + 1;
        catch ME
            warning('ieeglab_vis_elec:renderFailed', ...
                'Could not render %s: %s', f, ME.message);
        end
    end
    if nRendered == 0
        error('ieeglab_vis_elec:noSurfaceRendered', ...
            'None of the %d selected surface file(s) could be rendered. See the warnings above.', numel(surf_files));
    end

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
    % Fallback if no Freesurfer file is available: dipfit standard_BEM (smooth)
    dipfit_root = fileparts(which('dipfitdefs'));
    assert(~isempty(dipfit_root), 'dipfit not found on path. Enable dipfit in EEGLAB.');
    std_dir = fullfile(dipfit_root, 'standard_BEM');
    f_sccn  = fullfile(std_dir, 'standard_vol_SCCN.mat');
    f_std   = fullfile(std_dir, 'standard_vol.mat');

    if exist(f_sccn,'file'), S = load(f_sccn); else, S = load(f_std); end
    assert(isfield(S,'vol'), 'standard_BEM volume not found in dipfit.');

    bnd = S.vol.bnd; if numel(bnd) > 1, bnd = bnd(end); end
    if isfield(bnd,'pnt'), Vb = double(bnd.pnt); else, Vb = double(bnd.pos); end
    if isfield(bnd,'tri'), Fb = double(bnd.tri); else, Fb = double(bnd.face); end

    patch('Faces',Fb,'Vertices',Vb, ...
        'FaceColor',[0.75 0.80 0.90],'EdgeColor','none','FaceAlpha',0.12);

    E = [[EEG.chanlocs.X]' [EEG.chanlocs.Y]' [EEG.chanlocs.Z]'];
    if median(abs(E(:))) < 2, E = E*1000; end  % units only
    s = scatter3(E(:,1), E(:,2), E(:,3), 'o', 'filled');
    s.SizeData = 10; s.MarkerFaceColor = [.9 .5 .5]; s.MarkerEdgeColor = [0 0 0];

    try icadefs; set(gcf, 'color', BACKCOLOR); catch; end  % eeglab color
    axis equal off 
    % vis3d tight
    camlight headlight; camlight right
    lighting gouraud; material dull
    view([-135 20]);
    title('Visualization using standard BEM template (fallback when no Freesurfer pial surface files is detected)','Interpreter','none');

end

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
