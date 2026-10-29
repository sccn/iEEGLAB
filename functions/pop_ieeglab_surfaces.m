function [EEG, com] = pop_ieeglab_surfaces(EEG, varargin)
% pop_ieeglab_surfaces() - Choose, once, the brain every brain figure is drawn on.
%
% Usage:
%   EEG = pop_ieeglab_surfaces(EEG)                                     % dialog
%   EEG = pop_ieeglab_surfaces(EEG, 'surf_files', {'pial.L.surf.gii', 'pial.R.surf.gii'})
%   EEG = pop_ieeglab_surfaces(EEG, 'brain', 'template')                % template brain
%   EEG = pop_ieeglab_surfaces(EEG, 'brain', 'auto')                    % forget the choice
%
% Menu: iEEGLAB > Choose brain surfaces, before Visualize electrodes. Visualize
% electrodes, Plot electrode values on brain and the band-power brain map then
% draw on this brain and never ask. Without a choice they draw on the brain
% that matches the space of the coordinates (ieeglab_brain_surfaces): the
% participant's pial surfaces, or fsaverage's for template coordinates; failing
% both, the template brain.
%
% Inputs:
%   'surf_files'  pial surface files, one per hemisphere: GIfTI (pial.L.surf.gii)
%                 or FreeSurfer (surf/lh.pial). They must be in the same space as
%                 the electrode coordinates. Bare names are looked for next to
%                 the dataset.
%   'brain'       'template' (the template brain), or 'auto' (no choice: the
%                 brain that matches the coordinates)
%
% The surfaces are stored on the dataset (EEG.ieeglab.opt.surf_files) and
% remembered for the same participant of the same dataset in later MATLAB
% sessions (ieeglab_remembered_surfaces). The template brain is stored on the
% dataset (EEG.ieeglab.opt.brain = 'template').
%
% Outputs:
%   EEG  the dataset with the choice stored
%   com  the command that repeats it ('' when the dialog was closed)
%
% Cedric Cannard, iEEGLAB, 2026

com = '';
ieeglab_require_data(EEG, 'pop_ieeglab_surfaces');
if mod(numel(varargin), 2) ~= 0
    error('pop_ieeglab_surfaces:badArgs', 'Options come in name-value pairs.');
end
files = {}; brain = '';
for k = 1:2:numel(varargin)
    switch lower(char(varargin{k}))
        case 'surf_files', files = varargin{k+1};
        case 'brain',      brain = lower(char(varargin{k+1}));
        otherwise
            error('pop_ieeglab_surfaces:unknownOption', 'Unknown option ''%s''.', char(varargin{k}));
    end
end
if ischar(files) || isstring(files), files = cellstr(files); end
if ~isempty(brain) && ~any(strcmp(brain, {'template', 'auto'}))
    error('pop_ieeglab_surfaces:brain', '''brain'' is ''template'' or ''auto'', not ''%s''.', brain);
end
if ~isempty(brain) && ~isempty(files)
    error('pop_ieeglab_surfaces:both', 'Give ''surf_files'' or ''brain'', not both.');
end

if isempty(files) && isempty(brain)
    [files, brain] = local_dialog(EEG);
    if isempty(files) && isempty(brain), return; end   % closed: dataset unchanged
end

if ~isfield(EEG, 'ieeglab') || ~isstruct(EEG.ieeglab), EEG.ieeglab = struct(); end
if ~isfield(EEG.ieeglab, 'opt') || ~isstruct(EEG.ieeglab.opt), EEG.ieeglab.opt = struct(); end
o = rmfield(EEG.ieeglab.opt, intersect(fieldnames(EEG.ieeglab.opt), {'surf_files', 'brain'}));
switch brain
    case 'template'
        o.brain = 'template';
        com = 'EEG = pop_ieeglab_surfaces(EEG, ''brain'', ''template'');';
        fprintf('[surfaces] Every brain figure of this dataset is drawn on the template brain.\n');
    case 'auto'
        ieeglab_remembered_surfaces(EEG, {});          % forgotten for this participant too
        com = 'EEG = pop_ieeglab_surfaces(EEG, ''brain'', ''auto'');';
        [~, bi] = ieeglab_brain_surfaces(EEG, 'check_fit', false);
        fprintf('[surfaces] No choice: %s\n', bi.message);
    otherwise
        files = local_resolve_paths(files, EEG.filepath);
        if isempty(files)
            error('pop_ieeglab_surfaces:notFound', 'None of the surface files was found (see the warnings above).');
        end
        o.surf_files = files;
        ieeglab_remembered_surfaces(EEG, files);
        com = sprintf('EEG = pop_ieeglab_surfaces(EEG, ''surf_files'', {%s});', ...
            strjoin(strcat('''', strrep(files, '''', ''''''), ''''), ', '));
        [~, nm, ext] = cellfun(@fileparts, files, 'UniformOutput', false);
        fprintf('[surfaces] Every brain figure of this dataset is drawn on %s.\n', strjoin(strcat(nm, ext), ', '));
end
EEG.ieeglab.opt = o;
end

function [files, brain] = local_dialog(EEG)
% The current choice, else the brain that matches the coordinates, is offered;
% other pial files or the template brain can be picked instead.
files = {}; brain = '';
prev = ieeglab_remembered_surfaces(EEG);
startDir = local_start_dir(EEG, prev);
[best, bi] = ieeglab_brain_surfaces(EEG);
offered = prev; how = 'the brain surfaces chosen earlier';
if isempty(offered), offered = best; how = bi.how; end
% Instructions go to the command window first (file-window titles are not
% shown on macOS), then into the dialogs.
bar = repmat('=', 1, 78);
fprintf('\n%s\n iEEGLAB > Choose brain surfaces: the brain every brain figure is drawn on\n%s\n', bar, bar);
fprintf(' BRAIN: %s\n', bi.message);
fprintf([' TO SELECT OTHERS: pial surface files, one per hemisphere: GIfTI (pial.L.surf.gii,\n' ...
         '   pial.R.surf.gii) or FreeSurfer (surf/lh.pial, surf/rh.pial); Ctrl+click (Cmd+click on\n' ...
         '   a Mac) selects both. In a BIDS dataset: derivatives/freesurfer/sub-<label>. They must\n' ...
         '   be in the same space as the electrode coordinates. The file window opens in:\n' ...
         '   %s\n'], startDir);
if ~isempty(offered)
    [~, nm, ext] = cellfun(@fileparts, offered, 'UniformOutput', false);
    fprintf(' FOUND (%s):\n%s', how, sprintf('   %s\n', offered{:}));
    fprintf('%s\n\n', bar);
    choice = questdlg(sprintf('Draw every brain figure on %s?\n\n%s\n\nin %s', how, ...
        strjoin(strcat(nm, ext), newline), fileparts(offered{1})), 'iEEGLAB: brain surfaces', ...
        'Use these', 'Select other files', 'Template brain', 'Use these');
else
    fprintf('%s\n\n', bar);
    dflt = 'Select pial files';
    if bi.template_ok, dflt = 'Template brain'; end
    choice = questdlg(sprintf(['%s\n\nSelect pial files: one per hemisphere, GIfTI (pial.L.surf.gii) or ' ...
        'FreeSurfer (lh.pial); Ctrl+click (Cmd+click on a Mac) selects both.'], bi.message), ...
        'iEEGLAB: brain surfaces', 'Select pial files', 'Template brain', dflt);
end
switch choice
    case 'Use these'
        files = offered;
    case 'Template brain'
        brain = 'template';
    case {'Select other files', 'Select pial files'}
        fprintf('[surfaces] File window open: select the PIAL files (.gii), both hemispheres (Ctrl/Cmd+click).\n');
        [fn, fp] = uigetfile( ...
            {'*.gii;*.pial;*.white', 'Surface files (GIfTI *.gii, FreeSurfer lh.pial / rh.pial)'; ...
             '*.*', 'All files (*.*)'}, ...
            'SELECT PIAL SURFACES (.gii): pial.L and pial.R, Ctrl/Cmd+click for both', ...
            startDir, 'MultiSelect', 'on');
        if ~isequal(fn, 0)
            % Full paths: meshes outside EEG.filepath must still load
            files = fullfile(fp, cellstr(fn));
        end
end
end

function s = local_subject(EEG)
s = '';
if isfield(EEG,'subject') && ~isempty(EEG.subject), s = char(string(EEG.subject)); end
if isempty(s), s = char(EEG.filename); end   % BIDS file name: sub-<label>_...
tok = regexp(s, '^sub-([A-Za-z0-9]+)', 'tokens', 'once');
if ~isempty(tok), s = tok{1}; end
end

function d = local_start_dir(EEG, prev)
% Where the file picker opens: last folder used, else the subject's
% derivatives/freesurfer folder of the BIDS dataset, else the dataset folder
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

function out = local_resolve_paths(sf, filepath)
% Bare file names or full paths, resolved against the dataset folder
out = {};
for k = 1:numel(sf)
    p = char(sf{k});
    if exist(p, 'file') == 2
        out{end+1} = p; %#ok<AGROW>
        continue
    end
    [~, nm, ext] = fileparts(p);
    cand = fullfile(filepath, [nm ext]);
    if exist(cand, 'file') == 2
        out{end+1} = cand; %#ok<AGROW>
    else
        warning('pop_ieeglab_surfaces:surfaceNotFound', ...
            'Surface file not found as "%s" nor in the dataset folder as "%s".', p, cand);
    end
end
end
