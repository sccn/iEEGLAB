function [vals, info] = ieeglab_topoplot(EEG, what, varargin)
% ieeglab_topoplot() - Plot one value per electrode on the brain: the iEEG
%                      counterpart of EEGLAB's topoplot.
%
% Usage:
%   ieeglab_topoplot(EEG, 120)                        % trial-mean amplitude at 120 ms
%   ieeglab_topoplot(EEG, [100 150])                  % mean over 100-150 ms
%   ieeglab_topoplot(EEG, values)                     % any per-channel vector (nbchan x 1)
%   ieeglab_topoplot(EEG, 'in_degree')                % from the CCEP connectivity matrix
%   ieeglab_topoplot(EEG, 'n1_amplitude', 'site','RA1-RA2')   % responses to one site
%   ieeglab_topoplot(EEG, 'n1_latency',   'site','RA1-RA2')
%   [vals, info] = ieeglab_topoplot(..., 'draw', false)       % compute only
%
% Options:
%   'draw'        default true
%   'site'        stimulation site. Required for the n1_* / crp_* values; for a
%                 latency or window it restricts the average to that site's trials
%   'sig_only'    per-site metrics: show only significant responses; tested pairs
%                 without one are drawn grey like untested ones. Default true, as
%                 in ieeglab_plot_ccep_matrix
%   'surf_files'  surfaces to draw under the electrodes (as in ieeglab_vis_elec);
%                 defaults to the ones chosen earlier, stored in EEG.ieeglab.opt
%   'clim'        colour limits
%   'title'       figure title
%
% Why not topoplot: EEGLAB's topoplot interpolates values across a 2D scalp.
% Depth and grid electrodes are not on a scalp, and interpolating between a
% contact in the hippocampus and one in white matter produces a picture of
% nothing. The honest equivalent is the electrodes themselves, in 3D,
% coloured by the value. ieeglab_is_ieeg() identifies datasets that should be
% routed here instead of topoplot.
%
% Channels without coordinates are skipped; channels with a NaN value are
% drawn small and grey, so "no value" stays distinct from "value near zero".
%
% Latency and window maps on CCEP data leave out, for every contact, the trials
% in which that contact was itself stimulated - its own stimulation artifact
% would otherwise dominate the map (N1 and CRP exclude them the same way).
%
% Cedric Cannard, iEEGLAB, 2026

p = inputParser;
p.addParameter('draw', true);
p.addParameter('site', '');
p.addParameter('surf_files', {});
p.addParameter('clim', []);
p.addParameter('title', '');
p.addParameter('sig_only', true);
p.parse(varargin{:});
o = p.Results;

C = EEG.nbchan;
labels = {EEG.chanlocs.labels};
usedSite = false;   % whether the values are specific to o.site (for the title)
ownStim = false;    % whether own-stimulation trials were left out

% ---------- values ----------
if isnumeric(what) && numel(what) == C && C > 2
    vals = double(what(:));
    desc = 'value';
elseif isnumeric(what) && any(numel(what) == [1 2])
    if ~isfield(EEG,'times') || numel(EEG.times) ~= size(EEG.data,2)
        error('ieeglab_topoplot:noTimes', 'A latency needs EEG.times matching the data.');
    end
    w = double(what(:))';
    t0 = double(EEG.times(1)); t1 = double(EEG.times(end));
    if any(w < t0 | w > t1)
        error('ieeglab_topoplot:latencyOutside', ...
            'Latency %s ms is outside the epoch [%g %g] ms.', mat2str(w), t0, t1);
    end
    if numel(w) == 1
        [~, idx] = min(abs(double(EEG.times) - w));
        desc = sprintf('Amplitude at %g ms (\\muV)', EEG.times(idx));
    else
        idx = find(EEG.times >= min(w) & EEG.times <= max(w));
        if isempty(idx)
            error('ieeglab_topoplot:emptyWindow', 'No sample between %g and %g ms at %g Hz.', min(w), max(w), EEG.srate);
        end
        desc = sprintf('Mean amplitude %g-%g ms (\\muV)', min(w), max(w));
    end
    tr = 1:size(EEG.data,3);
    if ~isempty(o.site)
        if EEG.trials < 2
            error('ieeglab_topoplot:notEpoched', '''site'' needs epoched data.');
        end
        sites = ieeglab_epoch_sites(EEG);
        cs = ieeglab_canonical_site({char(o.site)}, labels);
        tr = find(sites == cs(1));
        if isempty(tr)
            error('ieeglab_topoplot:unknownSite', 'No epoch was stimulated at "%s". Sites: %s', ...
                char(o.site), strjoin(cellstr(unique(sites)), ', '));
        end
        usedSite = true;
    end
    D = double(EEG.data(:, idx, tr));
    if EEG.trials > 1 && strcmp(ieeglab_detect_mode(EEG), 'ccep')
        [~, stimIdx] = ieeglab_epoch_sites(EEG);
        for k = 1:numel(tr)
            sIdx = stimIdx{tr(k)};
            if ~isempty(sIdx), D(sIdx, :, k) = NaN; end
        end
        ownStim = true;
    end
    vals = reshape(mean(mean(D, 3, 'omitnan'), 2, 'omitnan'), [], 1);
elseif ischar(what) || isstring(what)
    [vals, desc] = local_named_value(EEG, char(what), o.site, o.sig_only);
    usedSite = ~isempty(o.site) && ~any(strcmpi(char(what), {'in_degree','indegree'}));
else
    error('ieeglab_topoplot:badInput', ...
        ['Second argument must be a latency in ms, a [start stop] window, a vector with one ' ...
         'value per channel (%d), or a named metric such as ''in_degree''.'], C);
end
vals = vals(:);

% ---------- coordinates ----------
xyz = nan(C, 3);
for c = 1:C
    cl = EEG.chanlocs(c);
    if isfield(cl,'X') && ~isempty(cl.X) && ~isempty(cl.Y) && ~isempty(cl.Z)
        xyz(c,:) = double([cl.X cl.Y cl.Z]);
    end
end
hasXYZ = all(isfinite(xyz), 2);
if ~any(hasXYZ)
    error('ieeglab_topoplot:noCoordinates', ...
        'No channel has 3D coordinates. Load an electrodes.tsv first (iEEGLAB > Load).');
end

info = struct('values', vals, 'xyz', xyz, 'labels', {labels}, 'has_xyz', hasXYZ, ...
    'has_value', isfinite(vals), 'description', desc);
if isempty(o.clim)
    v = vals(hasXYZ & isfinite(vals));
    if isempty(v), info.clim = [0 1];
    elseif all(v >= 0), info.clim = [0 max(v)];
    else, m = max(abs(v)); info.clim = [-m m];
    end
    if info.clim(1) == info.clim(2), info.clim = info.clim + [-1 1]; end
else
    info.clim = o.clim;
end
info.title = o.title;
if isempty(info.title)
    info.title = desc;
    if usedSite, info.title = sprintf('%s - stimulation of %s', desc, char(o.site)); end
    if ownStim, info.title = {info.title, 'each contact''s own stimulation trials excluded'}; end
end
if ~o.draw, return; end

% ---------- draw ----------
figure('Color','w', 'Name', 'iEEGLAB electrode values', 'NumberTitle','off');
hold on
surfFiles = o.surf_files;
if isempty(surfFiles) && isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'opt') && isfield(EEG.ieeglab.opt,'surf_files')
    surfFiles = EEG.ieeglab.opt.surf_files;
end
if ischar(surfFiles) || isstring(surfFiles), surfFiles = cellstr(surfFiles); end
for k = 1:numel(surfFiles)
    fk = char(surfFiles{k});
    if exist(fk,'file') ~= 2 && isfield(EEG,'filepath'), fk = fullfile(EEG.filepath, fk); end
    if exist(fk,'file') == 2 && exist('gifti','file') ~= 0
        try
            h = ieeg_RenderGifti(gifti(fk)); h.FaceAlpha = 0.08;
        catch ME
            warning('ieeglab_topoplot:surface', 'Could not render %s: %s', fk, ME.message);
        end
    end
end
nv = hasXYZ & ~isfinite(vals);
if any(nv)
    scatter3(xyz(nv,1), xyz(nv,2), xyz(nv,3), 12, [0.6 0.6 0.6], 'filled');
end
ok = hasXYZ & isfinite(vals);
scatter3(xyz(ok,1), xyz(ok,2), xyz(ok,3), 60, vals(ok), 'filled', 'MarkerEdgeColor', [0.2 0.2 0.2]);
colormap(parula); clim(info.clim);
cb = colorbar; ylabel(cb, desc);
axis equal off
if exist('ieeg_viewLight','file'), ieeg_viewLight(90, 0); else, view(3); end
title(info.title, 'Interpreter','tex', 'FontWeight','normal');
end

% ======================= helpers =======================

function [v, desc] = local_named_value(EEG, name, site, sigOnly)
C = EEG.nbchan;
v = nan(C,1);
labels = upper({EEG.chanlocs.labels});
M = [];
if isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'ccep_matrix'), M = EEG.ieeglab.ccep_matrix; end
if ~isempty(M)
    [ok, why] = ieeglab_matrix_current(EEG);
    if ~ok, error('ieeglab_topoplot:staleMatrix', 'The connectivity matrix is out of date: %s.', why); end
end
switch lower(name)
    case {'in_degree','indegree'}
        if isempty(M), error('ieeglab_topoplot:noMatrix', 'No CCEP matrix; run ieeglab_ccep_matrix first.'); end
        [tf, loc] = ismember(labels, upper(M.channels));
        v(tf) = M.in_degree(loc(tf));
        desc = 'In-degree (stimulation sites evoking a response)';
    case {'n1_amplitude','n1_latency','crp_tr','crp_explained_var','response'}
        if isempty(M), error('ieeglab_topoplot:noMatrix', 'No CCEP matrix; run ieeglab_ccep_matrix first.'); end
        if isempty(site)
            error('ieeglab_topoplot:noSite', '"%s" is per stimulation site; pass ''site'', e.g. ''%s''.', name, M.sites{1});
        end
        r = find(strcmpi(M.sites, site), 1);
        if isempty(r)
            cs = ieeglab_canonical_site({char(site)}, M.channels);
            r = find(strcmpi(M.sites, cs(1)), 1);
        end
        if isempty(r)
            error('ieeglab_topoplot:unknownSite', 'Site "%s" is not in the matrix. Sites: %s', site, strjoin(M.sites, ', '));
        end
        map = struct('n1_amplitude','amplitude_uv', 'n1_latency','latency_ms', ...
                     'crp_tr','tR_ms', 'crp_explained_var','explained_var', 'response','response');
        fld = map.(lower(name));
        if ~isfield(M, fld)
            error('ieeglab_topoplot:metricNotInMatrix', ...
                'The matrix was built from %s results and has no %s.', upper(M.source), fld);
        end
        [tf, loc] = ismember(labels, upper(M.channels));
        row = M.(fld)(r,:);
        if sigOnly && ~strcmp(fld, 'response')
            row(M.response(r,:) ~= 1) = NaN;    % tested, not significant: no N1 to show
        end
        v(tf) = row(loc(tf));
        descs = struct('amplitude_uv','N1 amplitude (\muV)', 'latency_ms','N1 latency (ms)', ...
            'tR_ms','CRP \tau_R (ms)', 'explained_var','CRP explained variance', 'response','Significant response');
        desc = descs.(fld);
        if sigOnly && ~strcmp(fld, 'response'), desc = [desc ', significant responses']; end
    otherwise
        error('ieeglab_topoplot:unknownMetric', ...
            'Unknown metric "%s". Use in_degree, n1_amplitude, n1_latency, crp_tr, crp_explained_var or response.', name);
end
end
