function [EEG, com] = pop_ieeglab_reref(EEG, varargin)
% pop_ieeglab_reref() - iEEG re-referencing: choose the method and its options.
%
% Usage:
%   [EEG, com] = pop_ieeglab_reref(EEG)                       % dialog
%   EEG = pop_ieeglab_reref(EEG, 'method', 'carla')           % CCEP data
%   EEG = pop_ieeglab_reref(EEG, 'method', 'carla', 'sensitive', true, 'neighbors', 2)
%   EEG = pop_ieeglab_reref(EEG, 'method', 'car')
%   EEG = pop_ieeglab_reref(EEG, 'method', 'varsubset', 'fraction', 0.25)
%   EEG = pop_ieeglab_reref(EEG, 'method', 'ica')                % rank estimated
%   EEG = pop_ieeglab_reref(EEG, 'method', 'ica', 'rank', 12)
%
% Options (name/value):
%   'method'     'carla'     - CAR by Least Anticorrelation, per stimulation site
%                              (Huang et al., 2024). CCEP data only. Default for CCEP.
%                'car'       - common average of all good channels. Default otherwise.
%                'varsubset' - lowest-covariance fixed fraction of channels
%                              (Ojeda Valencia et al., 2023). Kept for comparison.
%                'ica'       - ICA re-referencing (Michelmann et al., 2018): the
%                              components spread uniformly over the contacts (the
%                              reference and other shared signals) are removed.
%                              See ieeglab_icaref.
%   'timewin'    [start end] response window in ms used to rank channels
%                (carla, varsubset). Default [10 300].
%   'sensitive'  CARLA's sensitive cutoff (true/false). Default false.
%   'persite'    choose the reference separately for each stimulation site
%                (true/false). Default true.
%   'neighbors'  also leave out of the reference the contacts within this many
%                contacts of the stimulated pair, on the same lead (carla, car,
%                varsubset). Default 0. The CARLA publication scripts use 2.
%   'fraction'   share of channels kept by 'varsubset' (0-1). Default 0.25.
%   'rank'       number of ICA components ('ica'). Default [] = the effective rank
%                of the data, estimated; give it when known (e.g. after a common
%                average the rank is one less than the number of channels).
%   'p_broad'    chi-square p above which an ICA component counts as broad
%                ('ica'). Default 0.2.
%
% Runs on epoched data (iEEGLAB > Preprocess iEEG data, with segmentation).
% Clinician-bad channels and, for CCEP data, the stimulated pair of every
% epoch are always left out of the reference. When preprocessing already
% removed a baseline, the same baseline correction is applied
% again after re-referencing, which gives exactly the result of re-referencing
% before it, the order ieeglab_preprocess uses: a per-epoch constant changes
% neither the reference channels (CARLA ranks by covariance and correlation)
% nor the re-referenced signal once its baseline is removed.
%
% This is where further iEEG methods (bipolar, Laplacian) will be offered.
% The work is done by ieeglab_car and ieeglab_icaref.
%
% Cedric Cannard, iEEGLAB, 2026

com = '';
ieeglab_require_data(EEG, 'pop_ieeglab_reref');
if ~isfield(EEG, 'trials') || EEG.trials < 2
    error('pop_ieeglab_reref:notEpoched', ...
        ['Re-referencing needs epoched data. Run iEEGLAB > Preprocess iEEG data with ' ...
         'segmentation (epoching) first.']);
end
isCCEP = strcmp(ieeglab_detect_mode(EEG), 'ccep');
prev = struct();
if isfield(EEG, 'ieeglab') && isfield(EEG.ieeglab, 'opt') && isstruct(EEG.ieeglab.opt), prev = EEG.ieeglab.opt; end
g = struct('method', iff(isCCEP, 'carla', 'car'), 'timewin', local_prev(prev, 'car_timewin', [10 300]), ...
    'sensitive', local_prev(prev, 'car_sens', false), 'persite', local_prev(prev, 'car_persite', true), ...
    'fraction', local_prev(prev, 'car_fraction', 0.25), 'neighbors', local_prev(prev, 'car_neighbors', 0), ...
    'rank', [], 'p_broad', 0.2);

if nargin < 2
    % ---- dialog: one section per method; options of the other methods are greyed out
    M = local_methods(isCCEP);
    iSel = find(strcmp({M.name}, g.method), 1); if isempty(iSel), iSel = 1; end
    onFor = @(tag) iff(ismember(tag, M(iSel).uses), 'on', 'off');
    row = @(tag, label, ctrl) {{'style' 'text' 'string' label 'tag' ['lbl_' tag] 'enable' onFor(tag)}, ...
                               [ctrl {'tag' tag 'enable' onFor(tag)}]};
    head = @(s) {{'style' 'text' 'string' s 'fontweight' 'bold'}};
    uilist = [{{'style' 'text' 'string' 'Method'}, ...
               {'style' 'popupmenu' 'string' {M.label} 'value' iSel 'tag' 'method' ...
                'callback' @(h,~) local_update(h, M)}}, ...
              {{'style' 'text' 'string' M(iSel).desc 'tag' 'desc'}}];
    geom = {[1 1.3] 1};
    gv = [1 1.6];
    if isCCEP
        uilist = [uilist head('CARLA and lowest-covariance subset') ...
            row('timewin', '   Response window, ms (channels are ranked on it)', {'style' 'edit' 'string' sprintf('%g %g', g.timewin)}) ...
            head('CARLA') ...
            row('sensitive', '   Sensitive cutoff (smaller reference, stops at the first significant drop)', {'style' 'checkbox' 'string' '' 'value' double(g.sensitive)}) ...
            head('Lowest-covariance subset') ...
            row('fraction', '   Share of channels kept (0-1)', {'style' 'edit' 'string' num2str(g.fraction)}) ...
            head('Stimulated pair (CARLA, common average, subset)') ...
            row('neighbors', '   Also leave out N contacts on each side, same lead', {'style' 'edit' 'string' num2str(g.neighbors)})];
        geom = [geom {1 [1.3 1] 1 [1.3 1] 1 [1.3 1] 1 [1.3 1]}];
        gv = [gv 1 1 1 1 1 1 1 1];
    end
    uilist = [uilist head('ICA') ...
        row('rank', '   Number of components (empty = estimated rank)', {'style' 'edit' 'string' ''}) ...
        row('p_broad', '   A component is broad (removed) if chi-square p above', {'style' 'edit' 'string' num2str(g.p_broad)})];
    note = 'Clinician-bad channels are always left out of the reference.';
    if isCCEP, note = 'Clinician-bad channels and each epoch''s stimulated pair are always left out of the reference.'; end
    uilist = [uilist {{'style' 'text' 'string' note}}];
    geom = [geom {1 [1.3 1] [1.3 1] 1}];
    gv = [gv 1 1 1 1];
    [res, ~, ~, out] = inputgui('geometry', geom, 'geomvert', gv, 'uilist', uilist, ...
        'helpcom', 'pophelp(''pop_ieeglab_reref'')', 'title', 'iEEGLAB - iEEG re-referencing', 'minwidth', 640);
    if isempty(res), return; end
    g.method = M(out.method).name;
    if isfield(out, 'timewin')
        tw = sscanf(out.timewin, '%f');
        if ismember('timewin', M(out.method).uses) && (numel(tw) ~= 2 || tw(2) <= tw(1))
            error('pop_ieeglab_reref:badWindow', 'Enter the response window as "start end" in ms, e.g. 10 300.');
        end
        if numel(tw) == 2, g.timewin = tw(:)'; end
    end
    if isfield(out, 'sensitive'), g.sensitive = logical(out.sensitive); end
    if isfield(out, 'fraction'),  g.fraction = str2double(out.fraction); end
    if isfield(out, 'neighbors'), g.neighbors = str2double(out.neighbors); end
    g.rank = str2double(out.rank); if ~isfinite(g.rank), g.rank = []; end
    g.p_broad = str2double(out.p_broad);
else
    % ---- name/value pairs
    if mod(numel(varargin), 2), error('pop_ieeglab_reref:args', 'Options come in name/value pairs.'); end
    for k = 1:2:numel(varargin)
        name = lower(char(varargin{k}));
        if ~isfield(g, name), error('pop_ieeglab_reref:unknownOption', 'Unknown option ''%s''.', varargin{k}); end
        g.(name) = varargin{k+1};
    end
    g.method = lower(char(g.method));
end

if ~ismember(g.method, {'carla', 'car', 'varsubset', 'ica'})
    error('pop_ieeglab_reref:method', 'Unknown method ''%s'': use carla, car, varsubset or ica.', g.method);
end
if strcmp(g.method, 'ica') && ~(isscalar(g.p_broad) && g.p_broad > 0 && g.p_broad < 1)
    error('pop_ieeglab_reref:pBroad', 'p_broad must be between 0 and 1.');
end
if strcmp(g.method, 'carla') && ~isCCEP
    error('pop_ieeglab_reref:carlaNeedsCCEP', ...
        ['CARLA needs CCEP data: it ranks channels against each epoch''s stimulated pair, and the ' ...
         'event types of this dataset do not name contact pairs. Use ''car''.']);
end
if ~(isscalar(g.neighbors) && isfinite(g.neighbors) && g.neighbors >= 0 && g.neighbors == round(g.neighbors))
    error('pop_ieeglab_reref:neighbors', 'neighbors must be a whole number of contacts (0 = only the stimulated pair).');
end
if strcmp(g.method, 'varsubset') && ~(isscalar(g.fraction) && g.fraction > 0 && g.fraction <= 1)
    error('pop_ieeglab_reref:fraction', 'The share of channels kept must be between 0 and 1.');
end
hadBaseline = isfield(prev, 'apply_baseline') && isequal(prev.apply_baseline, true);
if isfield(EEG, 'ieeglab') && isfield(EEG.ieeglab, 'car') && ~isempty(EEG.ieeglab.car)
    warning('pop_ieeglab_reref:alreadyReferenced', ...
        ['These data were already re-referenced (%s). Re-referencing again stacks the two; to compare ' ...
         'methods, start each from the same preprocessed dataset.'], char(string(EEG.ref)));
end

% ---- run
opt = prev;
opt.car_method  = g.method;
opt.car_timewin = g.timewin;
opt.car_sens    = logical(g.sensitive);
opt.car_persite = logical(g.persite);
opt.car_fraction = g.fraction;
opt.car_neighbors = g.neighbors;
opt.bad_channels = [];          % bad channels reach ieeglab_car through chanlocs.status
if strcmp(g.method, 'ica')
    EEG = ieeglab_icaref(EEG, struct('rank', g.rank, 'p_broad', g.p_broad, 'verbose', true));
else
    EEG = ieeglab_car(EEG, opt);
end
if ~isfield(EEG, 'ieeglab'), EEG.ieeglab = struct(); end
if ~isfield(EEG.ieeglab, 'opt') || ~isstruct(EEG.ieeglab.opt), EEG.ieeglab.opt = struct(); end
for f = {'car_method', 'car_timewin', 'car_sens', 'car_persite', 'car_fraction', 'car_neighbors'}
    EEG.ieeglab.opt.(f{1}) = opt.(f{1});
end
if hadBaseline
    % same window, method and mode as in preprocessing (read from EEG.ieeglab.opt)
    EEG = ieeglab_rm_baseline(EEG);
end
EEG = eeg_checkset(EEG);

args = {'''method''', ['''' g.method '''']};
if ismember(g.method, {'carla', 'varsubset'}), args = [args {'''timewin''', mat2str(g.timewin)}]; end
if strcmp(g.method, 'carla'), args = [args {'''sensitive''', mat2str(logical(g.sensitive)), '''persite''', mat2str(logical(g.persite))}]; end
if strcmp(g.method, 'varsubset'), args = [args {'''fraction''', num2str(g.fraction)}]; end
if ~strcmp(g.method, 'ica') && g.neighbors > 0, args = [args {'''neighbors''', num2str(g.neighbors)}]; end
if strcmp(g.method, 'ica')
    if ~isempty(g.rank), args = [args {'''rank''', num2str(g.rank)}]; end
    args = [args {'''p_broad''', num2str(g.p_broad)}];
end
com = sprintf('EEG = pop_ieeglab_reref(EEG, %s);', strjoin(args, ', '));
end

function M = local_methods(isCCEP)
% Methods offered in the dialog, a one-line description, and the options each uses
M = struct( ...
    'name',  {'carla', 'car', 'varsubset', 'ica'}, ...
    'label', {'CARLA, per stimulation site (Huang et al., 2024)', 'Common average (all good channels)', ...
              'Lowest-covariance subset (Ojeda Valencia et al., 2023)', ...
              'ICA, local components only (Michelmann et al., 2018)'}, ...
    'desc',  {'For each stimulation site, averages only the channels that carry no evoked response.', ...
              'Averages all good channels.', ...
              'Averages the given share of channels with the lowest covariance in the response window.', ...
              'Removes the ICA components spread evenly over all contacts (the reference and other shared signals).'}, ...
    'uses',  {{'timewin','sensitive','neighbors'}, {'neighbors'}, {'timewin','fraction','neighbors'}, {'rank','p_broad'}});
if ~isCCEP, M = M(ismember({M.name}, {'car','ica'})); M(1).uses = {}; end
end

function local_update(h, M)
% Grey out the options the selected method does not use, and describe it
fig = ancestor(h, 'figure');
m = M(get(h, 'value'));
for tag = {'timewin','sensitive','fraction','neighbors','rank','p_broad'}
    st = iff(ismember(tag{1}, m.uses), 'on', 'off');
    set(findobj(fig, 'tag', tag{1}), 'enable', st);
    set(findobj(fig, 'tag', ['lbl_' tag{1}]), 'enable', st);
end
set(findobj(fig, 'tag', 'desc'), 'string', m.desc);
end

function v = local_prev(s, f, default)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = default; end
end

function out = iff(c, a, b)
if c, out = a; else, out = b; end
end
