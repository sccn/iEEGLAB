function [EEG, com] = pop_ieeglab_reref(EEG, varargin)
% pop_ieeglab_reref() - iEEG re-referencing: choose the method and its options.
%
% Usage:
%   [EEG, com] = pop_ieeglab_reref(EEG)                       % dialog
%   EEG = pop_ieeglab_reref(EEG, 'method', 'carla')           % CCEP data
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
    'fraction', local_prev(prev, 'car_fraction', 0.25), 'rank', [], 'p_broad', 0.2);

if nargin < 2
    % ---- dialog
    if isCCEP
        labels  = {'CARLA, per stimulation site (Huang et al., 2024)', 'Common average (all good channels)', ...
                   'Lowest-covariance subset (Ojeda Valencia et al., 2023)', ...
                   'ICA, local components only (Michelmann et al., 2018)'};
        methods = {'carla', 'car', 'varsubset', 'ica'};
    else
        labels  = {'Common average (all good channels)', 'ICA, local components only (Michelmann et al., 2018)'};
        methods = {'car', 'ica'};
    end
    note = 'Clinician-bad channels are always left out of the reference.';
    if isCCEP, note = 'Clinician-bad channels and each epoch''s stimulated pair are always left out of the reference.'; end
    uilist = { ...
        {'style' 'text' 'string' 'Method'}, {'style' 'popupmenu' 'string' labels 'tag' 'method'}, ...
        {'style' 'text' 'string' 'Response window for ranking, ms (CARLA, subset)'}, ...
        {'style' 'edit' 'string' sprintf('%g %g', g.timewin) 'tag' 'timewin'}, ...
        {'style' 'text' 'string' 'Choose the reference per stimulation site (CARLA)'}, ...
        {'style' 'checkbox' 'string' '' 'value' double(g.persite) 'tag' 'persite'}, ...
        {'style' 'text' 'string' 'Sensitive cutoff (CARLA)'}, ...
        {'style' 'checkbox' 'string' '' 'value' double(g.sensitive) 'tag' 'sensitive'}, ...
        {'style' 'text' 'string' 'Share of channels kept (subset method, 0-1)'}, ...
        {'style' 'edit' 'string' num2str(g.fraction) 'tag' 'fraction'}, ...
        {'style' 'text' 'string' 'Number of components (ICA; empty = estimated rank)'}, ...
        {'style' 'edit' 'string' '' 'tag' 'rank'}, ...
        {'style' 'text' 'string' 'Broad component if chi-square p above (ICA)'}, ...
        {'style' 'edit' 'string' num2str(g.p_broad) 'tag' 'p_broad'}, ...
        {'style' 'text' 'string' note} };
    geom = {[1.3 1] [1.3 1] [1.3 1] [1.3 1] [1.3 1] [1.3 1] [1.3 1] 1};
    [res, ~, ~, out] = inputgui(geom, uilist, 'pophelp(''pop_ieeglab_reref'')', 'iEEGLAB - iEEG re-referencing');
    if isempty(res), return; end
    g.method = methods{out.method};
    tw = sscanf(out.timewin, '%f');
    if numel(tw) ~= 2 || tw(2) <= tw(1)
        error('pop_ieeglab_reref:badWindow', 'Enter the response window as "start end" in ms, e.g. 10 300.');
    end
    g.timewin = tw(:)';
    g.persite = logical(out.persite);
    g.sensitive = logical(out.sensitive);
    g.fraction = str2double(out.fraction);
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
opt.bad_channels = [];          % bad channels reach ieeglab_car through chanlocs.status
if strcmp(g.method, 'ica')
    EEG = ieeglab_icaref(EEG, struct('rank', g.rank, 'p_broad', g.p_broad, 'verbose', true));
else
    EEG = ieeglab_car(EEG, opt);
end
if ~isfield(EEG, 'ieeglab'), EEG.ieeglab = struct(); end
if ~isfield(EEG.ieeglab, 'opt') || ~isstruct(EEG.ieeglab.opt), EEG.ieeglab.opt = struct(); end
for f = {'car_method', 'car_timewin', 'car_sens', 'car_persite', 'car_fraction'}
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
if strcmp(g.method, 'ica')
    if ~isempty(g.rank), args = [args {'''rank''', num2str(g.rank)}]; end
    args = [args {'''p_broad''', num2str(g.p_broad)}];
end
com = sprintf('EEG = pop_ieeglab_reref(EEG, %s);', strjoin(args, ', '));
end

function v = local_prev(s, f, default)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = default; end
end

function out = iff(c, a, b)
if c, out = a; else, out = b; end
end
