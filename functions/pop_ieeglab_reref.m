function [EEG, com] = pop_ieeglab_reref(EEG, varargin)
% pop_ieeglab_reref() - iEEG re-referencing: choose the method and its options.
%
% Usage:
%   [EEG, com] = pop_ieeglab_reref(EEG)                       % dialog
%   EEG = pop_ieeglab_reref(EEG, 'method', 'carla')           % CCEP data
%   EEG = pop_ieeglab_reref(EEG, 'method', 'carla', 'sensitive', false, 'neighbors', 0)  % CARLA.m alone
%   EEG = pop_ieeglab_reref(EEG, 'method', 'car')
%   EEG = pop_ieeglab_reref(EEG, 'method', 'varsubset', 'fraction', 0.2)
%   EEG = pop_ieeglab_reref(EEG, 'method', 'ica')                % rank estimated
%   EEG = pop_ieeglab_reref(EEG, 'method', 'ica', 'rank', 12)
%   EEG = pop_ieeglab_reref(EEG, 'method', 'bipolar')            % neighbouring contacts
%   EEG = pop_ieeglab_reref(EEG, 'method', 'montage', 'montage', 'my_montage.tsv')
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
%                'bipolar'   - each contact minus the next one on its lead
%                              ('RA1-RA2' = RA1 - RA2), as FieldTrip's montage
%                              for depth electrodes. See ieeglab_montage.
%                'montage'   - a custom montage: a .tsv/.csv file or a table
%                              with columns positive, negative (one contact, or
%                              several for a local average, as a Laplacian),
%                              name; or a FieldTrip montage struct. See
%                              ieeglab_montage.
%   'montage'    the custom montage ('montage' method): file name, table or
%                struct.
%   'spacing_factor'  bipolar: with coordinates, contacts farther apart than
%                this times the lead's median spacing are not paired (the end
%                of one grid row and the start of the next). Default 2; 0 = off.
%   'timewin'    [start end] response window in ms used to rank channels
%                (carla, varsubset). Default [10 300] for carla, [15 500] for
%                varsubset (each method's publication scripts).
%   'sensitive'  CARLA's sensitive cutoff (true/false). Default true, as the
%                CARLA authors use it on real CCEP data.
%   'persite'    choose the reference separately for each stimulation site
%                (true/false). Default true.
%   'neighbors'  also leave out of the reference the contacts within this many
%                contacts of the stimulated pair, on the same lead (carla, car,
%                varsubset). Default 2 for carla, as its publication scripts,
%                0 for the others.
%   'fraction'   share of channels kept by 'varsubset' (0-1). Default 0.2 (HAPwave).
%   'rank'       number of ICA components ('ica'). Default [] = the effective rank
%                of the data, estimated; give it when known (e.g. after a common
%                average the rank is one less than the number of channels).
%   'p_broad'    chi-square p above which an ICA component counts as broad
%                ('ica'). Default 0.2.
%
% Continuous data are offered what needs no epoch: common average, ICA,
% bipolar and custom montages. CARLA and the lowest-covariance subset rank
% channels on each stimulation site's epochs, so they come after segmentation
% (iEEGLAB > Preprocess iEEG data); on stimulation (CCEP) data the common
% average and ICA also wait for epochs, where each pulse's stimulated pair is
% known and left out. Once a montage is applied the channels are derived ones
% ('RA1-RA2'), and no further re-referencing is offered: start again from the
% referential data to compare methods.
% Clinician-bad channels and, for CCEP data, the stimulated pair of every
% epoch are always left out of the reference. When preprocessing already
% removed a baseline, the same baseline correction is applied
% again after re-referencing, which gives exactly the result of re-referencing
% before it, the order ieeglab_preprocess uses: a per-epoch constant changes
% neither the reference channels (CARLA ranks by covariance and correlation)
% nor the re-referenced signal once its baseline is removed.
%
% The baseline is applied again after a montage too: a montage is linear, but
% the default baseline is a median, and the median of a difference is not the
% difference of the medians.
% The work is done by ieeglab_car, ieeglab_icaref, ieeglab_montage and
% ieeglab_apply_montage.
%
% Cedric Cannard, iEEGLAB, 2026

com = '';
ieeglab_require_data(EEG, 'pop_ieeglab_reref');
isEpoched = isfield(EEG, 'trials') && EEG.trials >= 2;
if isfield(EEG, 'ieeglab') && isfield(EEG.ieeglab, 'montage') && ~isempty(EEG.ieeglab.montage)
    error('pop_ieeglab_reref:afterMontage', ...
        ['These data are already re-referenced to a %s montage (channels such as %s). Re-reference ' ...
         'the dataset before the montage instead.'], EEG.ieeglab.montage.method, EEG.chanlocs(1).labels);
end
isCCEP = strcmp(ieeglab_detect_mode(EEG), 'ccep');
prev = struct();
if isfield(EEG, 'ieeglab') && isfield(EEG.ieeglab, 'opt') && isstruct(EEG.ieeglab.opt), prev = EEG.ieeglab.opt; end
% timewin, sensitive, fraction, neighbors: [] = the method's published default, or
% the value of the previous run when it used the same method (local_fill)
g = struct('method', iff(isCCEP, 'carla', 'car'), 'timewin', [], 'sensitive', [], ...
    'persite', local_prev(prev, 'car_persite', true), 'fraction', [], 'neighbors', [], ...
    'rank', [], 'p_broad', 0.2, 'montage', [], 'spacing_factor', 2);

if nargin < 2
    % ---- dialog: one section per method; options of the other methods are greyed out
    M = local_methods(isCCEP, isEpoched);
    if ~any(strcmp({M.name}, g.method)), g.method = M(1).name; end
    iSel = find(strcmp({M.name}, g.method), 1); if isempty(iSel), iSel = 1; end
    g = local_fill(g, prev, M(iSel).name);
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
    if any(strcmp({M.name}, 'carla'))
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
    if any(strcmp({M.name}, 'ica'))
        uilist = [uilist head('ICA') ...
            row('rank', '   Number of components (empty = estimated rank)', {'style' 'edit' 'string' ''}) ...
            row('p_broad', '   A component is broad (removed) if chi-square p above', {'style' 'edit' 'string' num2str(g.p_broad)})];
        geom = [geom {1 [1.3 1] [1.3 1]}];
        gv = [gv 1 1 1];
    end
    browse = {'style' 'pushbutton' 'string' '...' 'tag' 'browse' 'enable' onFor('montage') 'callback' @local_browse};
    uilist = [uilist head('Bipolar') ...
        row('spacing_factor', '   Skip pairs farther apart than N x the lead''s median spacing (0 = off)', {'style' 'edit' 'string' num2str(g.spacing_factor)}) ...
        head('Custom montage') ...
        row('montage', '   File (.tsv/.csv; columns positive, negative, name)', {'style' 'edit' 'string' ''}) {browse}];
    geom = [geom {1 [1.3 1] 1 [1.3 0.85 0.15]}];
    gv = [gv 1 1 1 1];
    note = 'Clinician-bad channels are always left out of the reference and are never paired.';
    if isCCEP && isEpoched, note = 'Clinician-bad channels and each epoch''s stimulated pair are always left out of the reference.'; end
    if isCCEP && ~isEpoched
        note = ['Stimulation data: common average, CARLA and ICA are offered after segmentation (Preprocess), ' ...
                'where each epoch''s stimulated pair is known and left out.'];
    end
    uilist = [uilist {{'style' 'text' 'string' note}}];
    geom = [geom {1}];
    gv = [gv 1];
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
    if isfield(out, 'rank'), g.rank = str2double(out.rank); if ~isfinite(g.rank), g.rank = []; end; end
    if isfield(out, 'p_broad'), g.p_broad = str2double(out.p_broad); end
    g.spacing_factor = str2double(out.spacing_factor);
    if strcmp(g.method, 'montage'), g.montage = strtrim(out.montage); end
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
g = local_fill(g, prev, g.method);                  % each method's published settings

if ~ismember(g.method, {'carla', 'car', 'varsubset', 'ica', 'bipolar', 'montage'})
    error('pop_ieeglab_reref:method', 'Unknown method ''%s'': use carla, car, varsubset, ica, bipolar or montage.', g.method);
end
isMontage = ismember(g.method, {'bipolar', 'montage'});
if ~isEpoched && ismember(g.method, {'carla', 'varsubset'})
    error('pop_ieeglab_reref:notEpoched', ...
        ['%s ranks the channels on each stimulation site''s response window: it needs epoched data. Run ' ...
         'iEEGLAB > Preprocess iEEG data with segmentation (epoching) first.'], upper(g.method));
end
if ~isEpoched && isCCEP && ismember(g.method, {'car', 'ica'})
    error('pop_ieeglab_reref:notEpoched', ...
        ['On continuous stimulation data the stimulated pair of each pulse cannot be left out, so its ' ...
         'artifact would enter every channel''s reference. Segment the data first (iEEGLAB > Preprocess ' ...
         'iEEG data), or use a bipolar or custom montage.']);
end
if strcmp(g.method, 'montage') && (isempty(g.montage) || ((ischar(g.montage) || isstring(g.montage)) && strtrim(string(g.montage)) == ""))
    error('pop_ieeglab_reref:noMontage', 'Give the custom montage: a .tsv/.csv file, a table or a FieldTrip montage struct.');
end
if strcmp(g.method, 'bipolar') && ~(isscalar(g.spacing_factor) && isfinite(g.spacing_factor) && g.spacing_factor >= 0)
    error('pop_ieeglab_reref:spacing', 'spacing_factor must be 0 (no check) or a positive number.');
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
hadBaseline = isEpoched && isfield(prev, 'apply_baseline') && isequal(prev.apply_baseline, true);
if ~isMontage && isfield(EEG, 'ieeglab') && isfield(EEG.ieeglab, 'car') && ~isempty(EEG.ieeglab.car)
    % (a bipolar difference cancels any reference common to its two contacts)
    warning('pop_ieeglab_reref:alreadyReferenced', ...
        ['These data were already re-referenced (%s). Re-referencing again stacks the two; to compare ' ...
         'methods, start each from the same preprocessed dataset.'], char(string(EEG.ref)));
end

if isMontage
    % ---- run: montage
    if strcmp(g.method, 'bipolar')
        Mt = ieeglab_montage(EEG, 'bipolar', 'spacing_factor', g.spacing_factor);
    else
        Mt = ieeglab_montage(EEG, g.montage);
    end
    EEG = ieeglab_apply_montage(EEG, Mt);
    if ~isfield(EEG.ieeglab, 'opt') || ~isstruct(EEG.ieeglab.opt), EEG.ieeglab.opt = struct(); end
    EEG.ieeglab.opt.montage_method = Mt.method;
    EEG.ieeglab.opt.montage_source = Mt.source;
    if hadBaseline
        % a median (default) or trimmed-mean baseline of the new channels is not
        % the montage of the old ones' baselines
        EEG = ieeglab_rm_baseline(EEG);
    end
    EEG = eeg_checkset(EEG);
    if strcmp(g.method, 'bipolar')
        com = sprintf('EEG = pop_ieeglab_reref(EEG, ''method'', ''bipolar'', ''spacing_factor'', %s);', mat2str(g.spacing_factor));
    elseif ischar(g.montage) || isstring(g.montage)
        com = sprintf('EEG = pop_ieeglab_reref(EEG, ''method'', ''montage'', ''montage'', ''%s'');', strrep(char(g.montage), '''', ''''''));
    else
        com = 'EEG = pop_ieeglab_reref(EEG, ''method'', ''montage'', ''montage'', montage); % montage: the table or struct given';
    end
    return
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
if strcmp(g.method, 'varsubset'), args = [args {'''fraction''', mat2str(g.fraction)}]; end
if ~strcmp(g.method, 'ica')                         % always: the default depends on the method
    args = [args {'''neighbors''', num2str(g.neighbors)}];
end
if strcmp(g.method, 'ica')
    if ~isempty(g.rank), args = [args {'''rank''', num2str(g.rank)}]; end
    args = [args {'''p_broad''', mat2str(g.p_broad)}];
end
com = sprintf('EEG = pop_ieeglab_reref(EEG, %s);', strjoin(args, ', '));
end

function M = local_methods(isCCEP, isEpoched)
% Methods offered in the dialog, a one-line description, and the options each uses
M = struct( ...
    'name',  {'carla', 'car', 'varsubset', 'ica', 'bipolar', 'montage'}, ...
    'label', {'CARLA, per stimulation site (Huang et al., 2024)', 'Common average (all good channels)', ...
              'Lowest-covariance subset (Ojeda Valencia et al., 2023)', ...
              'ICA, local components only (Michelmann et al., 2018)', ...
              'Bipolar, neighbouring contacts on each lead', 'Custom montage (from a file)'}, ...
    'desc',  {'For each stimulation site, averages only the channels that carry no evoked response.', ...
              'Averages all good channels.', ...
              'Averages the given share of channels with the lowest covariance in the response window.', ...
              'Removes the ICA components spread evenly over all contacts (the reference and other shared signals).', ...
              'Each contact minus the next one on its lead: RA1-RA2 = RA1 - RA2.', ...
              'One new channel per row: positive contact minus negative contact(s), e.g. a bipolar pair or a local average.'}, ...
    'uses',  {{'timewin','sensitive','neighbors'}, {'neighbors'}, {'timewin','fraction','neighbors'}, {'rank','p_broad'}, ...
              {'spacing_factor'}, {'montage','browse'}});
% Continuous data: what needs no epoch. Stimulation data: only the fixed
% montages before segmentation, since the stimulated pair is known per epoch.
if ~isCCEP, M = M(~ismember({M.name}, {'carla','varsubset'})); M(1).uses = {}; end
if ~isEpoched && isCCEP, M = M(ismember({M.name}, {'bipolar','montage'})); end
end

function local_browse(h, ~)
% File picker for the custom montage
[f, p] = uigetfile({'*.tsv;*.csv;*.txt', 'Montage table (*.tsv, *.csv, *.txt)'}, 'Custom montage');
if isequal(f, 0), return; end
set(findobj(ancestor(h, 'figure'), 'tag', 'montage'), 'string', fullfile(p, f));
end

function local_update(h, M)
% Grey out the options the selected method does not use, and describe it
fig = ancestor(h, 'figure');
m = M(get(h, 'value'));
for tag = {'timewin','sensitive','fraction','neighbors','rank','p_broad','spacing_factor','montage','browse'}
    st = iff(ismember(tag{1}, m.uses), 'on', 'off');
    set(findobj(fig, 'tag', tag{1}), 'enable', st);
    set(findobj(fig, 'tag', ['lbl_' tag{1}]), 'enable', st);
end
set(findobj(fig, 'tag', 'desc'), 'string', m.desc);
[tw, nb] = ieeglab_car_defaults(m.name);   % switching method shows that method's defaults
set(findobj(fig, 'tag', 'timewin'), 'string', sprintf('%g %g', tw));
set(findobj(fig, 'tag', 'neighbors'), 'string', num2str(nb));
end

function g = local_fill(g, prev, method)
% Options not given: the previous run's value when it used the same method,
% else the method's published default (ieeglab_car_defaults).
[tw0, nb0] = ieeglab_car_defaults(method);
same = isfield(prev, 'car_method') && strcmpi(char(string(prev.car_method)), method);
pick = @(f, def) iff(same, local_prev(prev, f, def), def);
if isempty(g.timewin),   g.timewin   = pick('car_timewin', tw0); end
if isempty(g.neighbors), g.neighbors = pick('car_neighbors', nb0); end
if isempty(g.sensitive), g.sensitive = pick('car_sens', true); end
if isempty(g.fraction),  g.fraction  = pick('car_fraction', 0.2); end
end

function v = local_prev(s, f, default)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = default; end
end

function out = iff(c, a, b)
if c, out = a; else, out = b; end
end
