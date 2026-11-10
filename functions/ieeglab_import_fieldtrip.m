function [EEG, com] = ieeglab_import_fieldtrip(data, varargin)
% ieeglab_import_fieldtrip() - FieldTrip iEEG data (epochs) and electrode
%                              positions into an EEGLAB dataset for iEEGLAB.
%
% Usage:
%   [EEG, com] = ieeglab_import_fieldtrip                % dialogs (iEEGLAB > Import FieldTrip epochs)
%   EEG = ieeglab_import_fieldtrip('SubjectUCI29_data.mat', 'elec', 'SubjectUCI29_elec_acpc_f.mat')
%   EEG = ieeglab_import_fieldtrip(data, 'elec', elec_acpc_f, ...
%             'types', {'LPG*', 'ECOG'; 'LTG*', 'ECOG'; '*', 'SEEG'})
%
% data   FieldTrip raw data (output of ft_preprocessing: label, trial, time,
%        fsample, optional trialinfo, elec and chanunit), or a .mat file holding it.
%        Every trial must have the same time axis (epochs), and it must contain
%        time 0. A single trial is continuous data.
%
% Options (name/value):
%   'elec'       FieldTrip electrode structure (label, chanpos or elecpos, unit,
%                coordsys), a structure holding one in field elec (e.g. FieldTrip
%                data or freq), or a .mat file holding either. Default: data.elec
%                when present. Positions are matched to the channels by label
%                (case ignored): chanpos when it has one row per label (FieldTrip
%                pairs them), else elecpos.
%   'types'      channel types: n x 2 cell of {label pattern, type}, first match
%                wins (* and ? wildcards, case ignored), e.g. {'LPG*','ECOG'; '*','SEEG'}.
%                Default: none (EEGLAB's empty type).
%   'event_type' type (text or number) of the event placed at time 0 of every
%                epoch. Default 'stimulus'.
%   'verbose'    default true
%
% Events: FieldTrip's data have no event list. Each epoch gets one event at time
% 0, with the epoch's row of trialinfo in event field 'trialinfo' (FieldTrip's
% tutorials keep reaction times or conditions there). When time 0 falls between
% two samples, the event's latency is fractional. Continuous data (one trial) get
% no event, and their time starts at 0 (EEGLAB's convention).
%
% Positions are converted to mm from elec.unit: m, dm, cm, mm or um (also µm);
% another unit is an error. Without a unit, the unit is estimated from the size
% of the layout as FieldTrip's ft_estimate_units does (the diagonal of the box
% around the positions, against a head of about 0.1 to 0.2 m: 0.05 to 0.5 is m,
% 0.5 to 5 dm, 5 to 50 cm, 50 to 500 mm), and a warning names the unit assumed.
% The estimate can be wrong for a small layout (one depth lead in mm reads as
% cm): set elec.unit.
%
% Values: EEGLAB holds iEEG in uV. The unit FieldTrip gives the channels
% (data.chanunit, else data.elec.chanunit, data.hdr.chanunit, then elec.chanunit)
% is checked against the size of the values: the median absolute deviation of
% each channel around its median (up to 20000 samples), median over the
% channels, which for iEEG in uV lies between 0.5 and 500. V or mV of that size
% is converted to uV, and said so. A unit that gives a size outside that range
% is taken as a wrong label: the values are kept and a warning says so (FieldTrip
% often labels electrodes 'V' by default: SubjectUCI29 says 'V', its values are uV).
%
% Coordinates go to EEG.chanlocs X, Y, Z; their space (FieldTrip's coordsys, e.g.
% 'acpc') to EEG.ieeglab.coordsys. Brain surfaces of the same space are chosen
% with iEEGLAB > Choose brain surfaces (pop_ieeglab_surfaces). Data read from a
% file keep its folder and name in EEG.filepath and EEG.filename.
%
% History (com): files are written by name, structures by the name of the
% caller's variable, so re-running the history needs that variable in the
% workspace. When a structure is not a variable (e.g. s.elec), com is empty.
%
% Cedric Cannard, iEEGLAB, 2026

nameData = ''; nameElec = '';                   % the caller's variables, for the history
if nargin < 1 || isempty(data)
    EEG = []; com = '';
    [data, varargin] = local_dialog();
    if isempty(data), return; end
else
    if isstruct(data), try nameData = inputname(1); catch, end; end
    for j = 1:2:numel(varargin) - 1
        a = varargin{j};
        if (ischar(a) || isstring(a)) && strncmpi('elec', a, max(2, strlength(a))) && isstruct(varargin{j + 1})
            try nameElec = inputname(j + 2); catch, nameElec = ''; end
        end
    end
end
p = inputParser;
p.addParameter('elec', []);
p.addParameter('types', {});
p.addParameter('event_type', 'stimulus');
p.addParameter('verbose', true);
p.parse(varargin{:});
o = p.Results;
if isstring(o.event_type), o.event_type = char(o.event_type); end

src = '';
if ischar(data) || isstring(data)
    src = char(data);
    data = local_load_struct(src, {'trial', 'label'}, 'FieldTrip data (label, trial, time)');
end
if ~all(isfield(data, {'label', 'trial', 'time'}))
    error('ieeglab_import_fieldtrip:notRaw', 'data needs the fields label, trial and time (FieldTrip raw data).');
end
elec = o.elec;
if ischar(elec) || isstring(elec)
    elec = local_load_elec(char(elec));
elseif ~isempty(elec)
    elec = local_elec(elec, 'elec');
elseif isfield(data, 'elec')
    if isempty(local_pos(data.elec))
        warning('ieeglab_import_fieldtrip:noPositions', ...
            'data.elec has no positions (chanpos or elecpos, one row per label): imported without coordinates.');
    else
        elec = data.elec;
    end
end

% ---------------- epochs ----------------
nTr = numel(data.trial);
T = numel(data.time{1});
if any(cellfun(@numel, data.time) ~= T) || any(cellfun(@(x) max(abs(x - data.time{1})), data.time) > 1e-9)
    error('ieeglab_import_fieldtrip:unequalTrials', ...
        'The trials have different time axes: cut them to one window first (ft_redefinetrial, toilim).');
end
C = numel(data.label);
t = double(data.time{1});
if isfield(data, 'fsample'), fs = double(data.fsample); else, fs = 1 / mean(diff(t)); end
lat0 = [];
if nTr > 1                                      % epochs: the event goes at time 0
    if t(1) > 1e-3 / fs || t(end) < -1e-3 / fs
        error('ieeglab_import_fieldtrip:noTimeZero', ['The epochs (%g to %g ms) do not contain time 0, where ' ...
            'the event of each epoch goes: shift their time axis or cut them around the event (ft_redefinetrial).'], ...
            1000 * t(1), 1000 * t(end));
    end
    lat0 = 1 - t(1) * fs;                       % sample of time 0, fractional when between samples
    if abs(lat0 - round(lat0)) < 1e-6, lat0 = round(lat0); end
end

% ---------------- units of the values ----------------
[cu, cuFrom] = local_chanunit(data, elec);
gain = ones(C, 1);                              % to uV
unitMsg = {}; unitWarn = {};
if ~isempty(cuFrom)
    [gain, unitMsg, unitWarn] = local_to_uv(cu, local_size(data, T, nTr), cuFrom);
end
X = zeros(C, T, nTr, 'single');
if all(gain == 1)
    for k = 1:nTr, X(:, :, k) = single(data.trial{k}); end
else
    for k = 1:nTr, X(:, :, k) = single(gain .* double(data.trial{k})); end
end

EEG = eeg_emptyset;
EEG.setname = 'FieldTrip import';
if ~isempty(src)
    [~, EEG.setname] = fileparts(src);
    [EEG.filepath, EEG.filename] = local_file(src);
end
EEG.data = X;
EEG.nbchan = C;
EEG.pnts = T;
EEG.trials = nTr;
EEG.srate = fs;
if nTr > 1
    EEG.xmin = t(1);
    EEG.xmax = t(end);
    EEG.times = t * 1000;
else                                            % continuous: time starts at 0
    EEG.xmin = 0;
    EEG.xmax = (T - 1) / fs;
    EEG.times = (0:T - 1) / fs * 1000;
end
EEG.ref = 'unknown';
EEG.chanlocs = struct('labels', cellstr(data.label(:))', 'type', '', 'X', [], 'Y', [], 'Z', []);

% ---------------- channel types ----------------
for r = 1:size(o.types, 1)
    pat = ['^' regexprep(regexptranslate('escape', char(o.types{r, 1})), {'\\\*', '\\\?'}, {'.*', '.'}) '$'];
    for c = find(cellfun(@isempty, {EEG.chanlocs.type}))
        if ~isempty(regexp(data.label{c}, pat, 'once', 'ignorecase'))
            EEG.chanlocs(c).type = upper(char(o.types{r, 2}));
        end
    end
end

% ---------------- coordinates, by label ----------------
coordsys = '';
hasPos = false(C, 1);
unitNote = '';
if ~isempty(elec)
    pos = local_pos(elec);
    [scale, unitNote] = local_mm(elec, pos);
    if isfield(elec, 'coordsys'), coordsys = char(elec.coordsys); end
    [found, loc] = ismember(lower(data.label), lower(elec.label));
    for c = find(found(:))'
        xyz = pos(loc(c), :) * scale;
        if all(isfinite(xyz))
            EEG.chanlocs(c).X = xyz(1); EEG.chanlocs(c).Y = xyz(2); EEG.chanlocs(c).Z = xyz(3);
            hasPos(c) = true;
        end
    end
end

% ---------------- one event at time 0 of each epoch ----------------
if nTr > 1
    ti = [];
    if isfield(data, 'trialinfo'), ti = data.trialinfo; end
    ev = struct('type', {}, 'latency', {}, 'epoch', {}, 'trialinfo', {}, 'urevent', {});
    for k = 1:nTr
        info = [];
        if ~isempty(ti), info = ti(k, :); end
        ev(k) = struct('type', o.event_type, 'latency', (k - 1) * T + lat0, 'epoch', k, 'trialinfo', info, 'urevent', k);
    end
    EEG.event = ev;
    EEG.urevent = rmfield(ev, {'epoch', 'urevent'});
end
EEG.ieeglab = struct('coordsys', coordsys, 'source', 'fieldtrip');
EEG = eeg_checkset(EEG, 'eventconsistency');

if o.verbose
    if nTr > 1
        fprintf('[FieldTrip import] %d channels, %d epochs of %.3f s (%g to %g ms) at %g Hz; %d channels with coordinates%s.\n', ...
            C, nTr, T / fs, EEG.times(1), EEG.times(end), fs, nnz(hasPos), local_space(coordsys));
    else
        fprintf(['[FieldTrip import] %d channels, continuous (one trial, no event added), %.3f s at %g Hz; ' ...
            '%d channels with coordinates%s.\n'], C, T / fs, fs, nnz(hasPos), local_space(coordsys));
    end
    for k = 1:numel(unitMsg), fprintf('[FieldTrip import] %s\n', unitMsg{k}); end
    fprintf('[FieldTrip import] Next: iEEGLAB > Choose brain surfaces, re-referencing, Frequency analysis.\n');
end
for k = 1:numel(unitWarn), warning('ieeglab_import_fieldtrip:chanunit', '%s', unitWarn{k}); end
if ~isempty(unitNote), warning('ieeglab_import_fieldtrip:unitAssumed', '%s', unitNote); end
if ~isempty(elec) && ~all(hasPos)
    miss = data.label(~hasPos);
    more = ''; if numel(miss) > 20, more = ' ...'; end
    warning('ieeglab_import_fieldtrip:noPosition', '%d channel(s) have no position in elec: %s%s', ...
        numel(miss), strjoin(reshape(miss(1:min(end, 20)), 1, []), ' '), more);
end
com = local_history(src, nameData, nameElec, o);
end

% ===================== local helpers =====================

function com = local_history(src, nameData, nameElec, o)
% A history line that re-runs the import: files by name, structures by the
% caller's variable name; empty when a structure has no name
com = '';
if ~isempty(src), args = {local_q(src)}; elseif ~isempty(nameData), args = {nameData}; else, return; end
if ischar(o.elec) || isstring(o.elec)
    args{end+1} = ['''elec'', ' local_q(o.elec)];
elseif ~isempty(o.elec)
    if isempty(nameElec), return; end
    args{end+1} = ['''elec'', ' nameElec];
end
if ~isempty(o.types)
    rows = cellfun(@(a, b) [local_q(a) ', ' local_q(b)], o.types(:, 1), o.types(:, 2), 'UniformOutput', false);
    args{end+1} = ['''types'', {' strjoin(rows', '; ') '}'];
end
if ~isequal(o.event_type, 'stimulus')
    if isnumeric(o.event_type) || islogical(o.event_type), v = mat2str(o.event_type); else, v = local_q(o.event_type); end
    args{end+1} = ['''event_type'', ' v];
end
com = sprintf('EEG = ieeglab_import_fieldtrip(%s);', strjoin(args, ', '));
end

function s = local_q(x)
s = ['''' strrep(char(x), '''', '''''') ''''];
end

function [file, args] = local_dialog()
% The data file, then the electrode file, channel types and event type
file = ''; args = {};
[f, p] = uigetfile({'*.mat', 'FieldTrip data (*.mat)'}, 'SELECT THE FIELDTRIP DATA (.mat, output of ft_preprocessing)');
if isequal(f, 0), return; end
browse = ['[f_, p_] = uigetfile(''*.mat'', ''SELECT THE ELECTRODE POSITIONS (.mat with a FieldTrip elec)''); ' ...
          'if ~isequal(f_, 0), set(findobj(gcbf, ''tag'', ''elec''), ''string'', fullfile(p_, f_)); end; clear f_ p_;'];
uilist = { ...
    {'style' 'text' 'string' f 'fontweight' 'bold'}, ...
    {'style' 'text' 'string' 'Electrode positions: .mat holding a FieldTrip elec (empty: the data''s own elec, if any)'}, ...
    {'style' 'edit' 'string' '' 'tag' 'elec'}, {'style' 'pushbutton' 'string' 'Browse' 'callback' browse}, ...
    {'style' 'text' 'string' 'Channel types: label pattern and type, separated by ; (e.g. LPG* ECOG; LTG* ECOG; * SEEG)'}, ...
    {'style' 'edit' 'string' '' 'tag' 'types'}, ...
    {'style' 'text' 'string' 'Event at time 0 of each epoch, type'}, {'style' 'edit' 'string' 'stimulus' 'tag' 'evtype'} };
uigeom = {1, 1, [1.6 0.4], 1, 1, [1.2 0.8]};
[res, ~, ~, out] = inputgui('geometry', uigeom, 'uilist', uilist, 'minwidth', 600, ...
    'title', 'iEEGLAB - Import FieldTrip epochs', 'helpcom', 'pophelp(''ieeglab_import_fieldtrip'');');
if isempty(res), return; end
file = fullfile(p, f);
if ~isempty(strtrim(out.elec)), args = [args {'elec', strtrim(out.elec)}]; end
types = cell(0, 2);
rows = strtrim(strsplit(out.types, ';'));
for r = rows(~cellfun(@isempty, rows))
    tok = strsplit(r{1});
    if numel(tok) ~= 2
        error('ieeglab_import_fieldtrip:types', 'Channel types: "%s" is not a label pattern and a type, e.g. LPG* ECOG.', r{1});
    end
    types(end + 1, :) = tok; %#ok<AGROW>
end
if ~isempty(types), args = [args {'types', types}]; end
ev = strtrim(out.evtype);
if ~isempty(ev) && ~strcmp(ev, 'stimulus'), args = [args {'event_type', ev}]; end
end

function s = local_load_struct(file, need, what)
% The one variable of a .mat file that has the fields in need
v = load(file);
fn = fieldnames(v);
for k = 1:numel(fn)
    if isstruct(v.(fn{k})) && all(isfield(v.(fn{k}), need)), s = v.(fn{k}); return; end
end
error('ieeglab_import_fieldtrip:notFound', '%s holds no %s.', file, what);
end

function elec = local_load_elec(file)
% The variable of a .mat file that is an elec with positions, or holds one in field elec
v = load(file);
fn = fieldnames(v);
for k = 1:numel(fn)
    s = v.(fn{k});
    if ~isempty(local_pos(s)) || (isstruct(s) && isscalar(s) && isfield(s, 'elec') && ~isempty(local_pos(s.elec)))
        elec = local_elec(s, file); return
    end
end
error('ieeglab_import_fieldtrip:noPositions', ['%s holds no electrode positions: a FieldTrip elec (label, ' ...
    'chanpos or elecpos with one row per label), or a structure holding one in field elec.'], file);
end

function elec = local_elec(s, what)
% s itself when it has positions, else its field elec (FieldTrip data, freq)
if ~isempty(local_pos(s)), elec = s; return; end
if isstruct(s) && isscalar(s) && isfield(s, 'elec') && ~isempty(local_pos(s.elec)), elec = s.elec; return; end
error('ieeglab_import_fieldtrip:noPositions', ['%s has no electrode positions: a FieldTrip elec needs label and ' ...
    'chanpos or elecpos with one row per label, or a structure must hold one in field elec.'], what);
end

function pos = local_pos(s)
% Positions with one row per label: chanpos (FieldTrip pairs it with label), else elecpos
pos = [];
if ~isstruct(s) || ~isscalar(s) || ~isfield(s, 'label') || isempty(s.label), return; end
for f = {'chanpos', 'elecpos'}
    if isfield(s, f{1}) && isnumeric(s.(f{1})) && isequal(size(s.(f{1})), [numel(s.label) 3])
        pos = double(s.(f{1})); return
    end
end
end

function [scale, note] = local_mm(elec, pos)
% mm per unit of the positions: elec.unit, else estimated from the size of the
% layout as FieldTrip's ft_estimate_units does (the note says which was assumed)
units = {'m', 1000; 'dm', 100; 'cm', 10; 'mm', 1; 'um', 1e-3};
note = '';
u = '';
if isfield(elec, 'unit') && ~isempty(elec.unit), u = strtrim(char(elec.unit)); end
if isempty(u)
    ok = all(isfinite(pos), 2);
    d = 0;
    if nnz(ok) > 1, d = norm(max(pos(ok, :), [], 1) - min(pos(ok, :), [], 1)); end
    if d > 0
        u = units{min(max(round(log10(d) + 1.8), 1), 4), 1};
        note = sprintf(['elec has no unit: assumed %s, from the size of the layout (diagonal %.3g; a head is ' ...
            'about 0.1 to 0.2 m, as in FieldTrip''s ft_estimate_units). Set elec.unit if this is wrong.'], u, d);
    else
        u = 'mm';
        note = 'elec has no unit, and too few positions to estimate it: assumed mm. Set elec.unit if this is wrong.';
    end
end
k = find(strcmp(regexprep(lower(u), ['^[' char([181 956]) ']'], 'u'), units(:, 1)), 1);
if isempty(k)
    error('ieeglab_import_fieldtrip:unit', 'elec.unit is ''%s'', not a unit of length known here (m, dm, cm, mm, um).', u);
end
scale = units{k, 2};
end

function [u, from] = local_chanunit(data, elec)
% The unit of each channel's values as FieldTrip labels it ('' when not given),
% from data.chanunit, else data.elec, data.hdr or elec, matched by label
u = repmat({''}, numel(data.label), 1); from = '';
cand = {data, 'data'};
if isfield(data, 'elec'), cand(end+1, :) = {data.elec, 'data.elec'}; end
if isfield(data, 'hdr'), cand(end+1, :) = {data.hdr, 'data.hdr'}; end
if ~isempty(elec), cand(end+1, :) = {elec, 'elec'}; end
for k = 1:size(cand, 1)
    s = cand{k, 1};
    if ~isstruct(s) || ~isscalar(s) || ~isfield(s, 'chanunit') || isempty(s.chanunit), continue; end
    try cu = cellstr(s.chanunit); catch, continue; end
    if isscalar(cu)
        u(:) = cu;
    elseif isfield(s, 'label') && numel(s.label) == numel(cu)
        [ok, loc] = ismember(lower(data.label(:)), lower(s.label(:)));
        if ~any(ok), continue; end
        u(ok) = cu(loc(ok));
    else
        continue
    end
    from = [cand{k, 2} '.chanunit']; return
end
end

function sz = local_size(data, T, nTr)
% Median absolute deviation of each channel around its median, from up to 20000 samples
N = T * nTr;
i = unique(round(linspace(1, N, min(N, 20000))));
k = ceil(i / T); j = i - (k - 1) * T;
Y = zeros(numel(data.label), numel(i));
for tr = unique(k)
    m = k == tr;
    Y(:, m) = double(data.trial{tr}(:, j(m)));
end
sz = median(abs(Y - median(Y, 2, 'omitnan')), 2, 'omitnan');
end

function [gain, msg, warn] = local_to_uv(cu, sz, from)
% Factor to uV per channel: the label's, when the values have the size of iEEG in that unit
volt = {'V', 1e6; 'mV', 1e3; 'uV', 1};
band = [0.5 500];                               % median absolute deviation of iEEG, in uV
gain = ones(numel(cu), 1); msg = {}; warn = {};
key = regexprep(strtrim(cu), ['^[' char([181 956]) ']'], 'u');
for k = 1:size(volt, 1)
    ch = strcmpi(key, volt{k, 1});
    if ~any(ch), continue; end
    s = median(sz(ch), 'omitnan');
    fits = s * [volt{:, 2}] >= band(1) & s * [volt{:, 2}] < band(2);
    if fits(k)
        gain(ch) = volt{k, 2};
        if volt{k, 2} ~= 1
            msg{end+1} = sprintf('%s is %s (%d channels): values converted to uV (x %g).', from, volt{k, 1}, nnz(ch), volt{k, 2}); %#ok<AGROW>
        end
    elseif any(fits)
        warn{end+1} = sprintf(['%s is %s (%d channels), but the values have the size of %s (median absolute ' ...
            'deviation %.3g): kept as they are; the label looks wrong.'], from, volt{k, 1}, nnz(ch), volt{fits, 1}, s); %#ok<AGROW>
    else
        warn{end+1} = sprintf(['%s is %s (%d channels), but the values (median absolute deviation %.3g) do not have ' ...
            'the size of iEEG in V, mV or uV: kept as they are.'], from, volt{k, 1}, nnz(ch), s); %#ok<AGROW>
    end
end
end

function [folder, name] = local_file(file)
% Absolute folder and name of the .mat file that load() read
[folder, base, ext] = fileparts(file);
if isempty(ext), ext = '.mat'; end              % load('x') reads x.mat
name = [base ext];
d = dir(fullfile(folder, name));
if ~isempty(d), folder = d(1).folder; return; end
w = which(name);
if ~isempty(w), folder = fileparts(w); end
end

function s = local_space(cs)
if isempty(cs), s = ''; else, s = sprintf(' (%s space)', cs); end
end
