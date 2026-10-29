function M = ieeglab_montage(EEG, src, varargin)
% ieeglab_montage() - Build a montage: bipolar pairs of neighbouring contacts,
%                     or a custom montage read from a file, a table or a
%                     FieldTrip montage. Apply it with ieeglab_apply_montage.
%
% Usage:
%   M = ieeglab_montage(EEG, 'bipolar')
%   M = ieeglab_montage(EEG, 'bipolar', 'spacing_factor', 0)    % no distance check
%   M = ieeglab_montage(EEG, 'my_montage.tsv')                  % custom, from a file
%   M = ieeglab_montage(EEG, T)                                 % custom, from a table
%   M = ieeglab_montage(EEG, ftMontage)                         % FieldTrip montage struct
%
% Bipolar: on each lead (the label without its trailing number: RA1, RA2, ...
% form lead RA), every contact minus the next one, 'RA1-RA2' = RA1 - RA2, as
% FieldTrip's tutorial montage for depth electrodes (Stolk et al., 2018). Only
% consecutive numbers are paired, so a removed contact leaves a gap rather than
% a pair spanning two spacings. Clinician-bad contacts are never paired.
% Contacts with no neighbour, and channels that are not iEEG (scalp EEG, ECG),
% are left out of the bipolar montage.
%
% Options for 'bipolar' (name/value):
%   'spacing_factor'  with electrode coordinates, a pair farther apart than this
%                     times the median distance between neighbouring contacts of
%                     its lead is not formed: on a grid numbered row by row, the
%                     last contact of a row and the first of the next are not
%                     neighbours. Default 2; 0 = no check.
%   'types'           channel types paired. Default {'SEEG','ECOG','DBS','IEEG'}
%                     plus channels with no type.
%
% Custom montage: one row per new channel, columns
%   positive  the contact (required)
%   negative  empty or n/a = the contact as it is; one contact = their
%             difference; several contacts separated by spaces, commas or
%             semicolons = the contact minus their mean (a local average)
%   name      optional; default 'positive-negative', or the contact's name
% read from a .tsv, .csv or .txt file (tab or comma separated, header row) or
% given as a table. 'anode'/'cathode' and 'reference' are accepted for
% positive/negative. A FieldTrip montage (labelold, labelnew, tra) is used as
% it is: new = tra * old.
%
% Output M: labelold (contacts used), labelnew (new channels), tra (new x old
% weights), method ('bipolar' or 'custom'), source, and skipped (a table of the
% contacts or rows left out and why).
%
% Cedric Cannard, iEEGLAB, 2026

if nargin < 2 || isempty(src), src = 'bipolar'; end
labels = string({EEG.chanlocs.labels});
labels = strtrim(labels(:));
if (ischar(src) || isstring(src)) && strcmpi(src, 'bipolar')
    M = local_bipolar(EEG, labels, varargin{:});
    return
end

% ---------------- custom ----------------
source = '';
if isstruct(src) && all(isfield(src, {'labelold', 'labelnew', 'tra'}))
    M = struct('labelold', {cellstr(string(src.labelold(:)'))}, 'labelnew', {cellstr(string(src.labelnew(:)'))}, ...
        'tra', double(src.tra), 'method', 'custom', 'source', 'montage struct', ...
        'skipped', table(strings(0,1), strings(0,1), 'VariableNames', {'item', 'reason'}));
    if ~isequal(size(M.tra), [numel(M.labelnew) numel(M.labelold)])
        error('ieeglab_montage:badTra', 'tra must be numel(labelnew) x numel(labelold) (%d x %d); it is %d x %d.', ...
            numel(M.labelnew), numel(M.labelold), size(M.tra, 1), size(M.tra, 2));
    end
    return
end
if ischar(src) || isstring(src)
    source = char(src);
    if ~isfile(source)
        error('ieeglab_montage:noFile', 'Montage file not found: %s', source);
    end
    T = local_read_table(source);
elseif istable(src)
    T = src;
    source = 'table';
else
    error('ieeglab_montage:badInput', ...
        'Give ''bipolar'', a montage file (.tsv/.csv), a table, or a FieldTrip montage (labelold, labelnew, tra).');
end

vn = lower(strtrim(string(T.Properties.VariableNames)));
pos = local_column(T, vn, {'positive', 'anode', 'channel', 'contact'});
if isempty(pos)
    error('ieeglab_montage:noPositive', ...
        'The montage needs a column named positive (or anode): the contact of each new channel.');
end
neg = local_column(T, vn, {'negative', 'cathode', 'reference', 'ref'});
nam = local_column(T, vn, {'name', 'label', 'labelnew'});
nR = numel(pos);
if isempty(neg), neg = strings(nR, 1); end
if isempty(nam), nam = strings(nR, 1); end

Lu = upper(labels);
used = false(numel(labels), 1);
rows = zeros(nR, numel(labels));
newNames = strings(nR, 1);
keep = true(nR, 1);
skipItem = strings(0, 1); skipWhy = strings(0, 1);
for r = 1:nR
    p = strtrim(pos(r));
    ip = find(Lu == upper(p), 1);
    nTok = local_split_negative(neg(r));
    [inN, iN] = ismember(upper(nTok), Lu);
    if isempty(ip) || ~all(inN)
        missing = [p(isempty(ip)) nTok(~inN)];
        keep(r) = false;
        skipItem(end+1, 1) = local_row_name(nam(r), p, nTok); %#ok<AGROW>
        skipWhy(end+1, 1) = "not in the dataset: " + strjoin(missing, ', '); %#ok<AGROW>
        continue
    end
    rows(r, ip) = 1;
    if ~isempty(iN)
        rows(r, iN) = rows(r, iN) - 1 / numel(iN);
    end
    used([ip; iN(:)]) = true;
    newNames(r) = local_row_name(nam(r), labels(ip), labels(iN)');
end
if ~any(keep)
    error('ieeglab_montage:nothingMatched', 'No row of the montage names contacts of this dataset.');
end
if numel(unique(upper(newNames(keep)))) < nnz(keep)
    error('ieeglab_montage:duplicateNames', 'Two rows of the montage give a channel the same name.');
end
M = struct('labelold', {cellstr(labels(used))'}, 'labelnew', {cellstr(newNames(keep))'}, ...
    'tra', rows(keep, used), 'method', 'custom', 'source', source, ...
    'skipped', table(skipItem, skipWhy, 'VariableNames', {'item', 'reason'}));
end

% ========================== bipolar ==========================

function M = local_bipolar(EEG, labels, varargin)
p = inputParser;
p.addParameter('spacing_factor', 2, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('types', {'SEEG', 'ECOG', 'DBS', 'IEEG'});
p.parse(varargin{:});
o = p.Results;

C = numel(labels);
type = strings(C, 1);
if isfield(EEG.chanlocs, 'type')
    type = upper(strtrim(string(cellfun(@(x) char(string(x)), {EEG.chanlocs.type}, 'UniformOutput', false))));
    type = type(:);
end
bad = false(C, 1);
if isfield(EEG.chanlocs, 'status')
    bad = cellfun(@(x) ~isempty(x) && strcmpi(char(string(x)), 'bad'), {EEG.chanlocs.status})';
end
xyz = nan(C, 3);
if isfield(EEG.chanlocs, 'X')
    for c = 1:C
        cl = EEG.chanlocs(c);
        if ~isempty(cl.X) && ~isempty(cl.Y) && ~isempty(cl.Z), xyz(c, :) = double([cl.X cl.Y cl.Z]); end
    end
end
[lead, num] = local_lead_number(labels);

skipItem = strings(0, 1); skipWhy = strings(0, 1);
iEEGtype = ismember(type, upper(string(o.types))) | type == "" | type == "N/A";
ok = iEEGtype & ~bad & lead ~= "";
for c = find(~ok)'
    skipItem(end+1, 1) = labels(c); %#ok<AGROW>
    if ~iEEGtype(c), skipWhy(end+1, 1) = "channel type " + type(c) + " is not iEEG"; %#ok<AGROW>
    elseif bad(c),   skipWhy(end+1, 1) = "marked bad"; %#ok<AGROW>
    else,            skipWhy(end+1, 1) = "label is not a lead name followed by a number"; %#ok<AGROW>
    end
end

pairs = zeros(0, 2);
paired = false(C, 1);
for L = unique(lead(ok))'
    idx = find(ok & lead == L);
    [~, o2] = sort(num(idx)); idx = idx(o2);
    cand = [idx(1:end-1) idx(2:end)];
    cand = cand(num(cand(:, 2)) == num(cand(:, 1)) + 1, :);    % consecutive numbers only
    if isempty(cand), continue; end
    if o.spacing_factor > 0 && size(cand, 1) >= 2
        d = sqrt(sum((xyz(cand(:, 1), :) - xyz(cand(:, 2), :)).^2, 2));
        med = median(d, 'omitnan');
        far = d > o.spacing_factor * med;                      % NaN distances are kept
        for k = find(far)'
            skipItem(end+1, 1) = labels(cand(k, 1)) + "-" + labels(cand(k, 2)); %#ok<AGROW>
            skipWhy(end+1, 1) = sprintf("%.1f mm apart, more than %g x the lead's median %.1f mm", ...
                d(k), o.spacing_factor, med); %#ok<AGROW>
        end
        cand = cand(~far, :);
    end
    pairs = [pairs; cand]; %#ok<AGROW>
    paired(cand(:)) = true;
end
for c = find(ok & ~paired)'
    skipItem(end+1, 1) = labels(c); %#ok<AGROW>
    skipWhy(end+1, 1) = "no neighbouring contact on its lead"; %#ok<AGROW>
end
if isempty(pairs)
    error('ieeglab_montage:noPairs', ...
        ['No bipolar pair could be formed: contacts need labels that are a lead name followed by a ' ...
         'number (RA1, RA2, ...), with consecutive numbers on the same lead.']);
end

used = unique(pairs(:))';
[~, a] = ismember(pairs(:, 1), used);
[~, b] = ismember(pairs(:, 2), used);
nP = size(pairs, 1);
tra = zeros(nP, numel(used));
tra(sub2ind(size(tra), (1:nP)', a)) = 1;
tra(sub2ind(size(tra), (1:nP)', b)) = -1;
M = struct('labelold', {cellstr(labels(used))'}, ...
    'labelnew', {cellstr(labels(pairs(:, 1)) + "-" + labels(pairs(:, 2)))'}, ...
    'tra', tra, 'method', 'bipolar', 'source', '', ...
    'skipped', table(skipItem, skipWhy, 'VariableNames', {'item', 'reason'}));
end

% ========================== helpers ==========================

function [lead, num] = local_lead_number(names)
% Lead = the label without punctuation and trailing number; same rule as the
% neighbour exclusion of ieeglab_car.
names = upper(regexprep(strtrim(string(names)), '[^A-Za-z0-9]', ''));
lead = strings(size(names)); num = nan(size(names));
tok = regexp(names, '^(.*[A-Z])(\d+)$', 'tokens', 'once');
for i = 1:numel(names)
    if ~isempty(tok{i})
        lead(i) = tok{i}(1); num(i) = str2double(tok{i}(2));
    end
end
end

function T = local_read_table(f)
[~, ~, ext] = fileparts(f);
args = {'FileType', 'text', 'TextType', 'string', 'Delimiter', {'\t', ','}};
if strcmpi(ext, '.tsv'), args{end} = '\t'; elseif strcmpi(ext, '.csv'), args{end} = ','; end
opts = detectImportOptions(f, args{:});
opts = setvartype(opts, opts.VariableNames, 'string');
T = readtable(f, opts);
end

function v = local_column(T, vn, names)
v = [];
k = find(ismember(vn, names), 1);
if isempty(k), return; end
v = strtrim(string(T{:, k}));
v(ismissing(v)) = "";
end

function tok = local_split_negative(s)
s = strtrim(string(s));
if s == "" || any(strcmpi(s, ["n/a", "na", "none", "nan"])), tok = strings(1, 0); return; end
tok = strtrim(split(s, {' ', ',', ';'}))';
tok = tok(tok ~= "");
end

function s = local_row_name(given, p, neg)
s = strtrim(string(given));
if s ~= "", return; end
s = string(p);
if ~isempty(neg), s = strjoin([s string(neg(:))'], '-'); end
end
