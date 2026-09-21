function EEG = get_elec_coor(EEG, elecs, opt)
% get_elec_coor() - Copy electrode XYZ coordinates from a BIDS electrodes table
%                   onto EEG.chanlocs, matching by label.
%
% Usage:
%   EEG = get_elec_coor(EEG, elecs)
%   EEG = get_elec_coor(EEG, elecs, opt)
%
% opt.allow_positional_match (default false)
%   If no labels match, fall back to assigning coordinates by ROW ORDER.
%   Off by default because a wrong assumption here mislabels every channel and
%   still reports a 100% match. Only enable when you have verified that TSV row
%   order corresponds to channel order.
%GET_ELEC_COOR  Copy X/Y/Z from table into EEG.chanlocs with robust matching.
%
% Attempts:
%   1) Direct case-insensitive exact label match.
%   2) If 0% matched: normalize BOTH EEG and TSV labels by removing "EEG"/"eeg"
%      and tidying separators (without mutating visible labels), then retry.
%   3) If still 0% matched: overwrite EEG.chanlocs(i).labels from TSV order, then retry.
%
% Also:
%   - Converts X/Y/Z from strings/cells to doubles safely.
%   - Prints brief sanity previews (label changes, coordinates).

assert(istable(elecs), 'elecs must be a table');

% ---------- Column indices (case-insensitive) ----------
ilab = 1;
vars = elecs.Properties.VariableNames;
ix = find(strcmpi(vars,'x'), 1);
iy = find(strcmpi(vars,'y'), 1);
iz = find(strcmpi(vars,'z'), 1);
assert(~isempty(ix) && ~isempty(iy) && ~isempty(iz), ...
    'The electrodes table must have X, Y, Z columns (case-insensitive).');

% TSV label lookup
tsv_labels    = string(elecs{:, ilab});
tsv_labels_lc = lower(strtrim(tsv_labels));

% ---------- Attempt 1: direct label match ----------
if nargin < 3 || isempty(opt), opt = struct(); end

[EEG, matchedMask] = copy_xyz_if_matched(EEG, elecs, tsv_labels, tsv_labels_lc, ix, iy, iz);
report_match('[get_elec_coor] Attempt 1 (as-is)', EEG, matchedMask);

% ---------- Attempt 2: normalized comparison (remove "EEG"/tidy) w/o relabel ----------
if sum(matchedMask) == 0
    eeg_labels        = string({EEG.chanlocs.labels}.');
    eeg_labels_norm   = normalize_label(eeg_labels);        % remove "EEG", tidy
    tsv_labels_norm   = normalize_label(tsv_labels);        % same transform on TSV
    tsv_labels_norm_l = lower(strtrim(tsv_labels_norm));

    [EEG, matchedMask] = copy_xyz_if_matched_norm( ...
        EEG, elecs, tsv_labels, tsv_labels_norm_l, eeg_labels_norm, ix, iy, iz);

    report_match('[get_elec_coor] Attempt 2 (normalized "no-EEG")', EEG, matchedMask);
end

% ---------- Attempt 3: overwrite labels from TSV order and retry ----------
% DANGEROUS and therefore opt-in. This assumes TSV row order equals channel
% order. When that assumption is wrong it attaches the wrong name AND the wrong
% coordinates to every channel, and then reports "100% matched" - a silently
% wrong result that looks like a success. It is now gated behind
% opt.allow_positional_match and refuses outright unless the counts agree.
if sum(matchedMask) == 0
    allowPositional = isstruct(opt) && isfield(opt,'allow_positional_match') && opt.allow_positional_match;

    if ~allowPositional
        error('get_elec_coor:noLabelMatch', ...
            ['None of the %d electrode names in the .tsv match the %d channel labels ' ...
             'in the dataset, so coordinates cannot be assigned.\n\n' ...
             'TSV names (first 5):  %s\n' ...
             'Dataset labels (first 5): %s\n\n' ...
             'Fix the labels so they correspond. Only if you are certain that TSV row ' ...
             'order matches channel order, re-run with ' ...
             'get_elec_coor(EEG, elecs, struct(''allow_positional_match'', true)) - ' ...
             'this assigns coordinates by position and will silently mislabel every ' ...
             'channel if the assumption is wrong.'], ...
             height(elecs), EEG.nbchan, ...
             strjoin(cellstr(tsv_labels(1:min(5,end)))', ', '), ...
             strjoin({EEG.chanlocs(1:min(5,end)).labels}, ', '));
    end

    if height(elecs) ~= EEG.nbchan
        error('get_elec_coor:positionalCountMismatch', ...
            ['Positional matching was requested but the .tsv has %d rows and the dataset ' ...
             'has %d channels. Positional matching is only defensible when the counts ' ...
             'are equal.'], height(elecs), EEG.nbchan);
    end

    warning('get_elec_coor:positionalMatch', ...
        ['Assigning electrode coordinates BY ROW ORDER, not by label. Verify the result ' ...
         'against the anatomy before trusting any analysis built on it.']);

    oldLabels = string({EEG.chanlocs.labels}.');
    m = min(EEG.nbchan, height(elecs));

    % Overwrite labels in order
    for i = 1:m
        EEG.chanlocs(i).labels = char(tsv_labels(i));
    end

    % Sanity print: first 20 (or all if <20)
    k = min(20, m);
    fprintf('[get_elec_coor] Attempt 3 relabel: first %d channel(s) from TSV order (old -> new):\n', k);
    for i = 1:k
        fprintf('  [%3d] "%s" -> "%s"\n', i, oldLabels(i), string(EEG.chanlocs(i).labels));
    end
    if m > k
        fprintf('  ... and %d more relabeled.\n', m - k);
    end

    % Re-attempt coordinate copy after relabel
    [EEG, matchedMask] = copy_xyz_if_matched(EEG, elecs, tsv_labels, tsv_labels_lc, ix, iy, iz);
    report_match('[get_elec_coor] After Attempt 3 (TSV relabel)', EEG, matchedMask);
end

% ---------- Report channels left without a position ----------
noXYZ = arrayfun(@(ch) ~isfield(ch,'X') || isempty(ch.X), EEG.chanlocs);
if any(noXYZ) && ~all(noXYZ)
    warning('get_elec_coor:unpositioned', ...
        ['%d of %d channels have no position in the electrodes file and are left without ' ...
         'coordinates: %s. Preprocessing can remove them (remove_no_coords).'], ...
        nnz(noXYZ), numel(noXYZ), strjoin({EEG.chanlocs(noXYZ).labels}, ', '));
end

% ---------- Coordinate sanity preview ----------
preview_n = min(5, EEG.nbchan);
fprintf('[get_elec_coor] First %d coordinate rows after assignment:\n', preview_n);
for ii = 1:preview_n
    c = EEG.chanlocs(ii);
    fprintf('  [%3d] %-12s  X=%s  Y=%s  Z=%s\n', ii, c.labels, ...
        num2str_safe(c, 'X'), num2str_safe(c, 'Y'), num2str_safe(c, 'Z'));
end

end

% ---- helpers ----
function [EEG, matchedMask] = copy_xyz_if_matched(EEG, elecs, tsv_labels, tsv_labels_lc, ix, iy, iz)
keys = lower(strtrim(string({EEG.chanlocs.labels}.')));
[EEG, matchedMask] = assign_xyz(EEG, elecs, tsv_labels, tsv_labels_lc, keys, ix, iy, iz);
end

function [EEG, matchedMask] = copy_xyz_if_matched_norm(EEG, elecs, tsv_labels_orig, tsv_labels_norm_l, eeg_labels_norm, ix, iy, iz)
% Compare normalized labels, assign XYZ from the original TSV row. Does NOT
% change visible labels.
keys = lower(strtrim(string(eeg_labels_norm)));
[EEG, matchedMask] = assign_xyz(EEG, elecs, tsv_labels_orig, tsv_labels_norm_l, keys, ix, iy, iz);
end

function [EEG, matchedMask] = assign_xyz(EEG, elecs, tsv_labels, tsv_keys, keys, ix, iy, iz)
% Every channel's position is cleared first: a channel the TSV does not
% position must not keep coordinates from an earlier load, an importer or a
% template, possibly in another space. A row counts only when all three
% coordinates are finite (no half-filled positions). With duplicate names the
% first COMPLETE row is used and the duplication is reported.
matchedMask = false(EEG.nbchan,1);
dupNames = {}; incomplete = {};
for i = 1:EEG.nbchan
    EEG.chanlocs(i).X = []; EEG.chanlocs(i).Y = []; EEG.chanlocs(i).Z = [];
    hits = find(tsv_keys == keys(i));
    if isempty(hits), continue; end
    if numel(hits) > 1, dupNames{end+1} = char(tsv_labels(hits(1))); end %#ok<AGROW>
    chosen = 0;
    for h = hits(:)'
        xyz = [toNumericScalar(elecs{h, ix}), toNumericScalar(elecs{h, iy}), toNumericScalar(elecs{h, iz})];
        if all(isfinite(xyz)), chosen = h; break; end
    end
    if chosen == 0
        incomplete{end+1} = char(EEG.chanlocs(i).labels); %#ok<AGROW>
        continue
    end
    EEG.chanlocs(i).matched_elec_label = char(tsv_labels(chosen));
    EEG.chanlocs(i).X = xyz(1); EEG.chanlocs(i).Y = xyz(2); EEG.chanlocs(i).Z = xyz(3);
    matchedMask(i) = true;
end
if ~isempty(dupNames)
    warning('get_elec_coor:duplicateNames', ...
        'The electrodes file lists these names more than once; the first complete row was used: %s', ...
        strjoin(unique(dupNames), ', '));
end
if ~isempty(incomplete)
    warning('get_elec_coor:incompleteCoordinates', ...
        'These channels have a row but not three finite coordinates, so they get no position: %s', ...
        strjoin(incomplete, ', '));
end
end

function lab = normalize_label(lab)
% Remove any "EEG"/"eeg" substring and tidy separators; keeps original case of others.
    lab = regexprep(string(lab), '(?i)eeg', '');   % strip "EEG" case-insensitive
    lab = strip_separators(lab);                   % trim ends, collapse runs
end

function out = strip_separators(lbl)
% Trim leading/trailing separators/spaces and collapse repeats to a single underscore.
    out = regexprep(string(lbl), '(^[\s_\-]+|[\s_\-]+$)', '');
    out = regexprep(out, '[\s_\-]{2,}', '_');
    out = strtrim(out);
end

function report_match(tag, EEG, matchedMask)
nTotal     = EEG.nbchan;
nMatched   = sum(matchedMask);
nUnmatched = nTotal - nMatched;
pctMatched = 100 * nMatched / max(nTotal,1);
fprintf('%s: Matched %d / %d (%.1f%%). Unmatched: %d (%.1f%%).\n', ...
    tag, nMatched, nTotal, pctMatched, nUnmatched, 100 - pctMatched);
end

function val = toNumericScalar(x)
% Returns a double scalar or NaN. Handles cell, string, char, missing.
    if iscell(x), x = x{1}; end
    if isstring(x), x = char(x); end
    if ischar(x)
        x = strrep(x, ',', '.');    % allow comma decimals
        x = str2double(strtrim(x));
    end
    if isnumeric(x) && isscalar(x)
        val = double(x);
    else
        val = NaN;
    end
end

function s = num2str_safe(c, fld)
    if isfield(c, fld) && ~isempty(c.(fld))
        s = num2str(c.(fld));
    else
        s = '[]';
    end
end
