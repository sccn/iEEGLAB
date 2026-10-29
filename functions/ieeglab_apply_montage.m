function EEG = ieeglab_apply_montage(EEG, M, verbose)
% ieeglab_apply_montage() - Re-reference to a montage built by ieeglab_montage:
%                           new channels = M.tra * contacts.
%
% Usage:
%   EEG = ieeglab_apply_montage(EEG, ieeglab_montage(EEG, 'bipolar'))
%   EEG = ieeglab_apply_montage(EEG, M, false)          % quiet
%
% Works on continuous and epoched data. A montage row that uses a contact
% missing from the dataset is dropped with a warning; contacts of the dataset
% that no row uses are left out of the result.
%
% Each new channel keeps the channel information of its first positive
% contact, with:
%   labels   the montage's name ('RA1-RA2')
%   X, Y, Z  the midpoint of the two contacts of a bipolar pair (as FieldTrip),
%            otherwise the mean of the positive contacts; empty without
%            coordinates. Spherical fields, which iEEG coordinates do not use,
%            are emptied rather than left describing the old contact.
%   status   'bad' if any contact it uses is marked bad
%   ref      'bipolar' or 'montage', as EEG.ref (the contacts subtracted from
%            each channel are in EEG.ieeglab.montage)
%   contacts the contacts it uses, positive first ('RA3, RA2, RA4')
%
% Results computed on the old channels (CCEP statistics, connectivity matrix,
% spectra) and any ICA decomposition no longer describe the data and are
% removed. EEG.ieeglab.montage keeps the montage, restricted to the contacts
% used; stimulation sites resolve through chanlocs.contacts
% (ieeglab_label_contacts), so a channel that uses a stimulated contact is
% treated as stimulated, whatever its name.
%
% Cedric Cannard, iEEGLAB, 2026

if nargin < 3, verbose = true; end
labels = strtrim(string({EEG.chanlocs.labels}));
[found, idx] = ismember(upper(strtrim(string(M.labelold))), upper(labels));
tra = double(M.tra);
dropRow = any(tra(:, ~found) ~= 0, 2) | ~any(tra ~= 0, 2);
if any(dropRow)
    warning('ieeglab_apply_montage:missingContacts', ...
        '%d montage channel(s) use contacts not in the dataset (or none) and are left out: %s', ...
        nnz(dropRow), strjoin(string(M.labelnew(dropRow)), ', '));
end
tra = tra(~dropRow, :);
newLabels = string(M.labelnew(~dropRow));
if isempty(newLabels)
    error('ieeglab_apply_montage:empty', 'No channel of the montage can be formed from this dataset.');
end

% A: (n_new x n_chan), the montage in the dataset's channel order
C = EEG.nbchan;
A = zeros(numel(newLabels), C);
A(:, idx(found)) = tra(:, found);
inUse = any(A ~= 0, 1);

% ---------------- data: one product for all samples and epochs ----------------
% X: (n_chan x n_time*n_epochs) -> (n_new x n_time*n_epochs); the data keep their class
sz = size(EEG.data);
X = reshape(EEG.data, C, []);
Y = cast(A, class(X)) * X;
EEG.data = reshape(Y, [numel(newLabels) sz(2:end)]);

% ---------------- channel information ----------------
old = EEG.chanlocs;
isBad = false(1, C);
if isfield(old, 'status')
    isBad = cellfun(@(x) ~isempty(x) && strcmpi(char(string(x)), 'bad'), {old.status});
end
[oldContacts, oldOwner] = ieeglab_label_contacts(old);   % contacts behind the old channels
spherical = intersect(fieldnames(old), {'theta', 'radius', 'sph_theta', 'sph_phi', 'sph_radius', 'sph_theta_besa', 'sph_phi_besa'});
if strcmp(M.method, 'bipolar'), refName = 'bipolar'; else, refName = 'montage'; end
for r = numel(newLabels):-1:1
    pos = find(A(r, :) > 0);
    neg = find(A(r, :) < 0);
    if isempty(pos), pos = neg(1); neg = neg(2:end); end   % a row of negative weights only
    cl = old(pos(1));
    cl.labels = char(newLabels(r));
    xyzP = local_xyz(old(pos));
    xyzN = local_xyz(old(neg));
    if isscalar(pos) && isscalar(neg) && all(isfinite([xyzP xyzN]))
        xyz = (xyzP + xyzN) / 2;                     % bipolar pair: midpoint
    else
        xyz = mean(xyzP, 1);
    end
    if all(isfinite(xyz))
        cl.X = xyz(1); cl.Y = xyz(2); cl.Z = xyz(3);
    elseif isfield(cl, 'X')
        cl.X = []; cl.Y = []; cl.Z = [];
    end
    for f = spherical(:)'
        cl.(f{1}) = [];
    end
    if isfield(cl, 'status')
        if any(isBad([pos neg])), cl.status = 'bad'; else, cl.status = 'good'; end
    end
    if isfield(cl, 'urchan'), cl.urchan = []; end
    cl.ref = refName;                                % eeg_checkset copies chanlocs(1).ref into EEG.ref
    cs = strings(0, 1);
    for w = [pos neg]
        cs = [cs; oldContacts(oldOwner == w)]; %#ok<AGROW>
    end
    cl.contacts = char(strjoin(unique(cs, 'stable'), ', '));   % stimulation sites resolve through these
    newLocs(r) = cl; %#ok<AGROW>
end
EEG.chanlocs = newLocs;
EEG.nbchan = numel(newLabels);

% ---------------- what no longer describes the data ----------------
for f = {'icaweights', 'icasphere', 'icawinv', 'icachansind', 'icaact'}
    if isfield(EEG, f{1}), EEG.(f{1}) = []; end
end
if ~isfield(EEG, 'ieeglab') || ~isstruct(EEG.ieeglab), EEG.ieeglab = struct(); end
stale = intersect(fieldnames(EEG.ieeglab), {'stats', 'ccep_matrix', 'spectral', 'baseline_not_flat'});
if ~isempty(stale), EEG.ieeglab = rmfield(EEG.ieeglab, stale); end

M.labelold = cellstr(labels(inUse));
M.labelnew = cellstr(newLabels)';
M.tra = A(:, inUse);
EEG.ieeglab.montage = M;
EEG.ref = refName;

if verbose
    fprintf('[montage] %s: %d contacts -> %d channels.\n', M.method, nnz(inUse), EEG.nbchan);
    notUsed = labels(~inUse);
    if ~isempty(notUsed)
        fprintf('[montage] Left out (%d): %s\n', numel(notUsed), strjoin(notUsed, ', '));
    end
    if isfield(M, 'skipped') && height(M.skipped) > 0
        for k = 1:height(M.skipped)
            fprintf('[montage]   %-16s %s\n', M.skipped.item(k), M.skipped.reason(k));
        end
    end
end
end

function xyz = local_xyz(cl)
% (n x 3) coordinates; NaN rows for contacts without them
xyz = nan(numel(cl), 3);
if ~isfield(cl, 'X'), return; end
for k = 1:numel(cl)
    if ~isempty(cl(k).X) && ~isempty(cl(k).Y) && ~isempty(cl(k).Z)
        xyz(k, :) = double([cl(k).X cl(k).Y cl(k).Z]);
    end
end
end
