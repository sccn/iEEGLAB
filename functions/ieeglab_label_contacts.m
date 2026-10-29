function [contacts, owner] = ieeglab_label_contacts(labels)
% ieeglab_label_contacts() - The recording contacts behind each channel label.
%
% Usage:
%   [contacts, owner] = ieeglab_label_contacts(EEG.chanlocs)
%   [contacts, owner] = ieeglab_label_contacts({EEG.chanlocs.labels})
%
% A channel made by ieeglab_apply_montage lists its contacts in
% chanlocs.contacts ('RA3, RA2, RA4'), whatever its name; give the chanlocs to
% use them. Otherwise, a label made of two or more contact names joined by '-'
% (a bipolar or other derived channel, e.g. 'RA1-RA2') stands for those
% contacts, and any other label stands for itself. A piece counts as a
% contact name when it is letters followed by a number ('RA1', "LA'3",
% 'ROP 2'), so labels such as 'EKG-1' or 'C3-Ref' are left whole.
%
% Outputs:
%   contacts - string column, every contact named by the labels, in label order
%   owner    - column of the same length: the label (channel) each contact
%              belongs to
%
% On a montage of single contacts, contacts equals the labels and owner is
% 1:n, so code that matches stimulation sites against contacts behaves exactly
% as before; on a bipolar montage, a stimulated contact reaches every channel
% that uses it.
%
% Cedric Cannard, iEEGLAB, 2026

known = {};
if isstruct(labels)                                  % chanlocs (or a dataset)
    if isfield(labels, 'chanlocs'), labels = labels.chanlocs; end
    if isfield(labels, 'contacts'), known = {labels.contacts}; end
    labels = {labels.labels};
end
labels = strtrim(string(labels(:)));
n = numel(labels);
contacts = strings(0, 1);
owner = zeros(0, 1);
contactRe = '^[A-Za-z][A-Za-z'' ]*\d+$';
for i = 1:n
    if ~isempty(known) && ~isempty(known{i})
        parts = strtrim(split(string(known{i}), ','));
        parts = parts(parts ~= "");
        contacts = [contacts; parts(:)]; %#ok<AGROW>
        owner = [owner; repmat(i, numel(parts), 1)]; %#ok<AGROW>
        continue
    end
    parts = strtrim(split(labels(i), '-'));
    if numel(parts) >= 2 && all(~cellfun(@isempty, regexp(cellstr(parts), contactRe, 'once')))
        contacts = [contacts; parts(:)]; %#ok<AGROW>
        owner = [owner; repmat(i, numel(parts), 1)]; %#ok<AGROW>
    else
        contacts(end+1, 1) = labels(i); %#ok<AGROW>
        owner(end+1, 1) = i; %#ok<AGROW>
    end
end
end
