function [site, stimIdx] = ieeglab_canonical_site(types, labels)
% ieeglab_canonical_site() - Order-independent stimulation-site name for each
%                            event type, and the stimulated contacts it names.
%
% Usage:
%   [site, stimIdx] = ieeglab_canonical_site({EEG.event.type}, {EEG.chanlocs.labels})
%
% For a type that names a stimulated pair - two or more tokens that all look
% like contact labels, at least one of them in the montage - the site is the
% tokens in the montage's own spelling, sorted and joined with '-': 'ROP4-ROP2'
% and 'rop2-ROP4' both become 'ROP2-ROP4'. Any other type (a condition name,
% 'boundary') is returned unchanged, so ordinary hyphenated conditions such as
% 'house-face' are never reordered or merged.
%
% Labels are compared trimmed and case-insensitively, with a second pass that
% ignores punctuation (sEEG prime conventions: LA'3 vs LA3).
%
% Outputs:
%   site    - string array, one per type
%   stimIdx - cell array of column vectors: montage indices of the named contacts
%
% Every step that needs a site - rare-condition removal, re-referencing, N1,
% CRP, the connectivity matrix, trial rejection - goes through this function.
%
% Cedric Cannard, iEEGLAB, 2026

if ischar(types) || isstring(types), types = cellstr(types); end
n = numel(types);
site = strings(1, n);
stimIdx = repmat({zeros(0,1)}, 1, n);
labels = cellstr(string(labels));
L = upper(strtrim(string(labels)));
Lbare = regexprep(L, '[^\w]', '');
contactRe = '^[A-Za-z][A-Za-z'' ]*\d+$';

for i = 1:n
    t = types{i};
    if isnumeric(t) || islogical(t), t = num2str(t); end
    t = strtrim(string(t));
    if ismissing(t) || t == "", continue; end
    tok = ieeglab_site_tokens(t, labels);
    if isempty(tok), site(i) = t; continue; end

    ut = upper(strtrim(tok));
    [tf, loc] = ismember(ut, L);
    if ~all(tf)
        [tf2, loc2] = ismember(regexprep(ut, '[^\w]', ''), Lbare);
        loc(~tf & tf2) = loc2(~tf & tf2);
        tf = tf | tf2;
    end
    v = unique(loc(tf));
    stimIdx{i} = v(:);

    looksLikeContact = tf | ~cellfun(@isempty, regexp(cellstr(tok), contactRe, 'once'));
    if numel(tok) >= 2 && all(looksLikeContact) && any(tf)
        canon = tok;
        canon(tf) = strtrim(string(labels(loc(tf))));   % the montage's own spelling
        site(i) = strjoin(sort(canon), '-');
    else
        site(i) = t;
    end
end
end
