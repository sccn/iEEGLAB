function [sites, stimIdx] = ieeglab_epoch_sites(EEG)
% ieeglab_epoch_sites() - Stimulation site and stimulated contacts per epoch.
%
% Usage:
%   [sites, stimIdx] = ieeglab_epoch_sites(EEG)
%
% Outputs:
%   sites   - 1 x trials string. The canonical, order-independent site name
%             ('ROP2-ROP4' whether recorded as ROP2-ROP4 or ROP4-ROP2) for CCEP
%             epochs; the raw event type otherwise; "" when unknown.
%   stimIdx - 1 x trials cell of column vectors: indices of the stimulated
%             contacts that are present in the montage.
%
% The site is read from EEG.epoch(i).eventtype - the event the epoch is locked
% to, i.e. the one at latency 0 - which pop_epoch always populates and which
% is aligned to the data by construction. An event type is treated as a
% stimulation pair only when it has two or more tokens, all of which look like
% contact labels and at least one of which is in the montage. That keeps
% ordinary hyphenated condition names ('house-face') from being reordered and
% merged, while still naming a site correctly after one of its contacts has
% been removed as a bad channel.
%
% Cedric Cannard, iEEGLAB, 2026

N = EEG.trials;
sites   = strings(1, N);
stimIdx = repmat({zeros(0,1)}, 1, N);
if N < 1 || ~isfield(EEG,'epoch') || isempty(EEG.epoch) || numel(EEG.epoch) ~= N
    return
end

labels     = upper(strtrim(string({EEG.chanlocs.labels})));
labelsBare = regexprep(labels, '[^\w]', '');
contactRe  = '^[A-Za-z][A-Za-z'']*\d+$';

for i = 1:N
    ty = local_epoch_type(EEG.epoch(i));
    if ty == "", continue; end
    tok = ieeglab_site_tokens(ty);
    if isempty(tok), sites(i) = ty; continue; end

    utok = upper(tok);
    [tf, loc] = ismember(utok, labels);
    if ~any(tf)
        % Second chance: ignore punctuation differences such as LA'3 vs LA3
        [tf, loc] = ismember(regexprep(utok, '[^\w]', ''), labelsBare);
    end
    v = unique(loc(tf));
    stimIdx{i} = v(:);

    looksLikeContact = tf | ~cellfun(@isempty, regexp(cellstr(tok), contactRe, 'once'));
    if numel(tok) >= 2 && all(looksLikeContact) && any(tf)
        sites(i) = strjoin(sort(tok), '-');
    else
        sites(i) = ty;
    end
end
end

function ty = local_epoch_type(ep)
% Event type at latency 0 for this epoch. pop_epoch stores eventtype as a cell
% when several events fall inside the epoch window.
ty = "";
if ~isfield(ep,'eventtype') || isempty(ep.eventtype), return; end
t = ep.eventtype;
if iscell(t)
    k = 1;
    if isfield(ep,'eventlatency') && iscell(ep.eventlatency) && numel(ep.eventlatency) == numel(t)
        [~, k] = min(cellfun(@(x) abs(double(x(1))), ep.eventlatency));
    end
    t = t{k};
end
if isnumeric(t), t = num2str(t); end
ty = strtrim(string(t));
end
