function [sites, stimIdx] = ieeglab_epoch_sites(EEG)
% ieeglab_epoch_sites() - Stimulation site and stimulated contacts per epoch.
%
% Usage:
%   [sites, stimIdx] = ieeglab_epoch_sites(EEG)
%
% Outputs:
%   sites   - 1 x trials string. The canonical, order-independent site name
%             (see ieeglab_canonical_site) for CCEP epochs; the raw event type
%             otherwise; "" when unknown.
%   stimIdx - 1 x trials cell of column vectors: indices of the stimulated
%             contacts present in the montage.
%
% The site is read from the event each epoch is locked to (latency 0), which
% pop_epoch always records in EEG.epoch and which is aligned to the data by
% construction. N1, CRP, re-referencing, the connectivity matrix and trial
% rejection all take the stimulated contacts from here, so they agree even when
% labels differ in case, whitespace or punctuation from the event names.
%
% Cedric Cannard, iEEGLAB, 2026

N = EEG.trials;
sites   = strings(1, N);
stimIdx = repmat({zeros(0,1)}, 1, N);
if N < 1 || ~isfield(EEG,'epoch') || isempty(EEG.epoch) || numel(EEG.epoch) ~= N
    return
end
types = cell(1, N);
for i = 1:N
    types{i} = char(local_epoch_type(EEG.epoch(i)));
end
[sites, stimIdx] = ieeglab_canonical_site(types, {EEG.chanlocs.labels});
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
