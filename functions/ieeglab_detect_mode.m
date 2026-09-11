function [mode, info] = ieeglab_detect_mode(EEG)
% ieeglab_detect_mode() - Decide whether a dataset is CCEP, event-related, or
%                         continuous, from its events and channel labels.
%
% Usage:
%   mode = ieeglab_detect_mode(EEG)
%   [mode, info] = ieeglab_detect_mode(EEG)
%
% Returns one of:
%   'ccep'       - single-pulse electrical stimulation. Event types name a pair
%                  of recording contacts, e.g. 'RA1-RA2'. Methods that assume a
%                  known stimulated pair (CARLA, CRP, N1 detection) apply.
%   'erp'        - events exist but do not name contact pairs: an ordinary
%                  stimulus- or response-locked design.
%   'continuous' - no events.
%
% info is a struct with the evidence behind the decision: n_events,
% n_types, n_pair_types, frac_pair, and example_types.
%
% Why this exists: several methods in this plugin are defined only for CCEP
% data. CARLA in particular selects a reference by finding the channels least
% anticorrelated with the evoked response to a KNOWN stimulated pair; applied
% to a visual-task dataset it is not the published method and its output has no
% established interpretation. Detecting the mode lets the pipeline offer those
% methods only where they mean something.
%
% Cedric Cannard, iEEGLAB, 2026

info = struct('n_events',0, 'n_types',0, 'n_pair_types',0, 'frac_pair',0, ...
              'example_types',{{}});

if ~isfield(EEG,'event') || isempty(EEG.event)
    mode = 'continuous';
    return
end

types = {EEG.event.type};
num = cellfun(@isnumeric, types);
types(num) = cellfun(@num2str, types(num), 'UniformOutput', false);
types = string(types);

% 'boundary' events mark discontinuities, not conditions
types = types(~strcmpi(types, 'boundary'));
if isempty(types)
    mode = 'continuous';
    return
end
info.n_events = numel(types);
uTypes = unique(types);
info.n_types = numel(uTypes);
info.example_types = cellstr(uTypes(1:min(5,end)))';

if ~isfield(EEG,'chanlocs') || isempty(EEG.chanlocs)
    mode = 'erp';
    return
end
labels = upper(string({EEG.chanlocs.labels}));

% A type is a stimulation pair when it splits on a separator into two or more
% tokens that are BOTH recording-channel labels. Requiring both tokens to match
% is what stops ordinary hyphenated condition names ('face-house') registering
% as stimulation sites.
isPair = false(1, numel(uTypes));
for i = 1:numel(uTypes)
    parts = upper(strtrim(ieeglab_site_tokens(uTypes(i), {EEG.chanlocs.labels})));
    if numel(parts) >= 2
        isPair(i) = all(ismember(parts, labels));
    end
end

info.n_pair_types = nnz(isPair);
% Weight by how many EVENTS carry a pair type, not just how many type names do
info.frac_pair = mean(ismember(types, uTypes(isPair)));

if info.frac_pair >= 0.5
    mode = 'ccep';
else
    mode = 'erp';
end
end
