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
% n_types, n_pair_types, pair_types (the event types that name a pair),
% frac_pair, and example_types.
%
% Why this exists: several methods in this plugin are defined only for CCEP
% data. CARLA in particular selects a reference by finding the channels least
% anticorrelated with the evoked response to a KNOWN stimulated pair; applied
% to a visual-task dataset it is not the published method and its output has no
% established interpretation. Detecting the mode lets the pipeline offer those
% methods only where they mean something.
%
% Cedric Cannard, iEEGLAB, 2026

info = struct('n_events',0, 'n_types',0, 'n_pair_types',0, 'pair_types',strings(1,0), ...
              'frac_pair',0, 'example_types',{{}});

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
% Contacts, not channel labels: after bipolar re-referencing a channel is
% 'RA1-RA2' and the contacts RA1 and RA2 are what the event types name.
contacts = ieeglab_label_contacts(EEG.chanlocs);
labels = upper(contacts(:))';
leads = ieeglab_lead_number(contacts);
leads = leads(leads ~= "");

% A type is a stimulation pair when it splits on a separator into two or more
% tokens that are each a recording contact, or a contact name on one of the
% dataset's leads (a stimulated contact removed as bad: 'RA1-RA2' with RA1 and
% RA2 gone is still a pulse). Requiring every token to qualify is what stops
% ordinary hyphenated condition names ('face-house') registering as
% stimulation sites.
isPair = false(1, numel(uTypes));
for i = 1:numel(uTypes)
    parts = upper(strtrim(ieeglab_site_tokens(uTypes(i), cellstr(contacts))));
    if numel(parts) >= 2
        [pl, pn] = ieeglab_lead_number(parts);
        isPair(i) = all(ismember(parts, labels) | (ismember(pl, leads) & ~isnan(pn)));
    end
end

info.n_pair_types = nnz(isPair);
info.pair_types = uTypes(isPair);
% Weight by how many EVENTS carry a pair type, not just how many type names do
info.frac_pair = mean(ismember(types, uTypes(isPair)));

if info.frac_pair >= 0.5
    mode = 'ccep';
else
    mode = 'erp';
end
end
