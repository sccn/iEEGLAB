function mask = ieeglab_events_using(types, contacts)
% ieeglab_events_using() - Which events stimulate any of the given contacts.
%
% Usage:
%   mask = ieeglab_events_using(EEG.event, {'RA1','ROP3'})
%   mask = ieeglab_events_using({EEG.event.type}, removedLabels)
%
% Returns a logical row, one entry per event, true when the event's
% stimulation site (e.g. 'RA1-RA2') includes one of the contacts by EXACT
% token match, case-insensitive.
%
% This replaces contains(types, contacts), which several functions used and
% which is a SUBSTRING test: removing contact 'RA1' also matched 'RA10-RA9'
% and silently dropped those trials. Only events naming at least two contacts
% (a stimulation pair) can match, so an ordinary condition name that happens
% to equal a channel label is never touched.
%
% Cedric Cannard, iEEGLAB, 2026

if isstruct(types)
    if isempty(types) || ~isfield(types,'type'), mask = false(1,0); return; end
    types = {types.type};
end
if ischar(types) || isstring(types), types = cellstr(types); end
mask = false(1, numel(types));
if isempty(types) || isempty(contacts), return; end

if ischar(contacts), contacts = {contacts}; end
c = upper(strtrim(string(contacts(:)')));
c = c(c ~= "");
if isempty(c), return; end

for i = 1:numel(types)
    t = types{i};
    if isnumeric(t), t = num2str(t); end
    tok = upper(ieeglab_site_tokens(t));
    mask(i) = numel(tok) >= 2 && any(ismember(tok, c));
end
end
