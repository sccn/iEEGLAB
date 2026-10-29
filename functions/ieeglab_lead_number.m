function [lead, num] = ieeglab_lead_number(names)
% ieeglab_lead_number() - Lead name and contact number of contact labels.
%
% Usage:
%   [lead, num] = ieeglab_lead_number({'RA1', 'RA2', "LA'3", 'ROP 12', 'EKG'})
%   % lead = ["RA" "RA" "LA'" "ROP" ""], num = [1 2 3 12 NaN]
%
% The lead is the label without its trailing number, upper case, without
% spaces or punctuation except the prime: sEEG labels write the electrode of
% the other hemisphere with a prime (A'1 is not on lead A). Labels that do not
% end in a number get lead "" and number NaN.
%
% The one rule for "same lead" used by bipolar montages (ieeglab_montage), the
% neighbour exclusion of the common average (ieeglab_car) and the recognition of
% stimulation sites whose contacts were removed (ieeglab_detect_mode).
%
% Cedric Cannard, iEEGLAB, 2026

names = strrep(strtrim(string(names)), char(8217), '''');     % typographic prime
names = upper(regexprep(names, '[^A-Za-z0-9'']', ''));
lead = strings(size(names)); num = nan(size(names));
tok = regexp(names, '^(.*[A-Z''])(\d+)$', 'tokens', 'once');
if ~iscell(tok), tok = {tok}; end
for i = 1:numel(names)
    if ~isempty(tok{i})
        lead(i) = tok{i}(1); num(i) = str2double(tok{i}(2));
    end
end
end
