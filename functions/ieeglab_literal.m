function s = ieeglab_literal(opt, exclude)
% ieeglab_literal() - Re-runnable MATLAB literal of an options struct, for
%                     EEGLAB history (com) strings.
%
% Usage:
%   s = ieeglab_literal(opt)
%   s = ieeglab_literal(opt, {'events','plot'})   % fields to leave out
%
% Every iEEGLAB function that returns a history string builds it with this, so
% that eval(com) re-runs exactly the call that was made:
%   - text of any type (char, string scalar) is written as a char literal with
%     embedded quotes doubled ('o''brien');
%   - string arrays and cellstr become {{'a','b'}} - the double braces keep a
%     cell-valued field inside struct();
%   - logical and numeric values up to 1000 elements use mat2str;
%   - tables, structs, function handles and other objects are skipped, since
%     they cannot be written as a literal (derived fields such as the events
%     table are re-read from files named in the options).
%
% Cedric Cannard, iEEGLAB, 2026

if nargin < 2, exclude = {}; end
s = 'struct()';
if isempty(opt) || ~isstruct(opt) || ~isscalar(opt), return; end
f = fieldnames(opt);
parts = {};
for i = 1:numel(f)
    if any(strcmp(f{i}, exclude)), continue; end
    lit = local_value(opt.(f{i}));
    if isempty(lit), continue; end
    parts{end+1} = sprintf('''%s'',%s', f{i}, lit); %#ok<AGROW>
end
s = ['struct(' strjoin(parts, ', ') ')'];
end

function lit = local_value(v)
lit = '';
if isstring(v)
    if isscalar(v), v = char(v); else, v = cellstr(v); end
end
if ischar(v)
    if size(v,1) > 1
        v = cellstr(v);
    else
        lit = ['''' strrep(v, '''', '''''') ''''];
        return
    end
end
if iscell(v)
    if isempty(v), lit = '{{}}'; return; end
    isText = cellfun(@(x) ischar(x) || (isstring(x) && isscalar(x)), v(:));
    if all(isText)
        q = cellfun(@(x) ['''' strrep(char(x), '''', '''''') ''''], v(:)', 'UniformOutput', false);
        lit = ['{{' strjoin(q, ',') '}}'];
    end
    return
end
if isstruct(v) && isscalar(v)
    lit = ieeglab_literal(v);          % nested options, e.g. event_filters
    return
end
if (islogical(v) || isnumeric(v)) && ismatrix(v) && numel(v) <= 1000
    if isempty(v)
        lit = '[]';
    else
        lit = mat2str(v);
    end
end
end
