function tok = ieeglab_site_tokens(val, labels)
% ieeglab_site_tokens() - Split a stimulation-site string into contact tokens.
%
% Usage:
%   tok = ieeglab_site_tokens('RA1-RA2')                  % -> ["RA1" "RA2"]
%   tok = ieeglab_site_tokens("LA'3 / LA'4")              % -> ["LA'3" "LA'4"]
%   tok = ieeglab_site_tokens('ROP 1-ROP 2', labels)      % -> ["ROP 1" "ROP 2"]
%
% The string is split on the pair separators - + / | , ; first. A piece that
% is itself a montage label (when labels are given) is kept whole, so labels
% with internal spaces such as 'ROP 1' survive; any other piece is further
% split on whitespace. Apostrophes are kept (sEEG labels such as LA'3).
% BIDS missing markers (n/a, nan, none, undefined, missing) are removed before
% splitting: splitting first turned 'n/a' into the tokens 'n' and 'a'.
%
% This is the single tokenizer used wherever a stimulation site is parsed, so
% loading, preprocessing, re-referencing and statistics agree on what a site is.
%
% Cedric Cannard, iEEGLAB, 2026

tok = strings(1,0);
if nargin < 1 || isempty(val), return; end
if istable(val), val = val{1,1}; end
% Convert char BEFORE any (:) indexing: val(:) on a char vector yields one
% character per element, which turned 'RA1-RA2' into the tokens R A 1 R A 2.
if ischar(val) || iscategorical(val) || iscell(val), val = string(val); end
if isnumeric(val) || islogical(val)
    val = double(val(:)');
    val = val(~isnan(val));
    tok = string(val);
    return
end

s = char(strjoin(val(:)', ' '));
s = regexprep(s, '(?i)(^|\s|[-+/|,;])n/a(?=$|\s|[-+/|,;])', '$1 ');   % BIDS n/a, before splitting on '/'
pieces = regexp(s, '[-+/|,;]+', 'split');
pieces = strtrim(string(pieces));
pieces = pieces(pieces ~= "");
if isempty(pieces), return; end

L = strings(0,1);
if nargin >= 2 && ~isempty(labels)
    L = upper(strtrim(string(labels(:))));
end
for k = 1:numel(pieces)
    pc = pieces(k);
    if ~isempty(L) && any(upper(pc) == L)
        tok(end+1) = pc; %#ok<AGROW>
    else
        sub = regexp(char(pc), '\s+', 'split');
        sub = string(sub(~cellfun(@isempty, sub)));
        tok = [tok, sub(:)']; %#ok<AGROW>
    end
end
if isempty(tok), return; end
isMissing = ~cellfun(@isempty, regexpi(cellstr(tok), '^(nan|none|undefined|missing|n/?a|\?)$', 'once'));
tok = tok(~isMissing);
end
