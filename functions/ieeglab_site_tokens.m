function tok = ieeglab_site_tokens(val)
% ieeglab_site_tokens() - Split a stimulation-site string into contact tokens.
%
% Usage:
%   tok = ieeglab_site_tokens('RA1-RA2')      % -> ["RA1" "RA2"]
%   tok = ieeglab_site_tokens("LA'3 / LA'4")  % -> ["LA'3" "LA'4"]
%
% Separators are - + / | , ; and whitespace. Apostrophes are kept because
% sEEG labels such as LA'3 use them. BIDS missing markers (n/a, nan, none,
% undefined, missing) are removed BEFORE splitting: splitting first turned
% 'n/a' into the two tokens 'n' and 'a'.
%
% This is the single tokenizer used everywhere a stimulation site is parsed,
% so the load, preprocessing, referencing and statistics steps all agree on
% what a site is.
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

s = char(strjoin(string(val(:))', ' '));
s = regexprep(s, '(?i)n/a', ' ');                 % drop BIDS n/a before splitting on '/'
parts = regexp(s, '[-+/|,;\s]+', 'split');
parts = strtrim(string(parts));
parts = parts(parts ~= "");
if isempty(parts), return; end
isMissing = ~cellfun(@isempty, regexpi(cellstr(parts), '^(nan|none|undefined|missing|\?)$', 'once'));
tok = parts(~isMissing);
end
