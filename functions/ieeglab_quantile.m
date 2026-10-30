function q = ieeglab_quantile(x, p)
% ieeglab_quantile() - Quantile of the values in x, without a toolbox.
%
% Usage:
%   q = ieeglab_quantile(x, p)         % p between 0 and 1
%
% Same definition as MATLAB's quantile and prctile, which need the Statistics
% Toolbox in older releases: the sorted values sit at (i - 0.5) / n, with linear
% interpolation between them. NaN values are ignored.
%
% Cedric Cannard, iEEGLAB, 2026

x = sort(x(~isnan(x(:))));
n = numel(x);
if n == 0, q = NaN; return; end
r = n * p + 0.5;                    % fractional position in the sorted values
k = min(max(floor(r), 1), n);
q = x(k) + max(min(r - k, 1), 0) * (x(min(k + 1, n)) - x(k));
end
