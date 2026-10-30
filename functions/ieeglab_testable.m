function [ok, msg] = ieeglab_testable(p, pMin, group, correct, alpha)
% ieeglab_testable() - Which pairs can reach significance at all.
%
% Usage:
%   [ok, msg] = ieeglab_testable(p, pMin, group, correct, alpha)
%
% A sign-flip test has a smallest possible p: 2/2^K (or 1/2^K for a one-sided
% N1) when all sign patterns of K trials are used, 1/(n_perm + 1) when they are
% sampled. With few trials, or a large family of tests, that smallest p can
% stay above what the correction accepts, and the pair cannot come out
% significant whatever the data. Such a pair has no result: "no response" would
% read as a finding.
%
% Inputs (n x 1 each, except the last two):
%   p        p-values
%   pMin     the smallest p each pair could reach; 0 when it has no floor
%   group    the correction family of each pair
%   correct  function handle, padj = correct(p), the correction in use
%   alpha    significance level
%
% Outputs:
%   ok   (n x 1) logical. True when the pair comes out significant with all
%        pairs of the same floor at that floor. This is its best case: corrected
%        p-values never fall when a p-value rises.
%   msg  summary for a warning; empty when every pair is testable
%
% Cedric Cannard, iEEGLAB, 2026

ok = true(size(p));
[floors, ~, level] = unique(pMin);
for f = find(floors(:)' > 0)
    at = level == f;
    best = p;
    best(at) = floors(f);
    for g = unique(group(at))'
        in = group == g;
        padj = correct(best(in));
        ok(in & at) = padj(at(in)) < alpha;
    end
end
msg = '';
if ~all(ok)
    msg = sprintf(['%d of %d pairs cannot reach significance whatever the data: their smallest ' ...
        'possible p (%s) stays above what alpha = %g accepts after correction. They are marked ' ...
        'testable = false and count as not tested, not as "no response". More trials, or more ' ...
        'permutations (n_perm) when the smallest p is 1/(n_perm + 1), would make them testable.'], ...
        nnz(~ok), numel(ok), strjoin(compose('%.3g', unique(pMin(~ok))'), ', '), alpha);
end
end
