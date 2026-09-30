function sulc = ieeglab_surface_sulc(g, nIter)
% ieeglab_surface_sulc() - Sulcal depth estimated from a cortical mesh alone.
%
% Usage:
%   sulc = ieeglab_surface_sulc(g)
%   sulc = ieeglab_surface_sulc(g, nIter)
%
% Inputs:
%   g      - struct with .vertices (N x 3) and .faces (M x 3, 1-based)
%   nIter  - smoothing iterations for the reference surface (default 60)
%
% Output:
%   sulc   - N x 1, positive in sulci and negative on gyri, in mm, with the
%            same sign convention as FreeSurfer's ?h.sulc. Pass it to
%            ieeg_RenderGifti(g, sulc) to draw sulci darker than gyri.
%
% Each vertex is compared with a heavily smoothed copy of the surface: sulcal
% vertices lie below it (along the outward normal), gyral crowns above it.
% Drawing this as vertex colour shows the folding pattern without relying on
% lighting, which a transparent surface barely shows in MATLAB's current
% graphics.
%
% Cedric Cannard, iEEGLAB, 2026

if nargin < 2 || isempty(nIter), nIter = 60; end
V = double(g.vertices); F = double(g.faces); nv = size(V, 1);

A = sparse(F(:, [1 2 3]), F(:, [2 3 1]), 1, nv, nv);
A = double((A + A') > 0);
W = spdiags(1 ./ max(sum(A, 2), 1), 0, nv, nv) * A;   % neighbour average

Vs = V;
for k = 1:nIter, Vs = W * Vs; end

% outward vertex normals of the smoothed surface
Nf = cross(Vs(F(:,2),:) - Vs(F(:,1),:), Vs(F(:,3),:) - Vs(F(:,1),:), 2);
N = zeros(nv, 3);
for c = 1:3
    N(:, c) = accumarray(F(:), repmat(Nf(:, c), 3, 1), [nv 1]);
end
N = N ./ max(vecnorm(N, 2, 2), eps);
if mean(sum(N .* (Vs - mean(Vs, 1)), 2)) < 0, N = -N; end

sulc = sum((Vs - V) .* N, 2);
end
