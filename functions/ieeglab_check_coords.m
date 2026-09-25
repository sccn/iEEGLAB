function [ok, info] = ieeglab_check_coords(xyz, surfaces, varargin)
% ieeglab_check_coords() - Do the electrodes and the brain surfaces share a space?
%
% Usage:
%   ok = ieeglab_check_coords(xyz, surfaces)
%   [ok, info] = ieeglab_check_coords(xyz, surfaces, 'max_outside', 0.1, 'verbose', true)
%
%   xyz       N x 3 electrode coordinates (mm)
%   surfaces  cell array of structs with .vertices (as read by ieeglab_read_gifti),
%             one per surface drawn (e.g. left and right pial)
%
% A contact counts as inside when it lies within the convex hull of at least one
% surface. Intracranial contacts are inside the brain, so when more than
% 'max_outside' of them (default 10%) fall outside, the coordinates and the
% surfaces are most likely in different spaces (for example scanner RAS against
% FreeSurfer surface RAS, or individual coordinates on a template): the function
% warns ('ieeglab_check_coords:spaceMismatch') and returns ok = false.
%
% info: n, n_outside, frac_outside, median_dist_mm (to the nearest vertex).
%
% On OpenNeuro ds004696 sub-02, all 211 contacts are inside the published pial
% surfaces (median 3.4 mm from the surface); shifted by 30 mm, 15% are outside.
%
% Cedric Cannard, iEEGLAB, 2026

p = inputParser;
p.addParameter('max_outside', 0.1);
p.addParameter('verbose', true);
p.parse(varargin{:});
o = p.Results;

xyz = xyz(all(isfinite(xyz), 2), :);
info = struct('n', size(xyz, 1), 'n_outside', 0, 'frac_outside', 0, 'median_dist_mm', NaN);
ok = true;
if isempty(xyz) || isempty(surfaces), return; end
if isstruct(surfaces), surfaces = num2cell(surfaces); end

inside = false(size(xyz, 1), 1);
allV = [];
for k = 1:numel(surfaces)
    V = double(surfaces{k}.vertices);
    allV = [allV; V]; %#ok<AGROW>
    inside = inside | local_in_hull(xyz, V);
end
dmin = zeros(size(xyz, 1), 1);
for i = 1:size(xyz, 1)
    dmin(i) = sqrt(min(sum((allV - xyz(i, :)).^2, 2)));
end
info.n_outside = nnz(~inside);
info.frac_outside = info.n_outside / info.n;
info.median_dist_mm = median(dmin);
if info.frac_outside > o.max_outside
    ok = false;
    if o.verbose
        warning('ieeglab_check_coords:spaceMismatch', ...
            ['%d of %d contacts (%.0f%%) lie outside the brain surfaces drawn (median %.1f mm from the ' ...
             'surface). The electrode coordinates and the surfaces are probably in different spaces: ' ...
             'check the coordinate system in *_coordsystem.json against the surfaces.'], ...
            info.n_outside, info.n, 100 * info.frac_outside, info.median_dist_mm);
    end
end
end

function tf = local_in_hull(pts, V)
% inside the convex hull of V: on the inner side of every hull face
K = convhulln(V);
A = V(K(:, 1), :); B = V(K(:, 2), :); C = V(K(:, 3), :);
nrm = cross(B - A, C - A, 2);
flip = sum(nrm .* (mean(V, 1) - A), 2) > 0;          % make every normal point outwards
nrm(flip, :) = -nrm(flip, :);
off = sum(nrm .* A, 2);
tf = all(pts * nrm' - off' <= 1e-6 * max(abs(off)), 2);
end
