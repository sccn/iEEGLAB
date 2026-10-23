function [az, el] = ieeglab_brain_view(EEG, xyz)
% ieeglab_brain_view() - The direction a brain figure opens from.
%
% Usage:
%   [az, el] = ieeglab_brain_view(EEG, xyz)     % then ieeg_viewLight(az, el)
%
% Depth electrodes (sEEG, DBS): from the front and slightly above, so that both
% hemispheres and leads at any angle show. ECoG grids and strips lie on the
% cortex: seen face-on, along the normal of the plane that best fits the contacts,
% from outside the brain (a grid on the front-left of the frontal lobe, as in
% ds004080, is seen from the front-left; from the front or the side it looks
% jumbled). ECoG when ECOG channels outnumber SEEG and DBS ones (BIDS channel
% types). xyz: N x 3 contact coordinates (RAS, mm), NaN allowed.
%
% Cedric Cannard, iEEGLAB, 2026

az = 180; el = 20;                                        % front, slightly above
ty = {};
if isfield(EEG, 'chanlocs') && isfield(EEG.chanlocs, 'type')
    ty = upper(cellfun(@(x) char(string(x)), {EEG.chanlocs.type}, 'UniformOutput', false));
end
if nnz(strcmp(ty, 'ECOG')) <= nnz(strcmp(ty, 'SEEG')) + nnz(strcmp(ty, 'DBS')), return; end
xyz = xyz(all(isfinite(xyz), 2), :);
el = 0; az = 90;                                          % right lateral
if ~isempty(xyz) && median(xyz(:, 1)) < 0, az = -90; end % left hemisphere: left lateral
if size(xyz, 1) < 3, return; end
c = mean(xyz, 1);
[V, D] = eig(cov(xyz));
[~, k] = min(diag(D)); n = V(:, k)';                      % normal of the best-fitting plane
centre = [0 -20 10];                                      % about the middle of the brain (MNI/fsaverage, mm)
if dot(n, c - centre) < 0, n = -n; end                    % from outside the brain
el = asind(max(-1, min(1, n(3))));
az = atan2d(n(1), -n(2));                                 % MATLAB view: (0,0) from the back, (90,0) from the right
end
