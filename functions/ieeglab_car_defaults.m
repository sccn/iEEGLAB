function [timewin, neighbors] = ieeglab_car_defaults(method)
% ieeglab_car_defaults() - Default ranking window (ms) and number of neighbouring
%                          contacts left out of the reference, per re-referencing
%                          method, as in each method's publication scripts.
%
% Usage:
%   [timewin, neighbors] = ieeglab_car_defaults('carla')   % [10 300], 2
%
%   carla      [10 300] ms, 2 contacts on each side of the stimulated pair
%              (Huang et al., 2024: CARLA.m and applyCARLARealCCEPs.m)
%   varsubset  [15 500] ms, 0 (Ojeda Valencia et al., 2023: ccepPCC_loadAverageSubset)
%   otherwise  [10 300] ms, 0
%
% Cedric Cannard, iEEGLAB, 2026

switch lower(char(method))
    case 'carla'
        timewin = [10 300]; neighbors = 2;
    case 'varsubset'
        timewin = [15 500]; neighbors = 0;
    otherwise
        timewin = [10 300]; neighbors = 0;
end
end
