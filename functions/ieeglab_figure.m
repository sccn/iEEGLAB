function f = ieeglab_figure(varargin)
% ieeglab_figure() - New iEEGLAB figure: its own window, white background, no
%                    menu bar, toolbar or figure toolstrip.
%
% Usage:
%   f = ieeglab_figure('Name', 'CRP summary', 'Position', [100 100 900 360])
%
% Any figure property can be passed and overrides these defaults.
%
% Cedric Cannard, iEEGLAB, 2026

f = figure('Color', 'w', 'NumberTitle', 'off', 'WindowStyle', 'normal', ...
    'MenuBar', 'none', 'ToolBar', 'none', varargin{:});
end
