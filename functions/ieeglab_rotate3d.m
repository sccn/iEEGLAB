function ieeglab_rotate3d(fig)
% ieeglab_rotate3d() - Brain figures: rotate with the mouse, no toolbars.
%
% Usage:
%   ieeglab_rotate3d(gcf)
%
% Hides the figure and axes toolbars, turns on click-and-drag rotation, and
% keeps the lights behind the viewer while rotating so the surface never turns
% dark. Recent MATLAB versions no longer rotate a 3D plot by dragging unless
% rotation mode is switched on.
%
% Cedric Cannard, iEEGLAB, 2026

if nargin < 1, fig = gcf; end
set(fig, 'ToolBar', 'none', 'MenuBar', 'none');
for ax = findobj(fig, 'Type', 'axes')'
    axis(ax, 'vis3d');                                  % no rescaling while rotating
    try ax.Toolbar.Visible = 'off'; catch, end
end
h = rotate3d(fig);
h.ActionPostCallback = @(~, ev) local_follow(ev.Axes);
h.Enable = 'on';
end

function local_follow(ax)
for l = findobj(ax, 'Type', 'light')'
    camlight(l, 'headlight');
end
end
