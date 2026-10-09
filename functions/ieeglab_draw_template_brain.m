function scale = ieeglab_draw_template_brain(xyz)
% ieeglab_draw_template_brain() - Draw a transparent template brain (dipfit
%                                 standard BEM) in the current axes.
%
% Usage:
%   scale = ieeglab_draw_template_brain(xyz)
%
% xyz are the electrode coordinates to be drawn over it (N x 3, NaN allowed).
% scale is 1000 when they look like metres (median |value| < 2) and must be
% multiplied by it to match the template (mm), else 1. The template is in MNI
% space: electrodes in a subject's native space will not line up with it.
%
% Cedric Cannard, iEEGLAB, 2026

dipfit_root = fileparts(which('dipfitdefs'));
assert(~isempty(dipfit_root), 'dipfit not found on path. Enable dipfit in EEGLAB.');
std_dir = fullfile(dipfit_root, 'standard_BEM');
f_sccn  = fullfile(std_dir, 'standard_vol_SCCN.mat');
if exist(f_sccn,'file'), S = load(f_sccn); else, S = load(fullfile(std_dir, 'standard_vol.mat')); end
assert(isfield(S,'vol'), 'standard_BEM volume not found in dipfit.');

bnd = S.vol.bnd; if numel(bnd) > 1, bnd = bnd(end); end   % innermost boundary: brain
if isfield(bnd,'pnt'), Vb = double(bnd.pnt); else, Vb = double(bnd.pos); end
if isfield(bnd,'tri'), Fb = double(bnd.tri); else, Fb = double(bnd.face); end
patch('Faces',Fb,'Vertices',Vb, 'FaceColor',[0.75 0.80 0.90],'EdgeColor','none','FaceAlpha',0.12);

v = abs(xyz(all(isfinite(xyz),2), :));
scale = 1;
if ~isempty(v) && median(v(:)) < 2, scale = 1000; end   % metres -> mm
camlight headlight; camlight right
lighting gouraud; material dull
end
