function [opt, wasCanceled] = ieeglab_gui_load1(filepath, defaults)
% ieeglab_gui_load1() - Dialog to choose the BIDS sidecar files to load.
%
% Usage:
%   [opt, wasCanceled] = ieeglab_gui_load1(filepath)
%   [opt, wasCanceled] = ieeglab_gui_load1(filepath, defaults)
%
% Inputs:
%   filepath  - folder the file pickers open in
%   defaults  - struct with elec_tsv / events_tsv / channels_tsv used to
%               pre-fill the fields (ieeglab_load passes the BIDS siblings it
%               found next to the dataset, so the usual case is just "OK")
%
% Outputs:
%   opt       - struct with elec_tsv, events_tsv, channels_tsv ('' when not
%               chosen); [] when cancelled
%   wasCanceled - true when the dialog was cancelled or closed
%
% Cedric Cannard, iEEGLAB, 2025-2026

wasCanceled = false;
if nargin < 2 || ~isstruct(defaults), defaults = struct(); end
for f = {'elec_tsv','events_tsv','channels_tsv'}
    if ~isfield(defaults, f{1}) || isempty(defaults.(f{1})), defaults.(f{1}) = ''; end
end

if nargin >= 1 && ~isempty(filepath) && isfolder(filepath)
    start_dir = char(filepath);
else
    plugin_path = fileparts(which('eegplugin_ieeglab'));
    start_dir = pwd;
    if ~isempty(plugin_path)
        start_dir = fullfile(plugin_path, 'tutorial');
        if ~isfolder(start_dir), start_dir = plugin_path; end
    end
end

    function do_browse(src, ~, tag, what)
        [f, p] = uigetfile({'*.tsv', [what ' (*.tsv)']}, ['Select the ' what], start_dir);
        if isequal(f, 0), return; end
        set(findobj(ancestor(src,'figure'), 'tag', tag), 'string', fullfile(p, f));
        hFig = ancestor(src, 'figure'); try, figure(hFig); uistack(hFig, 'top'); drawnow; end %#ok<TRYNC>
    end

uilist = {
    {'style' 'text' 'string' 'Electrode locations (*_electrodes.tsv):' 'fontweight' 'bold'}
    {'style' 'edit' 'string' defaults.elec_tsv 'tag' 'elec_tsv' 'horizontalalignment' 'left'}
    {'style' 'pushbutton' 'string' 'Browse...' 'callback' @(s,e) do_browse(s,e,'elec_tsv','electrodes file')}
    {}
    {'style' 'text' 'string' 'Events (*_events.tsv):' 'fontweight' 'bold'}
    {'style' 'edit' 'string' defaults.events_tsv 'tag' 'events_tsv' 'horizontalalignment' 'left'}
    {'style' 'pushbutton' 'string' 'Browse...' 'callback' @(s,e) do_browse(s,e,'events_tsv','events file')}
    {}
    {'style' 'text' 'string' 'Channel status (*_channels.tsv, optional):' 'fontweight' 'bold'}
    {'style' 'edit' 'string' defaults.channels_tsv 'tag' 'channels_tsv' 'horizontalalignment' 'left'}
    {'style' 'pushbutton' 'string' 'Browse...' 'callback' @(s,e) do_browse(s,e,'channels_tsv','channels file')}
    {'style' 'text' 'string' '   Channels with status = bad in this file are marked bad (removable in Preprocess).' 'fontangle' 'italic'}
};
uigeom = { [.45 .45 .15] 1 [.45 .45 .15] 1 [.45 .45 .15] 1 };

[res, ~, ~, out] = inputgui(uigeom, uilist, 'pophelp(''ieeglab_load'')', ...
    'iEEGLAB - load BIDS sidecar files');
if isempty(res) || isempty(out)
    fprintf('iEEGLAB load dialog cancelled.\n');
    wasCanceled = true;
    opt = [];
    return
end

opt = struct();
opt.elec_tsv     = strtrim(char(out.elec_tsv));
opt.events_tsv   = strtrim(char(out.events_tsv));
opt.channels_tsv = strtrim(char(out.channels_tsv));
for f = {'elec_tsv','events_tsv','channels_tsv'}
    v = opt.(f{1});
    if ~isempty(v) && exist(v, 'file') ~= 2
        error('ieeglab_gui_load1:missingFile', 'File not found for %s: %s', strrep(f{1},'_',' '), v);
    end
end
end
