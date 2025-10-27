% eegplugin_ieeglab() - iEEGLAB plugin 
% 
% Analyzing intracranial EEG (iEEG) data with EEGLAB. 
%
% Copyright (C) - EEGLAB, Shwartz Center, UCSD, 2025
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function vers = eegplugin_ieeglab(fig, try_strings, catch_strings)

% Plugin version
vers = '1.0';

% Add paths to subfolders
plugin_path = fileparts(which('eegplugin_ieeglab.m'));
addpath(genpath(plugin_path));

% --- define callbacks
cb_load = [try_strings.no_check 'EEG = ieeglab_load(EEG);' catch_strings.new_and_hist];
cb_vis_elec = [try_strings.no_check 'EEG = ieeglab_vis_elec(EEG);' catch_strings.new_and_hist];
cb_preprocess = [try_strings.no_check 'EEG = ieeglab_preprocess(EEG);' catch_strings.new_and_hist];
cb_stats   = [try_strings.no_check '[EEG, LASTCOM] = ieeglab_stats_subject(EEG);'   catch_strings.new_and_hist];

% --- remove any existing copy to avoid duplicates on rehash
old = findobj(fig, 'Type', 'uimenu', 'Tag', 'menu_ieeglab');
if ~isempty(old), delete(old); end

% --- create a NEW TOP-LEVEL MENU on the EEGLAB menubar
menu_root = uimenu(fig, ...
    'Label',     'iEEGLAB', ...
    'Tag',       'menu_ieeglab', ...
    'Separator', 'on', ...
    'Position',  7);  % adjust position to place it where you like

% --- add items
uimenu(menu_root, 'Label', 'Load electrode coordinates and events', 'Callback', cb_load);
uimenu(menu_root, 'Label', 'Visualize electrodes', 'Callback', cb_vis_elec);
uimenu(menu_root, 'Label', 'Preprocess iEEG data', 'Callback', cb_preprocess);
uimenu(menu_root, 'Label', 'Within-subject statistics','Callback', cb_stats);

end
