% eegplugin_ieeglab() - iEEGLAB plugin 
% 
% Analyzing intracranial EEG (iEEG) data with EEGLAB. 
%
% Copyright (C) - EEGLAB, Swartz Center for Computational Neuroscience, UCSD, 2025-2026
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, see <https://www.gnu.org/licenses/>.

function vers = eegplugin_ieeglab(fig, try_strings, catch_strings)

% Plugin version
vers = '1.1.0';

if nargin < 3
    error('eegplugin_ieeglab:badCall', ...
        ['eegplugin_ieeglab requires 3 arguments and is called by eeglab(), not directly.\n' ...
         'To check your installation, run: ieeglab_check_install']);
end

% Add the plugin folder and its functions to the path.
plugin_path = fileparts(which('eegplugin_ieeglab.m'));
addpath(plugin_path);
addpath(fullfile(plugin_path, 'functions'));

% --- define callbacks
% Each callback assigns LASTCOM so that catch_strings.new_and_hist can store the
% modified dataset and refresh the main EEGLAB window (issue #4). Without it
% eeglab_new took neither branch and ALLEEG kept the un-preprocessed dataset.
cb_load       = [try_strings.no_check '[EEG, LASTCOM] = ieeglab_load(EEG);'           catch_strings.new_and_hist];
cb_vis_elec   = [try_strings.no_check 'EEG = ieeglab_vis_elec(EEG); LASTCOM = ''EEG = ieeglab_vis_elec(EEG);'';' catch_strings.new_and_hist];
cb_preprocess = [try_strings.no_check '[EEG, LASTCOM] = ieeglab_preprocess(EEG);'     catch_strings.new_and_hist];
cb_reref      = [try_strings.no_check '[EEG, LASTCOM] = pop_ieeglab_reref(EEG);'      catch_strings.new_and_hist];
cb_stats      = [try_strings.no_check '[EEG, LASTCOM] = ieeglab_stats_subject(EEG);'  catch_strings.new_and_hist];
cb_matrix     = [try_strings.no_check '[ok_, why_] = ieeglab_matrix_current(EEG); if ok_, ieeglab_plot_ccep_matrix(EEG.ieeglab.ccep_matrix); else, errordlg([''Cannot plot the connectivity matrix: '' why_ ''. Run iEEGLAB > CCEP analysis.''],''iEEGLAB''); end; clear ok_ why_; LASTCOM = '''';' catch_strings.add_to_hist];
cb_topo       = [try_strings.no_check 'LASTCOM = pop_ieeglab_topoplot(EEG);' catch_strings.add_to_hist];
cb_export     = [try_strings.no_check '[~, LASTCOM] = ieeglab_export(EEG);' catch_strings.add_to_hist];
cb_check      = 'ieeglab_check_install;';

% --- remove any existing copy to avoid duplicates on rehash
old = findobj(fig, 'Type', 'uimenu', 'Tag', 'menu_ieeglab');
if ~isempty(old), delete(old); end

% --- create a NEW TOP-LEVEL MENU on the EEGLAB menubar
menu_root = uimenu(fig, ...
    'Label',     'iEEGLAB', ...
    'Tag',       'menu_ieeglab', ...
    'Separator', 'on', ...
    'Position',  7);

uimenu(menu_root, 'Label', 'Load electrode coordinates and events', 'Callback', cb_load);
uimenu(menu_root, 'Label', 'Visualize electrodes', 'Callback', cb_vis_elec);
uimenu(menu_root, 'Label', 'Preprocess iEEG data', 'Callback', cb_preprocess);
uimenu(menu_root, 'Label', 'iEEG re-referencing', 'Callback', cb_reref);
uimenu(menu_root, 'Label', 'CCEP analysis (N1, CRP, connectivity)', 'Callback', cb_stats, 'Separator', 'on');
uimenu(menu_root, 'Label', 'Plot connectivity matrix', 'Callback', cb_matrix);
uimenu(menu_root, 'Label', 'Plot electrode values on brain', 'Callback', cb_topo);
uimenu(menu_root, 'Label', 'Export results (TSV / JSON / MAT)', 'Callback', cb_export, 'Separator', 'on');
uimenu(menu_root, 'Label', 'Check installation', 'Callback', cb_check, 'Separator', 'on');

% Also offer the electrode-value plot under EEGLAB's own Plot menu, next to the
% scalp topographies it replaces for intracranial data.
try
    plotMenu = findobj(fig, 'Type', 'uimenu', 'Tag', 'plot');
    if ~isempty(plotMenu)
        old = findobj(plotMenu(1), 'Tag', 'ieeglab_topo');
        if ~isempty(old), delete(old); end
        uimenu(plotMenu(1), 'Label', 'iEEG electrode values on brain (iEEGLAB)', ...
            'Tag', 'ieeglab_topo', 'Callback', cb_topo, 'Separator', 'on');
    end
catch
end

% Fail loudly at load time if a menu callback points at a function that does not
% exist, rather than at click time. ieeglab_stats_subject was missing for months
% because nothing checked this.
for f = {'ieeglab_load','ieeglab_vis_elec','ieeglab_preprocess','pop_ieeglab_reref','ieeglab_stats_subject', ...
         'ieeglab_plot_ccep_matrix','pop_ieeglab_topoplot','ieeglab_export'}
    if isempty(which(f{1}))
        warning('eegplugin_ieeglab:missingCallback', ...
            'iEEGLAB menu item calls %s, which is not on the path. That menu item will fail.', f{1});
    end
end

end
