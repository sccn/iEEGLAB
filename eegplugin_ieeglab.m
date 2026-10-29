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
vers = '1.2.0';

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
% Import returns nothing when cancelled: EEG is replaced only when a recording was read.
cb_import     = [try_strings.no_check '[EEGTMP_, LASTCOM] = pop_ieeglab_import; if ~isempty(LASTCOM), EEG = EEGTMP_; end; clear EEGTMP_;' catch_strings.new_and_hist];
cb_load       = [try_strings.no_check '[EEG, LASTCOM] = ieeglab_load(EEG);'           catch_strings.new_and_hist];
cb_vis_elec   = [try_strings.no_check 'EEG = ieeglab_vis_elec(EEG); LASTCOM = ''EEG = ieeglab_vis_elec(EEG);'';' catch_strings.new_and_hist];
cb_surfaces   = [try_strings.no_check '[EEG, LASTCOM] = pop_ieeglab_surfaces(EEG);'    catch_strings.new_and_hist];
cb_preprocess = [try_strings.no_check '[EEG, LASTCOM] = ieeglab_preprocess(EEG);'     catch_strings.new_and_hist];
cb_reref      = [try_strings.no_check '[EEG, LASTCOM] = pop_ieeglab_reref(EEG);'      catch_strings.new_and_hist];
cb_spectral   = [try_strings.no_check '[~, LASTCOM] = pop_ieeglab_spectral(EEG);'    catch_strings.add_to_hist];   % figures only: the dataset is unchanged
cb_stats      = [try_strings.no_check '[EEG, LASTCOM] = ieeglab_stats_subject(EEG);'  catch_strings.new_and_hist];
cb_matrix     = [try_strings.no_check 'LASTCOM = pop_ieeglab_ccep_matrix(EEG);' catch_strings.add_to_hist];
cb_crp        = [try_strings.no_check 'LASTCOM = ieeglab_plot_crp(EEG);' catch_strings.add_to_hist];
cb_topo       = [try_strings.no_check 'LASTCOM = pop_ieeglab_topoplot(EEG);' catch_strings.add_to_hist];
cb_heatmap    = [try_strings.no_check 'LASTCOM = ieeglab_plot_heatmap(EEG);' catch_strings.add_to_hist];
cb_export     = [try_strings.no_check '[~, LASTCOM] = ieeglab_export(EEG);' catch_strings.add_to_hist];
cb_check      = 'ieeglab_check_install;';

% --- remove any existing copy to avoid duplicates on rehash
old = findobj(fig, 'Type', 'uimenu', 'Tag', 'menu_ieeglab');
if ~isempty(old), delete(old); end

% --- create a NEW TOP-LEVEL MENU on the EEGLAB menubar, right of Help
top = findobj(allchild(fig), 'flat', 'Type', 'uimenu');
labels = get(top, 'Label'); if ~iscell(labels), labels = {labels}; end
iHelp = find(strcmpi(strtrim(strrep(labels, '&', '')), 'Help'), 1);
pos = numel(top) + 1;
if ~isempty(iHelp), pos = get(top(iHelp), 'Position') + 1; end
menu_root = uimenu(fig, ...
    'Label',     'iEEGLAB', ...
    'Tag',       'menu_ieeglab', ...
    'Position',  pos, ...
    'userdata',  'study:on');

% EEGLAB greys out items tagged 'startup:off' until a dataset is loaded, items
% without 'study:on' while a STUDY is selected, and 'continuous:off' items on
% continuous data. Check installation stays available in every state.
ondata  = 'startup:off';
onepoch = 'startup:off;continuous:off';
uimenu(menu_root, 'Label', 'Check installation', 'Callback', cb_check, 'userdata', 'study:on');
uimenu(menu_root, 'Label', 'Import raw recording (MEF3, BrainVision, EDF)', 'Callback', cb_import, 'Separator', 'on', 'userdata', 'study:on');
uimenu(menu_root, 'Label', 'Load electrode coordinates and events', 'Callback', cb_load, 'userdata', ondata);
% The brain is chosen once here; every brain figure then uses it without asking
uimenu(menu_root, 'Label', 'Choose brain surfaces', 'Callback', cb_surfaces, 'userdata', ondata);
uimenu(menu_root, 'Label', 'Visualize electrodes', 'Callback', cb_vis_elec, 'userdata', ondata);
uimenu(menu_root, 'Label', 'Preprocess iEEG data', 'Callback', cb_preprocess, 'userdata', ondata);
uimenu(menu_root, 'Label', 'iEEG re-referencing (common average, CARLA, ICA, bipolar, montage)', 'Callback', cb_reref, 'userdata', ondata);
uimenu(menu_root, 'Label', 'Frequency analysis (spectra, time-frequency, band power)', 'Callback', cb_spectral, 'userdata', ondata);
% Epoched data only (greyed out until Preprocess iEEG data has segmented the data)
uimenu(menu_root, 'Label', 'Plot average of epochs, every channel (heatmap)', 'Callback', cb_heatmap, 'Separator', 'on', 'userdata', onepoch);
uimenu(menu_root, 'Label', 'CCEP analysis (N1, CRP, connectivity)', 'Callback', cb_stats, 'userdata', onepoch);
uimenu(menu_root, 'Label', 'Plot CRP of one site and contact', 'Callback', cb_crp, 'userdata', onepoch);
uimenu(menu_root, 'Label', 'Plot connectivity matrix', 'Callback', cb_matrix, 'userdata', onepoch);
uimenu(menu_root, 'Label', 'Plot electrode values on brain', 'Callback', cb_topo, 'userdata', onepoch);
uimenu(menu_root, 'Label', 'Export results (TSV / JSON / MAT)', 'Callback', cb_export, 'Separator', 'on', 'userdata', onepoch);
local_docs_menu(menu_root, plugin_path, true);

% The same documentation in EEGLAB's Help menu, next to EEGLAB's own
if ~isempty(iHelp)
    old = findobj(top(iHelp), 'Tag', 'ieeglab_help');
    if ~isempty(old), delete(old); end
    local_docs_menu(top(iHelp), plugin_path, true);
end

% Also offer the electrode-value plot under EEGLAB's own Plot menu, next to the
% scalp topographies it replaces for intracranial data.
try
    plotMenu = findobj(fig, 'Type', 'uimenu', 'Tag', 'plot');
    if ~isempty(plotMenu)
        old = findobj(plotMenu(1), 'Tag', 'ieeglab_topo');
        if ~isempty(old), delete(old); end
        uimenu(plotMenu(1), 'Label', 'iEEG electrode values on brain (iEEGLAB)', ...
            'Tag', 'ieeglab_topo', 'Callback', cb_topo, 'Separator', 'on', 'userdata', onepoch);
    end
catch
end

% Fail loudly at load time if a menu callback points at a function that does not
% exist, rather than at click time. ieeglab_stats_subject was missing for months
% because nothing checked this.
for f = {'pop_ieeglab_import','ieeglab_load','pop_ieeglab_surfaces','ieeglab_vis_elec','ieeglab_preprocess','pop_ieeglab_reref','pop_ieeglab_spectral', ...
         'ieeglab_stats_subject','pop_ieeglab_ccep_matrix','pop_ieeglab_topoplot','ieeglab_export','ieeglab_plot_heatmap','ieeglab_plot_crp'}
    if isempty(which(f{1}))
        warning('eegplugin_ieeglab:missingCallback', ...
            'iEEGLAB menu item calls %s, which is not on the path. That menu item will fail.', f{1});
    end
end

end

function local_docs_menu(parent, plugin_path, separator)
% 'iEEGLAB documentation': the README and tutorial online, the menu-by-menu
% walkthrough shipped with the plugin, the tutorial script, and the help of
% every menu's function (the text each dialog's Help button shows).
m = uimenu(parent, 'Label', 'iEEGLAB documentation', 'Tag', 'ieeglab_help', ...
    'Separator', iff(separator, 'on', 'off'), 'userdata', 'study:on');
uimenu(m, 'Label', 'README and tutorial (web)', 'userdata', 'study:on', ...
    'Callback', 'web(''https://github.com/sccn/iEEGLAB#readme'', ''-browser'');');
walk = fullfile(plugin_path, 'iEEGLAB_GUI_walkthrough.html');
uimenu(m, 'Label', 'Menu-by-menu walkthrough', 'userdata', 'study:on', ...
    'Callback', sprintf('web(''%s'', ''-browser'');', strrep(walk, '''', '''''')));
uimenu(m, 'Label', 'Tutorial script (ieeglab_tutorial.m)', 'userdata', 'study:on', 'Callback', 'edit ieeglab_tutorial;');
fm = uimenu(m, 'Label', 'iEEGLAB functions', 'Separator', 'on', 'userdata', 'study:on');
fns = {'pop_ieeglab_import', 'Import raw recording'; 'ieeglab_load', 'Load electrode coordinates and events'; ...
       'pop_ieeglab_surfaces', 'Choose brain surfaces'; 'ieeglab_vis_elec', 'Visualize electrodes'; ...
       'ieeglab_preprocess', 'Preprocess iEEG data'; ...
       'pop_ieeglab_reref', 'iEEG re-referencing'; 'pop_ieeglab_spectral', 'Frequency analysis'; ...
       'ieeglab_stats_subject', 'CCEP analysis'; 'pop_ieeglab_ccep_matrix', 'Connectivity matrix'; ...
       'pop_ieeglab_topoplot', 'Electrode values on brain'; 'ieeglab_export', 'Export results'};
for k = 1:size(fns, 1)
    uimenu(fm, 'Label', sprintf('%s (%s)', fns{k, 2}, fns{k, 1}), 'userdata', 'study:on', ...
        'Callback', sprintf('pophelp(''%s'');', fns{k, 1}));
end
end

function out = iff(c, a, b)
if c, out = a; else, out = b; end
end
