function [ok, report] = ieeglab_check_install(mode)
% ieeglab_check_install() - Verify that iEEGLAB's dependencies are installed.
%
% Usage:
%   ieeglab_check_install                 % print a report
%   ok = ieeglab_check_install('quiet')   % just the boolean
%   ieeglab_check_install('assert')       % raise an error if anything required
%                                         % is missing (used by the menu callbacks)
%
% Checks the requirements in three tiers:
%   required - the plugin cannot work at all without these
%   feature  - a specific feature is unavailable, the rest still works
%   optional - only needed for particular data formats
%
% Output:
%   ok     - true when every REQUIRED dependency resolves
%   report - struct array with name, kind, found, resolved_path, needed_for, hint
%
% Cedric Cannard, iEEGLAB, 2026

if nargin < 1, mode = 'print'; end

D = { ...
  % name              kind        needed for                                  hint
  'eeglab',          'required', 'everything',                                'https://github.com/sccn/eeglab';
  'eeg_checkset',    'required', 'dataset handling',                          'Part of EEGLAB; run eeglab once to set the path.';
  'pop_select',      'required', 'channel and trial selection',               'Part of EEGLAB.';
  'pop_epoch',       'required', 'epoching',                                  'Part of EEGLAB.';
  'pop_resample',    'required', 'downsampling',                              'Part of EEGLAB.';
  'inputgui',        'required', 'all dialogs',                               'Part of EEGLAB (functions/guifunc).';
  'readtable',       'required', 'reading BIDS .tsv electrode and event files','Base MATLAB.';
  'pop_eegfiltnew',  'feature',  'filtering',                                 'Enable the firfilt plugin from the EEGLAB plugin manager.';
  'gifti',           'feature',  '3D glass-brain visualisation of surfaces',  'Install vistasoft: https://github.com/vistalab/vistasoft then addpath(genpath(vistasoft_dir)).';
  'trimmean',        'feature',  'trimmed-mean averaging in plots and CAR',   'Statistics and Machine Learning Toolbox.';
  'designfilt',      'feature',  'CARLA''s line-noise notch before ranking',  'Signal Processing Toolbox. CARLA still runs without it, slightly less robustly.';
  'dipfitdefs',      'feature',  'template-brain fallback when no surfaces',  'Enable the dipfit plugin from the EEGLAB plugin manager.';
  'niftiread',       'optional', 'reading MRI volumes',                       'Image Processing Toolbox.';
  'read_mef_header', 'optional', 'importing MEF3 (.mefd) recordings',         'Install matmef: https://github.com/MaxvandenBoom/matmef';
};

report = struct('name',{},'kind',{},'found',{},'resolved_path',{},'needed_for',{},'hint',{});
for i = 1:size(D,1)
    p = which(D{i,1});
    if isempty(p) && exist(D{i,1},'file') ~= 0, p = '(built-in)'; end
    report(i) = struct('name',D{i,1}, 'kind',D{i,2}, 'found',~isempty(p), ...
        'resolved_path',p, 'needed_for',D{i,3}, 'hint',D{i,4});
end

% The four menu callbacks must resolve, or the menu is broken
menuFcns = {'ieeglab_load','ieeglab_vis_elec','ieeglab_preprocess','ieeglab_stats_subject'};
for i = 1:numel(menuFcns)
    p = which(menuFcns{i});
    report(end+1) = struct('name',menuFcns{i}, 'kind','required', 'found',~isempty(p), ...
        'resolved_path',p, 'needed_for','iEEGLAB menu item', ...
        'hint','Part of iEEGLAB; re-run addpath(genpath(plugin_dir)).'); %#ok<AGROW>
end

isReq     = strcmp({report.kind}, 'required');
missReq   = isReq & ~[report.found];
missOther = ~isReq & ~[report.found];
ok = ~any(missReq);

switch lower(mode)
    case 'quiet'
        return

    case 'assert'
        if ~ok
            names = {report(missReq).name};
            error('ieeglab:missingDependency', ...
                ['iEEGLAB cannot run: %d required dependency/dependencies missing (%s).\n' ...
                 'Run ieeglab_check_install for details.'], ...
                numel(names), strjoin(names, ', '));
        end
        return

    otherwise
        fprintf('\n===== iEEGLAB installation check =====\n');
        fprintf('MATLAB %s on %s\n\n', version, computer);
        for kind = {'required','feature','optional'}
            sel = find(strcmp({report.kind}, kind{1}));
            fprintf('--- %s ---\n', upper(kind{1}));
            for i = sel
                mark = '  ok  ';
                if ~report(i).found, mark = ' MISS '; end
                fprintf('[%s] %-18s %s\n', mark, report(i).name, ...
                    iff(report(i).found, report(i).resolved_path, ['needed for ' report(i).needed_for]));
                if ~report(i).found
                    fprintf('                          -> %s\n', report(i).hint);
                end
            end
            fprintf('\n');
        end
        if ok && ~any(missOther)
            fprintf('All dependencies satisfied.\n\n');
        elseif ok
            fprintf(['Ready to use. %d non-essential dependency/dependencies missing - the\n' ...
                     'features listed above will be unavailable.\n\n'], nnz(missOther));
        else
            fprintf(2, 'NOT ready: %d required dependency/dependencies missing.\n\n', nnz(missReq));
        end
end

end

function s = iff(c, a, b)
if c, s = a; else, s = b; end
end
