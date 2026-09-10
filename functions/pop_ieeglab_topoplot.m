function com = pop_ieeglab_topoplot(EEG, what, varargin)
% pop_ieeglab_topoplot() - Dialog front end for ieeglab_topoplot: plot one value
%                          per electrode on the brain.
%
% Usage:
%   pop_ieeglab_topoplot(EEG)                      % dialog
%   pop_ieeglab_topoplot(EEG, 120)                 % same as ieeglab_topoplot(EEG, 120)
%   pop_ieeglab_topoplot(EEG, 'in_degree')
%
% This is the entry point EEGLAB's Plot menu uses for intracranial data, in
% place of a 2D scalp topography (see ieeglab_topoplot for why).
%
% Cedric Cannard, iEEGLAB, 2026

com = '';
if nargin >= 2
    ieeglab_topoplot(EEG, what, varargin{:});
    com = sprintf('pop_ieeglab_topoplot(EEG, %s);', local_lit(what));
    return
end

hasM = isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'ccep_matrix') && ~isempty(EEG.ieeglab.ccep_matrix);
kinds = {'Amplitude at a latency (ms)', 'Mean amplitude over a window (ms)'};
if hasM
    kinds = [kinds, {'In-degree (connectivity)', 'N1 amplitude, for one stimulation site', ...
                     'N1 latency, for one stimulation site'}];
    sites = EEG.ieeglab.ccep_matrix.sites(:)';
else
    sites = {'(run CCEP analysis to enable)'};
end
defLat = '20';
if isfield(EEG,'times') && ~isempty(EEG.times), defLat = num2str(round(min(max(20, EEG.times(1)), EEG.times(end)))); end

uilist = {
    {'style' 'text' 'string' 'What to show'}  {'style' 'popupmenu' 'string' kinds 'tag' 'kind'}
    {'style' 'text' 'string' 'Latency, or "start stop" window (ms)'} {'style' 'edit' 'string' defLat 'tag' 'lat'}
    {'style' 'text' 'string' 'Stimulation site'} {'style' 'popupmenu' 'string' sites 'tag' 'site'}
};
uilist = [uilist(:,1) uilist(:,2)]'; uilist = uilist(:)';
[res, ~, ~, out] = inputgui({[1 1] [1 1] [1 1]}, uilist, 'pophelp(''ieeglab_topoplot'')', ...
    'iEEGLAB - electrode values on brain');
if isempty(res), return; end

switch out.kind
    case {1, 2}
        w = str2num(out.lat); %#ok<ST2NM>
        if isempty(w) || (out.kind == 2 && numel(w) ~= 2)
            error('pop_ieeglab_topoplot:badLatency', 'Enter a latency in ms, or "start stop" for a window.');
        end
        if out.kind == 1, w = w(1); end
        ieeglab_topoplot(EEG, w);
        com = sprintf('pop_ieeglab_topoplot(EEG, %s);', mat2str(w));
    case 3
        ieeglab_topoplot(EEG, 'in_degree');
        com = 'pop_ieeglab_topoplot(EEG, ''in_degree'');';
    otherwise
        metric = 'n1_amplitude'; if out.kind == 5, metric = 'n1_latency'; end
        site = sites{out.site};
        ieeglab_topoplot(EEG, metric, 'site', site);
        com = sprintf('pop_ieeglab_topoplot(EEG, ''%s'', ''site'', ''%s'');', metric, site);
end
end

function s = local_lit(x)
if ischar(x) || isstring(x), s = ['''' char(x) '''']; else, s = mat2str(x); end
end
