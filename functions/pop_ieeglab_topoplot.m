function com = pop_ieeglab_topoplot(EEG, what, varargin)
% pop_ieeglab_topoplot() - Dialog front end for ieeglab_topoplot: plot one value
%                          per electrode on the brain.
%
% Usage:
%   pop_ieeglab_topoplot(EEG)                                  % dialog
%   pop_ieeglab_topoplot(EEG, 120)                             % as ieeglab_topoplot(EEG, 120)
%   pop_ieeglab_topoplot(EEG, 'n1_amplitude', 'site', 'RA1-RA2')
%
% The entry point EEGLAB's Plot menu uses for intracranial data, in place of a
% 2D scalp topography (see ieeglab_topoplot for why). The per-site entries are
% built from whatever the connectivity matrix contains: N1 amplitude and
% latency for an N1 matrix, CRP duration and explained variance for a CRP one.
%
% Cedric Cannard, iEEGLAB, 2026

com = '';
if nargin >= 2
    ieeglab_topoplot(EEG, what, varargin{:});
    com = local_com(what, varargin);
    return
end

kinds = {'Amplitude at a latency (ms)', 'Mean amplitude over a window (ms)'};
metrics = {'latency', 'window'};
sites = {'(none)'};
if isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'ccep_matrix') && ~isempty(EEG.ieeglab.ccep_matrix)
    M = EEG.ieeglab.ccep_matrix;
    sites = M.sites(:)';
    kinds{end+1} = 'In-degree (sites evoking a response)'; metrics{end+1} = 'in_degree';
    if isfield(M, 'amplitude_uv')
        kinds{end+1} = 'N1 amplitude, for one stimulation site'; metrics{end+1} = 'n1_amplitude';
        kinds{end+1} = 'N1 latency, for one stimulation site';   metrics{end+1} = 'n1_latency';
    end
    if isfield(M, 'tR_ms')
        kinds{end+1} = 'CRP response duration, for one stimulation site'; metrics{end+1} = 'crp_tr';
        kinds{end+1} = 'CRP explained variance, for one stimulation site'; metrics{end+1} = 'crp_explained_var';
    end
end
defLat = '20';
if isfield(EEG,'times') && ~isempty(EEG.times), defLat = num2str(round(min(max(20, EEG.times(1)), EEG.times(end)))); end

uilist = { ...
    {'style' 'text' 'string' 'What to show'}, {'style' 'popupmenu' 'string' kinds 'tag' 'kind'}, ...
    {'style' 'text' 'string' 'Latency, or "start stop" window (ms)'}, {'style' 'edit' 'string' defLat 'tag' 'lat'}, ...
    {'style' 'text' 'string' 'Stimulation site'}, {'style' 'popupmenu' 'string' sites 'tag' 'site'} };
[res, ~, ~, out] = inputgui({[1 1] [1 1] [1 1]}, uilist, 'pophelp(''ieeglab_topoplot'')', ...
    'iEEGLAB - electrode values on brain');
if isempty(res), return; end

m = metrics{out.kind};
switch m
    case {'latency', 'window'}
        w = str2num(out.lat); %#ok<ST2NM>
        if isempty(w) || (strcmp(m,'window') && numel(w) ~= 2)
            error('pop_ieeglab_topoplot:badLatency', 'Enter a latency in ms, or "start stop" for a window.');
        end
        if strcmp(m, 'latency'), w = w(1); end
        args = {};
        if ~strcmp(sites{1}, '(none)') && out.site >= 1, args = {'site', sites{out.site}}; end
        ieeglab_topoplot(EEG, w, args{:});
        com = local_com(w, args);
    case 'in_degree'
        ieeglab_topoplot(EEG, 'in_degree');
        com = local_com('in_degree', {});
    otherwise
        args = {'site', sites{out.site}};
        ieeglab_topoplot(EEG, m, args{:});
        com = local_com(m, args);
end
end

function com = local_com(what, args)
% History line that re-runs the same call, name/value pairs included (quotes
% doubled, so sEEG labels such as LA'3 do not break it). 'draw' is dropped.
parts = {local_lit(what)};
for k = 1:2:numel(args)-1
    if strcmpi(args{k}, 'draw'), continue; end
    parts{end+1} = local_lit(args{k});   %#ok<AGROW>
    parts{end+1} = local_lit(args{k+1}); %#ok<AGROW>
end
com = sprintf('pop_ieeglab_topoplot(EEG, %s);', strjoin(parts, ', '));
end

function s = local_lit(x)
if ischar(x) || (isstring(x) && isscalar(x))
    s = ['''' strrep(char(x), '''', '''''') ''''];
elseif iscell(x) || isstring(x)
    q = cellfun(@(v) ['''' strrep(char(v), '''', '''''') ''''], cellstr(x), 'UniformOutput', false);
    s = ['{' strjoin(q, ',') '}'];
else
    s = mat2str(x);
end
end
