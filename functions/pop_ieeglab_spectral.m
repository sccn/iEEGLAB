function [EEG, com] = pop_ieeglab_spectral(EEG, varargin)
% pop_ieeglab_spectral() - Frequency-domain analysis: power spectrum of each
%                          contact, or time-frequency power and band power
%                          (high-frequency activity by default), with figures.
%
% Usage:
%   [EEG, com] = pop_ieeglab_spectral(EEG)                         % dialog
%   EEG = pop_ieeglab_spectral(EEG, 'type', 'spectrum')
%   EEG = pop_ieeglab_spectral(EEG, 'type', 'timefreq', 'band', [70 150], 'norm', 'relchange')
%   EEG = pop_ieeglab_spectral(EEG, 'type', 'timefreq', 'plot', false)
%
% Options: those of ieeglab_spectral, plus
%   'plot'     draw the figures (ieeglab_plot_spectral). Default true
%   'channel'  contact for the time-frequency map. Default: average of the
%              good contacts
%
% From the menu this only draws: the dataset is left as it is. From the command
% line the result is also returned in EEG.ieeglab.spectral (see
% ieeglab_spectral for its fields), to plot again or export without
% recomputing:
%   ieeglab_plot_spectral(EEG, EEG.ieeglab.spectral)
%
% Cedric Cannard, iEEGLAB, 2026

com = '';
ieeglab_require_data(EEG, 'pop_ieeglab_spectral');
epoched = EEG.trials > 1;

if nargin < 2
    args = local_dialog(EEG, epoched);
    if isempty(args), return; end
else
    args = varargin;
end
if mod(numel(args), 2), error('pop_ieeglab_spectral:args', 'Options come in name/value pairs.'); end
names = lower(cellfun(@char, args(1:2:end), 'UniformOutput', false));
vals = args(2:2:end);
doPlot = true; channel = [];
k = find(strcmp(names, 'plot')); if ~isempty(k), doPlot = logical(vals{k(end)}); end
k = find(strcmp(names, 'channel')); if ~isempty(k), channel = vals{k(end)}; end
keep = ~ismember(names, {'plot', 'channel'});
cargs = reshape([args(1:2:end); args(2:2:end)], 1, []);
cargs = cargs(repelem(keep, 2));

S = ieeglab_spectral(EEG, cargs{:});
if ~isfield(EEG, 'ieeglab') || ~isstruct(EEG.ieeglab), EEG.ieeglab = struct(); end
EEG.ieeglab.spectral = S;
if doPlot
    ieeglab_plot_spectral(EEG, S, 'channel', channel);
end

parts = cell(1, numel(args) / 2);
for i = 1:numel(parts)
    parts{i} = sprintf('''%s'', %s', names{i}, local_literal(vals{i}));
end
com = sprintf('EEG = pop_ieeglab_spectral(EEG, %s);', strjoin(parts, ', '));
end

function args = local_dialog(EEG, epoched)
args = {};
fs = EEG.srate;
types = {'timefreq', 'spectrum'};
labels = {'Time-frequency power and band power (epoched data)', 'Power spectrum of each contact (Welch)'};
if ~epoched, types = types(2); labels = labels(2); end
conds = {'(all epochs)'};
if epoched
    s = unique(ieeglab_epoch_sites(EEG));
    s = s(s ~= "");
    conds = [conds cellstr(s)];
end
norms = {'db', 'relchange', 'ratio', 'none'};
fmax = min(200, floor(fs / 2) - 1);
band = [70 150];
if band(2) >= fs / 2, band = [max(1, round(fs / 4)) fmax]; end
bandWin = [0 round(double(EEG.times(end)))];
w = [1.2 0.8];
uigeom = {[0.6 1.4], 1, w, w, w, w, w, w, 1, w, w, 1, w, w, 1};
% defaults per analysis, shown again when the analysis changes
def.timefreq = {sprintf('5:5:%d', 5 * floor(fmax / 5)), '0.2'};
def.spectrum = {sprintf('1 %d', fmax), '1'};
onType = @(h, ~) local_type_changed(h, types, def);
uilist = { ...
    {'style' 'text' 'string' 'Analysis'}, {'style' 'popupmenu' 'string' labels 'tag' 'type' 'callback' onType}, ...
    {'style' 'text' 'string' 'Frequencies and windows' 'fontweight' 'bold'}, ...
    {'style' 'text' 'string' 'Frequencies, Hz (time-frequency: list, e.g. 5:5:200; spectrum: low high)'}, ...
        {'style' 'edit' 'string' def.(types{1}){1} 'tag' 'freqs'}, ...
    {'style' 'text' 'string' 'Window, s (time-frequency: fixed Hann window; spectrum: segment)'}, ...
        {'style' 'edit' 'string' def.(types{1}){2} 'tag' 'window'}, ...
    {'style' 'text' 'string' 'Or cycles per frequency (time-frequency; empty = fixed window)'}, ...
        {'style' 'edit' 'string' '' 'tag' 'cycles'}, ...
    {'style' 'text' 'string' 'Time step, ms (time-frequency)'}, {'style' 'edit' 'string' '10' 'tag' 'timestep'}, ...
    {'style' 'text' 'string' 'Baseline, ms (whole windows inside it only)'}, ...
        {'style' 'edit' 'string' sprintf('%g 0', round(double(EEG.times(1)))) 'tag' 'baseline'}, ...
    {'style' 'text' 'string' 'Normalization to the baseline'}, ...
        {'style' 'popupmenu' 'string' {'dB (10 log10 power/baseline)', 'Relative change (power/baseline - 1)', 'Ratio', 'None (uV^2)'} 'tag' 'norm'}, ...
    {'style' 'text' 'string' 'Band' 'fontweight' 'bold'}, ...
    {'style' 'text' 'string' 'Band, Hz (default: high-frequency activity)'}, {'style' 'edit' 'string' sprintf('%g %g', band) 'tag' 'band'}, ...
    {'style' 'text' 'string' 'Window for the brain map, ms'}, {'style' 'edit' 'string' sprintf('%g %g', bandWin) 'tag' 'bandwin'}, ...
    {'style' 'text' 'string' 'Epochs' 'fontweight' 'bold'}, ...
    {'style' 'text' 'string' 'Condition (event at time 0, or stimulation site)'}, {'style' 'popupmenu' 'string' conds 'tag' 'cond'}, ...
    {'style' 'checkbox' 'string' 'Plot' 'value' 1 'tag' 'plot'}, {'style' 'text' 'string' ''}, ...
    {'style' 'text' 'string' 'CCEP data: each contact''s own stimulation epochs are left out of its average.'} };
assert(numel(uilist) == sum(cellfun(@numel, uigeom)), 'pop_ieeglab_spectral: dialog geometry and controls disagree');
[res, ~, ~, out] = inputgui('geometry', uigeom, 'uilist', uilist, 'minwidth', 640, ...
    'title', 'iEEGLAB - Frequency analysis', 'helpcom', 'pophelp(''pop_ieeglab_spectral'');');
if isempty(res), return; end

type = types{out.type};
freqs = str2num(out.freqs); %#ok<ST2NM>
if isempty(freqs), error('pop_ieeglab_spectral:freqs', 'Enter the frequencies, e.g. 5:5:200 or 1 200.'); end
if strcmp(type, 'spectrum'), freqs = [min(freqs) max(freqs)]; end
args = {'type', type, 'freqs', freqs};
win = local_number(out.window, 'Window');
if ~isempty(win), args = [args {'window', win}]; end
if strcmp(type, 'timefreq')
    cyc = local_number(out.cycles, 'Cycles');
    if ~isempty(cyc), args = [args {'cycles', cyc}]; end
    args = [args {'timestep', str2double(out.timestep), 'baseline', sscanf(out.baseline, '%f')', ...
        'norm', norms{out.norm}, 'band_window', sscanf(out.bandwin, '%f')'}];
end
bd = sscanf(out.band, '%f')';
if numel(bd) == 2, args = [args {'band', bd}]; end
if out.cond > 1, args = [args {'condition', conds{out.cond}}]; end
args = [args {'plot', logical(out.plot)}];
end

function local_type_changed(h, types, def)
% Show the frequencies and window that suit the chosen analysis
fig = ancestor(h, 'figure');
d = def.(types{get(h, 'value')});
set(findobj(fig, 'tag', 'freqs'), 'string', d{1});
set(findobj(fig, 'tag', 'window'), 'string', d{2});
end

function x = local_number(s, what)
% One positive number from a dialog box; empty box = the default
s = strtrim(s);
x = [];
if isempty(s), return; end
x = str2double(s);
if ~(isscalar(x) && isfinite(x) && x > 0)
    error('pop_ieeglab_spectral:badNumber', '%s: enter one positive number, or leave it empty for the default (got "%s").', what, s);
end
end

function s = local_literal(v)
if ischar(v) || isstring(v)
    s = ['''' strrep(char(v), '''', '''''') ''''];
elseif islogical(v) && isscalar(v)
    s = mat2str(v);
elseif isnumeric(v) && isvector(v) && numel(v) > 2 && v(2) ~= v(1) ...
        && isequal(v(1):(v(2) - v(1)):v(end), v(:)')
    s = sprintf('%s:%s:%s', mat2str(v(1)), mat2str(v(2) - v(1)), mat2str(v(end)));   % a regular list, e.g. 5:5:200
else
    s = mat2str(v);
end
end
