function [EEG, T, com] = ieeglab_reject_trials(EEG, opt)
% ieeglab_reject_trials() - Flag, and optionally remove, outlier trials.
%
% Usage:
%   EEG         = ieeglab_reject_trials(EEG)
%   [EEG, T]    = ieeglab_reject_trials(EEG, opt)
%
% A trial is flagged when its amplitude is an outlier on MANY channels at once
% compared with the other trials of the same stimulation site (or condition).
% Whole-trial artifacts - movement, a failed or double pulse, an amplifier
% saturation - look like that; a genuinely large response on a few contacts
% does not, so it is left alone.
%
% Per channel and per site, each trial's amplitude (RMS or max |x| over the
% window) is converted to a robust z-score using the median and MAD across that
% site's trials. The stimulated contacts and channels marked bad are ignored.
%
% Options:
%   .window    [start stop] ms for the amplitude measure. Default: the epoch
%              minus the stimulation window [-2 15] ms, so the artifact itself
%              is never the reason a trial is rejected.
%   .measure   'rms' (default) or 'maxabs'
%   .z         robust z threshold. Default 5.
%   .frac      a trial is flagged when at least this fraction of the usable
%              channels exceed the threshold. Default 0.25.
%   .min_trials sites with fewer trials are not assessed (no stable median). Default 5.
%   .action    'remove' (default) or 'mark'. 'mark' sets EEG.reject.rejmanual,
%              EEGLAB's own trial-rejection field, so the flags show up in
%              EEGLAB's rejection tools.
%   .verbose   default true
%
% Output T: one row per trial - trial, site, n_channels_over, frac_over,
% max_abs_z, bad.
%
% Cedric Cannard, iEEGLAB, 2026

com = '';
if nargin < 2 || isempty(opt), opt = struct(); end
def = struct('window',[], 'measure','rms', 'z',5, 'frac',0.25, 'min_trials',5, ...
             'action','remove', 'verbose',true);
f = fieldnames(def);
for i = 1:numel(f)
    if ~isfield(opt,f{i}) || (isempty(opt.(f{i})) && ~isempty(def.(f{i}))), opt.(f{i}) = def.(f{i}); end
end
if EEG.trials < 2
    error('ieeglab_reject_trials:notEpoched', 'Trial rejection needs epoched data (this dataset has %d trial).', EEG.trials);
end
if ~any(strcmpi(opt.action, {'remove','mark'}))
    error('ieeglab_reject_trials:badAction', 'action must be ''remove'' or ''mark''.');
end

t = double(EEG.times(:))';
if isempty(opt.window)
    tmask = t < -2 | t > 15;
else
    tmask = t >= opt.window(1) & t <= opt.window(2);
end
if nnz(tmask) < 5
    error('ieeglab_reject_trials:shortWindow', 'The amplitude window selects only %d samples.', nnz(tmask));
end

[sites, stimIdx] = ieeglab_epoch_sites(EEG);
[uSites, ~, grp] = unique(sites);
C = EEG.nbchan; N = EEG.trials;
isBadCh = false(1, C);
if isfield(EEG.chanlocs,'status')
    isBadCh = cellfun(@(x) ~isempty(x) && strcmpi(char(x),'bad'), {EEG.chanlocs.status});
end

X = double(EEG.data(:, tmask, :));
switch lower(opt.measure)
    case 'rms',    A = squeeze(sqrt(mean(X.^2, 2, 'omitnan')));   % C x N
    case 'maxabs', A = squeeze(max(abs(X), [], 2));
    otherwise, error('ieeglab_reject_trials:badMeasure', 'measure must be ''rms'' or ''maxabs''.');
end
if C == 1, A = A(:)'; end

nOver = zeros(N,1); fracOver = zeros(N,1); maxZ = zeros(N,1); assessed = false(N,1);
for g = 1:numel(uSites)
    tr = find(grp == g);
    if numel(tr) < opt.min_trials, continue; end
    use = ~isBadCh;
    use(unique(vertcat(stimIdx{tr}))) = false;       % never judge by the stimulated contacts
    if ~any(use), continue; end
    Ag = log(max(A(use, tr), eps));                   % amplitudes are skewed; work in log
    med = median(Ag, 2, 'omitnan');
    mad1 = 1.4826 * median(abs(Ag - med), 2, 'omitnan');
    Z = (Ag - med) ./ max(mad1, eps);
    over = abs(Z) > opt.z;
    nOver(tr) = sum(over, 1)';
    fracOver(tr) = mean(over, 1)';
    maxZ(tr) = max(abs(Z), [], 1)';
    assessed(tr) = true;
end
bad = assessed & fracOver >= opt.frac;

T = table((1:N)', cellstr(sites(:)), nOver, fracOver, maxZ, bad, ...
    'VariableNames', {'trial','site','n_channels_over','frac_over','max_abs_z','bad'});
EEG.ieeglab.rejected_trials = T;

if opt.verbose
    fprintf('[trials] %d/%d trials flagged (robust z > %g on >= %.0f%% of channels, %s).', ...
        nnz(bad), N, opt.z, 100*opt.frac, opt.measure);
    if any(~assessed), fprintf(' %d trial(s) in sites with < %d trials not assessed.', nnz(~assessed), opt.min_trials); end
    fprintf('\n');
end

if any(bad)
    if strcmpi(opt.action, 'mark')
        if ~isfield(EEG,'reject') || ~isstruct(EEG.reject), EEG.reject = struct(); end
        EEG.reject.rejmanual = bad(:)';
        if ~isfield(EEG.reject,'rejmanualE') || size(EEG.reject.rejmanualE,2) ~= N
            EEG.reject.rejmanualE = zeros(C, N);
        end
    else
        EEG = pop_select(EEG, 'notrial', find(bad));
        EEG.ieeglab.rejected_trials = T;   % indices refer to the trials before removal
    end
end
com = sprintf('EEG = ieeglab_reject_trials(EEG, struct(''z'',%g,''frac'',%g,''action'',''%s''));', ...
    opt.z, opt.frac, opt.action);
end
