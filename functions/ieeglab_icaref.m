function [EEG, info] = ieeglab_icaref(EEG, opt)
% ieeglab_icaref() - ICA-based re-referencing of iEEG (Michelmann et al., 2018).
%
% Usage:
%   EEG = ieeglab_icaref(EEG)
%   [EEG, info] = ieeglab_icaref(EEG, opt)
%
% Independent component analysis unmixes the linear superposition of sources on
% the contacts, which is what re-referencing tries to undo. Following the
% procedure of Michelmann et al. (2018, steps 1-3 and 4a):
%   1. ICA (EEGLAB's runica) on the good contacts;
%   2. for every column of the absolute mixing matrix (the component's weights
%      on the contacts), a chi-square statistic of its departure from a uniform
%      distribution, with as many degrees of freedom as contacts;
%   3. components whose weights are uniform (p > opt.p_broad, default 0.2) are
%      "broad": the reference and other shared signals. They are discarded;
%   4. the local components are projected back to the contacts.
% The result keeps the channel representation; its rank drops by the number of
% broad components.
%
% Options (fields of opt):
%   .rank        number of components. Default [] = the effective rank of the
%                data used for the fit, estimated from the eigenvalues of its
%                covariance (a common average or a derivation lowers it). Give it
%                explicitly when it is known.
%   .p_broad     chi-square p-value above which a component counts as broad.
%                Default 0.2, as in the paper.
%   .weights     'normalized' (default) or 'raw'. The paper computes the
%                chi-square on the absolute mixing weights but does not say how
%                they are scaled, and on raw weights the statistic grows with the
%                data's units: on microvolt iEEG every component tests as local
%                (0 of 15 broad on the sEEG tutorial data). 'normalized' divides
%                each column by its mean first, so the test measures how unevenly
%                the component is spread over the contacts (chi2 = n x CV^2 of the
%                weights), whatever the units.
%   .fit_window  [start end] ms of each epoch used to estimate the unmixing.
%                Default for CCEP data: the pre-stimulus part, from the epoch
%                start to -50 ms, so the stimulation artifact cannot drive the
%                decomposition; for other data, the whole epoch. The unmixing
%                is then applied to every sample.
%   .extended    runica 'extended' option. Default 1 (as pop_runica).
%   .reproducible  true (default) runs runica with its random seed fixed
%                ('rndreset' 'no'), so the same data give the same result; false
%                lets runica seed from the clock, as by default in EEGLAB.
%   .bad_channels  indices or labels to leave out (untouched). Channels with
%                chanlocs.status 'bad' are always left out.
%   .verbose     default true
%
% Output info: rank, n_components, broad (logical per component), chi2, p,
% channels (the contacts used), fit_window_ms.
%
% Reference:
%   Michelmann, S., Treder, M. S., Griffiths, B., Kerren, C., Roux, F.,
%   Wimber, M., Rollings, D., Sawlani, V., Chelvarajah, R., Gollwitzer, S.,
%   Kreiselmeyer, G., Hamer, H., Bowman, H., Staresina, B., & Hanslmayr, S.
%   (2018). Data-driven re-referencing of intracranial EEG based on independent
%   component analysis (ICA). Journal of Neuroscience Methods, 307, 125-137.
%
% Cedric Cannard, iEEGLAB, 2026

if nargin < 2 || isempty(opt), opt = struct(); end
def = struct('rank', [], 'p_broad', 0.2, 'weights', 'normalized', 'fit_window', [], 'extended', 1, 'reproducible', true, ...
    'bad_channels', [], 'verbose', true);
for f = fieldnames(def)'
    if ~isfield(opt, f{1}), opt.(f{1}) = def.(f{1}); end
end
if isempty(which('runica'))
    error('ieeglab_icaref:noRunica', 'runica (EEGLAB) is not on the path. Start EEGLAB first.');
end
if ~isfield(EEG, 'trials') || EEG.trials < 2
    error('ieeglab_icaref:notEpoched', 'ICA re-referencing runs on epoched data.');
end

% ---- channels: good contacts only
labels = {EEG.chanlocs.labels};
bad = false(1, EEG.nbchan);
if isfield(EEG.chanlocs, 'status')
    bad = cellfun(@(x) ~isempty(x) && strcmpi(char(string(x)), 'bad'), {EEG.chanlocs.status});
end
if ~isempty(opt.bad_channels)
    if isnumeric(opt.bad_channels), bad(opt.bad_channels) = true;
    else, bad(ismember(lower(labels), lower(cellstr(opt.bad_channels)))) = true; end
end
X = double(EEG.data);
bad = bad | squeeze(any(any(~isfinite(X), 2), 3))';
good = find(~bad);
nGood = numel(good);
if nGood < 3
    error('ieeglab_icaref:tooFewChannels', 'ICA re-referencing needs at least 3 good contacts (%d here).', nGood);
end

% ---- samples used to estimate the unmixing
t = double(EEG.times);
if isempty(opt.fit_window)
    if strcmp(ieeglab_detect_mode(EEG), 'ccep')
        opt.fit_window = [t(1) -50];
    else
        opt.fit_window = [t(1) t(end)];
    end
end
fitIdx = t >= opt.fit_window(1) & t <= opt.fit_window(2);
if nnz(fitIdx) < 10
    error('ieeglab_icaref:fitWindow', 'The fit window [%g %g] ms holds fewer than 10 samples.', opt.fit_window);
end
Xfit = reshape(X(good, fitIdx, :), nGood, []);

% ---- effective rank of the data used for the fit
ev = sort(eig(cov(Xfit')), 'descend');
estRank = nnz(ev > ev(1) * 1e-7);
r = opt.rank;
if isempty(r), r = estRank; end
r = min(r, nGood);
if r > estRank
    warning('ieeglab_icaref:rankAboveEstimate', ...
        'Requested %d components but the data have an estimated rank of %d; ICA may be unstable.', r, estRank);
end
nSamp = size(Xfit, 2);
if nSamp < 20 * r^2 && opt.verbose
    fprintf('[ICA ref] %d samples for %d components (fewer than the 20 x rank^2 often recommended).\n', nSamp, r);
end

% ---- ICA
% runica works on small matrices, where multithreaded BLAS is much slower than
% one thread (about 100 times on a busy 16-thread machine); restore afterwards.
% runica switches the global generator to its own legacy state: restore it too.
prevRng = rng;
prevThreads = maxNumCompThreads(1);
cleanup = onCleanup(@() local_restore(prevRng, prevThreads));
args = {'extended', opt.extended, 'verbose', 'off', 'rndreset', iff(opt.reproducible, 'no', 'yes')};
if r < nGood, args = [args {'pca', r}]; end
[~, w, s] = evalc('runica(Xfit, args{:})');
U = w * s;                          % r x nGood unmixing
A = pinv(U);                        % nGood x r mixing: column k = component k on the contacts

% ---- broadness: chi-square of |mixing weights| against their mean, df = contacts
aw = abs(A);
switch lower(char(opt.weights))
    case 'normalized', aw = aw ./ mean(aw, 1);           % mean 1 per column: unit-free
    case 'raw'                                           % as the weights come (units of the data)
    otherwise, error('ieeglab_icaref:weights', 'weights must be ''normalized'' or ''raw''.');
end
e = mean(aw, 1);
chi2 = sum((aw - e).^2, 1) ./ e;
p = 1 - gammainc(chi2 / 2, nGood / 2);
broad = p > opt.p_broad;
if all(broad)
    warning('ieeglab_icaref:allBroad', 'All %d components are broad (p > %g); the data are left unchanged.', r, opt.p_broad);
    broad(:) = false;
    keepAll = true;
else
    keepAll = false;
end

% ---- back-project the local components onto the contacts, every sample
if ~keepAll
    P = A(:, ~broad) * U(~broad, :);        % nGood x nGood projector
    Y = reshape(X(good, :, :), nGood, []);
    X(good, :, :) = reshape(P * Y, nGood, EEG.pnts, EEG.trials);
end
EEG.data = cast(X, 'like', EEG.data);
EEG.ref = 'ICA re-referencing (Michelmann et al., 2018)';
info = struct('method', 'ica', 'rank', r, 'estimated_rank', estRank, 'n_components', r, ...
    'broad', broad, 'chi2', chi2, 'p', p, 'p_broad', opt.p_broad, 'weights', char(opt.weights), ...
    'mixing', A, 'channels', {labels(good)}, ...
    'fit_window_ms', opt.fit_window, 'reproducible', opt.reproducible);
if ~isfield(EEG, 'ieeglab'), EEG.ieeglab = struct(); end
EEG.ieeglab.car = info;
if opt.verbose
    fprintf(['[ICA ref] %d good contacts, rank %d (estimated %d); %d of %d components broad ' ...
        '(p > %g) removed; fit on %g to %g ms.\n'], nGood, r, estRank, nnz(broad), r, opt.p_broad, opt.fit_window);
end
end

function local_restore(prevRng, prevThreads)
rng(prevRng);
maxNumCompThreads(prevThreads);
end

function out = iff(c, a, b)
if c, out = a; else, out = b; end
end
