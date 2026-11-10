function out = validate_icaref_michelmann(dataDir)
% validate_icaref_michelmann() - iEEGLAB's ICA re-referencing (ieeglab_icaref)
%       against the authors' code (exampleICAreref.m; Michelmann et al., 2018),
%       on the depth leads of FieldTrip's intracranial tutorial dataset and on
%       simulated data with a known broad source.
%
% Usage:
%   out = validate_icaref_michelmann('C:\data\SubjectUCI29')
%
% dataDir holds SubjectUCI29_data.mat and SubjectUCI29_elec_acpc_f.mat
% (download.fieldtriptoolbox.org/tutorial/human_ecog/SubjectUCI29/). The
% authors' repository (github.com/smichelmann/ICA-rereferencing) has the
% script but no data; their method was shown on depth electrodes, so the 56
% contacts of the 7 depth leads are used here.
%
% The authors' script high-passes the epochs at 1 Hz, runs ICA through
% FieldTrip's ft_componentanalysis with its defaults (EEGLAB's runica, not
% extended; the sphering matrix divided by its norm, normalisesphere 'yes'),
% takes mixing = inv(unmixing), computes a chi-square statistic on each column
% of the absolute mixing matrix against chi2inv(0.8, contacts - 1) (lines
% 54-67: components above it are local, the others broad) and, for method 4a
% (lines 120-124, ft_rejectcomponent), projects the local components back.
% FieldTrip is not needed here: the decomposition is iEEGLAB's runica with
% FieldTrip's settings, and the authors' lines are applied to it after
% FieldTrip's sphere scaling. What is compared is the scaling, the test and
% the projection, not the ICA itself (runica's own result).
%   1. UCI29 depth leads, high-passed at 1 Hz (EEGLAB's FIR filter here; the
%      authors' script uses FieldTrip's default): iEEGLAB with 'weights',
%      'fieldtrip' and 'extended', 0 against the authors' lines; the counts of
%      broad components with the other scalings of the same decomposition; and
%      iEEGLAB's default (normalized weights, extended runica) on the
%      unfiltered data, with how much of the data it removes.
%   2. Simulation, 20 contacts: one source spread evenly over all contacts (a
%      reference) and 19 local ones. Both must find the broad source and give
%      the same re-referenced data, so the projection is compared where it
%      removes something.
% The authors apply the unmixing to the unfiltered data; iEEGLAB applies it to
% the data it was fitted on, so part 1 compares both on the high-passed data.
%
% Reference:
%   Michelmann, S., Treder, M. S., Griffiths, B., Kerren, C., Roux, F.,
%   Wimber, M., Rollings, D., Sawlani, V., Chelvarajah, R., Gollwitzer, S.,
%   Kreiselmeyer, G., Hamer, H., Bowman, H., Staresina, B., & Hanslmayr, S.
%   (2018). Data-driven re-referencing of intracranial EEG based on independent
%   component analysis (ICA). Journal of Neuroscience Methods, 307, 125-137.
%
% Cedric Cannard, iEEGLAB, 2026

t0 = tic;
f = @(n) fullfile(dataDir, ['SubjectUCI29_' n '.mat']);
EEG = ieeglab_import_fieldtrip(f('data'), 'elec', f('elec_acpc_f'), 'types', {'*', 'SEEG'}, 'verbose', false);
depth = ~startsWith({EEG.chanlocs.labels}, {'LPG', 'LTG'});
D = pop_select(EEG, 'channel', find(depth));
[~, Dhp] = evalc('pop_eegfiltnew(D, ''locutoff'', 1, ''plotfreqz'', 0)');

% ---------------- 1. UCI29, the authors' path ----------------
[R1, i1] = ieeglab_icaref(Dhp, struct('weights', 'fieldtrip', 'extended', 0, 'verbose', false));
c = local_authors(i1, reshape(double(Dhp.data), Dhp.nbchan, []));
Y1 = reshape(double(R1.data), Dhp.nbchan, []);
out.n_contacts = Dhp.nbchan;
out.n_components = size(i1.mixing, 2);
out.sphere_norm = i1.sphere_norm;
out.broad_ieeglab = nnz(i1.broad);
out.broad_authors = nnz(~c.local);
out.same_components = isequal(i1.broad(:), ~c.local(:));
out.chi2_max_rel_diff = max(abs(i1.chi2(:) - c.chi2(:))) / max(abs(c.chi2));
out.chi2_min_over_threshold = min(c.chi2) / c.threshold;
out.data_max_rel_diff = max(abs(Y1(:) - c.X4a(:))) / max(abs(c.X4a(:)));
aw = abs(i1.mixing);                                        % other scalings, same decomposition
out.broad_raw_same_decomposition = nnz(local_p(aw, out.n_contacts) > 0.2);
out.broad_normalized_same_decomposition = nnz(local_p(aw ./ mean(aw, 1), out.n_contacts) > 0.2);
fprintf(['1. UCI29 depth leads, %d contacts, %d components (1 Hz high-pass, runica not extended, FieldTrip''s ' ...
    'sphere scaling, norm %.3g): broad iEEGLAB %d, authors'' code %d, same components %d; chi2 max relative ' ...
    'difference %.1e (smallest chi2 %.2f x the threshold); re-referenced data %.1e. Same decomposition, raw ' ...
    'weights: %d broad; normalized: %d.\n'], out.n_contacts, out.n_components, out.sphere_norm, ...
    out.broad_ieeglab, out.broad_authors, out.same_components, out.chi2_max_rel_diff, ...
    out.chi2_min_over_threshold, out.data_max_rel_diff, out.broad_raw_same_decomposition, ...
    out.broad_normalized_same_decomposition);

[R2, i2] = ieeglab_icaref(D, struct('verbose', false));     % iEEGLAB's default, unfiltered data
X = reshape(double(D.data), D.nbchan, []);
Y2 = reshape(double(R2.data), D.nbchan, []);
out.default_broad = nnz(i2.broad);
out.default_removed_rms_rel = sqrt(mean((X(:) - Y2(:)).^2)) / sqrt(mean(X(:).^2));
fprintf(['   iEEGLAB''s default (normalized weights, extended runica, unfiltered data): %d of %d broad, ' ...
    'removing %.3f of the data''s RMS.\n'], out.default_broad, numel(i2.broad), out.default_removed_rms_rel);

% ---------------- 2. simulation with a known broad source ----------------
rs = RandStream('mt19937ar', 'Seed', 3);
nC = 20; T = 20000;
Sx = local_laplace(rs, nC, T);                              % sources, unit variance, super-Gaussian
M = 0.15 * randn(rs, nC, nC) + 3 * eye(nC);                 % local sources: mostly on their own contact
M(:, 1) = 4 * (1 + 0.02 * randn(rs, nC, 1));                % source 1 on every contact alike
Sim = eeg_emptyset();
Sim.data = M * Sx; Sim.nbchan = nC; Sim.pnts = T; Sim.trials = 1; Sim.srate = 1000;
Sim.xmin = 0; Sim.xmax = (T - 1) / 1000; Sim.times = (0:T-1);
Sim.chanlocs = struct('labels', compose('S%d', 1:nC), 'type', 'SEEG');
[R3, i3] = ieeglab_icaref(Sim, struct('weights', 'fieldtrip', 'verbose', false));
c3 = local_authors(i3, double(Sim.data));
[~, kBroad] = max(abs(corr_cols(i3.unmixing * double(Sim.data), Sx(1, :))));
out.sim_broad_ieeglab = find(i3.broad);
out.sim_broad_authors = find(~c3.local);
out.sim_same_components = isequal(i3.broad(:), ~c3.local(:));
out.sim_broad_is_the_reference = isequal(out.sim_broad_ieeglab, kBroad);
out.sim_data_max_rel_diff = max(abs(double(R3.data(:)) - c3.X4a(:))) / max(abs(c3.X4a(:)));
fprintf(['2. Simulation, %d contacts, one source on all of them: broad iEEGLAB %s, authors'' code %s (the ' ...
    'simulated reference: component %d), same components %d; re-referenced data %.1e.\n'], nC, ...
    mat2str(out.sim_broad_ieeglab), mat2str(out.sim_broad_authors), kBroad, out.sim_same_components, ...
    out.sim_data_max_rel_diff);
out.seconds = toc(t0);
fprintf('Done in %.0f s.\n', out.seconds);
end

function c = local_authors(info, X)
% The authors' lines on iEEGLAB's decomposition after FieldTrip's sphere
% scaling (unmixing / norm(sphere)); mixing as ft_componentanalysis builds it.
unmixing = info.unmixing / info.sphere_norm;
if size(unmixing, 1) == size(unmixing, 2) && rank(unmixing) == size(unmixing, 1)
    mixing = inv(unmixing);
else
    mixing = pinv(unmixing);
end
% exampleICAreref.m lines 54-67, as written
chi_val = sum(bsxfun(@rdivide, ((bsxfun(@minus, abs(mixing), mean(abs(mixing), 1))).^2), ...
    mean(abs(mixing), 1)), 1);
[chisort, s_indices] = sort(chi_val, 'descend');
c.threshold = local_chi2inv(0.8, size(mixing, 1) - 1);
crit_20 = numel(find(chisort > c.threshold));
c.local = false(1, size(mixing, 2));
c.local(s_indices(1:crit_20)) = true;
c.chi2 = chi_val;
% lines 120-124 (ft_rejectcomponent): the local components projected back
c.X4a = mixing(:, c.local) * (unmixing(c.local, :) * X);
end

function p = local_p(aw, n)
e = mean(aw, 1);
p = gammainc(sum((aw - e).^2, 1) ./ e / 2, (n - 1) / 2, 'upper');
end

function r = corr_cols(A, b)
% correlation of each row of A with the row vector b (base MATLAB)
A = A - mean(A, 2); b = b - mean(b);
r = (A * b') ./ (sqrt(sum(A.^2, 2)) * norm(b));
end

function S = local_laplace(rs, n, T)
u = rand(rs, n, T) - 0.5;
S = -sign(u) .* log(1 - 2 * abs(u)) / sqrt(2);
end

function x = local_chi2inv(p, df)
% chi2inv (Statistics and Machine Learning Toolbox, used by the authors' code),
% or the same through base MATLAB
if exist('chi2inv', 'file') == 2
    x = chi2inv(p, df);
else
    x = 2 * gammaincinv(p, df / 2);
end
end
