function tests = test_ieeglab_icaref
% ICA re-referencing (ieeglab_icaref, Michelmann et al., 2018) on simulated
% contacts: local sources that spread to their neighbours, plus a reference
% signal added equally to every contact. Headless.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(mfilename('fullpath')));
tc.assumeNotEmpty(which('eeglab'), 'EEGLAB is not on the MATLAB path.');
if isempty(which('runica')), evalc('eeglab nogui'); end
addpath(root); addpath(fullfile(root, 'functions'));
end

function test_removes_the_shared_reference(tc)
[E, S, r] = local_sim(8, 0);
X0 = reshape(double(E.data), E.nbchan, []);
[~, E2, info] = evalc('ieeglab_icaref(E, struct(''verbose'', false))');
X1 = reshape(double(E2.data), E2.nbchan, []);
before = arrayfun(@(k) abs(corr(X0(k,:)', r)), 1:E.nbchan);
after  = arrayfun(@(k) abs(corr(X1(k,:)', r)), 1:E.nbchan);
own    = arrayfun(@(k) corr(X1(k,:)', S(k,:)'), 1:size(S, 1));   % contacts with a source
tc.verifyGreaterThanOrEqual(nnz(info.broad), 1, 'The reference must be found as a broad component.');
tc.verifyGreaterThan(min(before), 0.4, 'Simulation check: every contact carries the reference.');
tc.verifyLessThan(max(after), 0.1, 'The reference must be gone from every contact.');
tc.verifyGreaterThan(min(own), 0.9, 'Each contact must still carry its own local source.');
tc.verifyEqual(E2.ref, 'ICA re-referencing (Michelmann et al., 2018)');
end

function test_rank_is_estimated(tc)
% An average-referenced input has one dimension fewer than it has contacts
[E, ~, ~] = local_sim(8, 0);
X = double(E.data); X = X - mean(X, 1); E.data = single(X);
[~, ~, info] = evalc('ieeglab_icaref(E, struct(''verbose'', false))');
tc.verifyEqual(info.estimated_rank, 7);
tc.verifyEqual(info.rank, 7);
[~, ~, info] = evalc('ieeglab_icaref(E, struct(''rank'', 6, ''verbose'', false))');
tc.verifyEqual(info.rank, 6, 'An explicit rank must be used as given.');
end

function test_broadness_does_not_depend_on_units(tc)
% The same data in volts instead of microvolts must give the same broad
% components with normalized weights (the default).
[E, ~, ~] = local_sim(8, 0);
[~, ~, a] = evalc('ieeglab_icaref(E, struct(''verbose'', false))');
E.data = E.data * 1e-6;
[~, ~, b] = evalc('ieeglab_icaref(E, struct(''verbose'', false))');
tc.verifyEqual(b.broad, a.broad);
end

function test_bad_channels_are_left_out(tc)
[E, ~, ~] = local_sim(8, 0);
E.chanlocs(3).status = 'bad';
X0 = double(E.data(3, :, :));
[~, E2, info] = evalc('ieeglab_icaref(E, struct(''verbose'', false))');
tc.verifyEqual(double(E2.data(3, :, :)), X0, 'A bad contact must be left untouched.');
tc.verifyFalse(ismember('C3', info.channels));
end

function test_same_through_the_dialog_function(tc)
[E, ~, ~] = local_sim(8, 0);
[~, A] = evalc('ieeglab_icaref(E, struct(''verbose'', false))');
[~, B, com] = evalc('pop_ieeglab_reref(E, ''method'', ''ica'')');
tc.verifyEqual(double(B.data), double(A.data), 'AbsTol', 1e-4);
tc.verifyTrue(contains(com, '''method'', ''ica'''));
end

function test_stimulated_contacts_stay_out_of_the_others(tc)
% An artifact on the stimulated contacts, in the epochs of their site, must
% reach no other contact.
[E, ~, r] = local_sim(8, 0);
E = local_ccep(E);
tc.assumeEqual(ieeglab_detect_mode(E), 'ccep');
[~, A] = evalc('ieeglab_icaref(E, struct(''verbose'', false))');
post = E.times > 0;
tr = 1:2:E.trials;                                   % the epochs of site C2-C3
Bad = E;
Bad.data([2 3], post, tr) = Bad.data([2 3], post, tr) + 5000;
[~, B] = evalc('ieeglab_icaref(Bad, struct(''verbose'', false))');
others = [1 4:8];
tc.verifyEqual(double(B.data(others, :, :)), double(A.data(others, :, :)), 'AbsTol', 1e-2);
gap = double(B.data([2 3], post, tr)) - double(A.data([2 3], post, tr));
tc.verifyEqual(gap, 5000 * ones(size(gap)), 'AbsTol', 1e-2, 'They keep their own artifact.');
% The reference still leaves the epochs in which contacts were left out
Y = reshape(double(A.data(others, :, tr)), numel(others), []);
R = reshape(r, E.pnts, E.trials);
R = reshape(R(:, tr), [], 1);
tc.verifyLessThan(max(abs(corr(Y', R))), 0.25);
end

function test_blanked_samples_cost_no_contact(tc)
% Blanking by NaN writes NaN into every contact: not a reason to leave them out
[E, ~, ~] = local_sim(8, 0);
E = local_ccep(E);
gap = E.times >= 0 & E.times < 10;
E.data(:, gap, :) = NaN;
[~, A, info] = evalc('ieeglab_icaref(E, struct(''verbose'', false))');
tc.verifyNumElements(info.channels, 8);
tc.verifyTrue(all(isnan(A.data(:, gap, :)), 'all'));
tc.verifyTrue(all(isfinite(A.data(:, ~gap, :)), 'all'));
end

% -------------------------------------------------------------------------
function E = local_ccep(E)
% Two stimulation sites taking turns, as pop_epoch records them
sites = {'C2-C3', 'C5-C6'};
i0 = find(E.times >= 0, 1);
E.event = struct('type', {}, 'latency', {}, 'epoch', {});
for n = 1:E.trials
    s = sites{mod(n - 1, 2) + 1};
    E.event(n) = struct('type', s, 'latency', (n - 1) * E.pnts + i0, 'epoch', n);
    E.epoch(n).event = n; E.epoch(n).eventtype = {s}; E.epoch(n).eventlatency = {0};
end
end

function [E, S, r] = local_sim(nCh, seed)
% nCh-1 local sources (Laplacian, non-Gaussian) on contacts 1..nCh-1, each
% strongest on its own contact and decaying to the neighbours, plus one
% reference signal added equally to all contacts: as many sources as contacts,
% so ICA can separate them.
rng(seed);
srate = 500; nT = 500; nTr = 40; t = (0:nT-1) / srate * 1000 - 400;
n = nT * nTr; nS = nCh - 1;
S = sign(randn(nS, n)) .* -log(rand(nS, n)) * 10;             % Laplacian, 10 uV
S = filter(ones(1,3)/3, 1, S, [], 2);                          % mild smoothing
r = filter(ones(1,5)/5, 1, sign(randn(1, n)) .* rand(1, n).^3 * 30, [], 2)';
A = zeros(nCh, nS);
for k = 1:nS, A(:, k) = 0.3 .^ abs((1:nCh)' - k); end          % spread to neighbours
X = A * S + ones(nCh, 1) * r';
E = eeg_emptyset;
E.data = single(reshape(X, nCh, nT, nTr));
E.srate = srate; E.nbchan = nCh; E.pnts = nT; E.trials = nTr;
E.xmin = t(1)/1000; E.xmax = t(end)/1000; E.times = t; E.ref = 'common';
E.chanlocs = struct('labels', arrayfun(@(k) sprintf('C%d', k), 1:nCh, 'UniformOutput', false), ...
    'status', repmat({'good'}, 1, nCh));
end
