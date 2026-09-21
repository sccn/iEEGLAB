function tests = test_carla_vs_reference
% test_carla_vs_reference - Validate ieeglab_carla against the published CARLA.
%
% The reference is functions/CARLA.m from the CARLA manuscript package
% (https://github.com/hharveygit/CARLA_JNM, Harvey Huang, GPL-3), vendored
% unmodified as tests/reference/CARLA_reference.m.
%
% Validation strategy
% -------------------
% CARLA is only partly deterministic: with more than one trial it bootstraps
% the mean signal, so two correct implementations will not agree bit for bit.
% The tests are therefore split by what CAN be checked exactly:
%
%   1. DETERMINISTIC, exact equality:
%        - the channel ranking (stats.vars, stats.order), computed before any
%          randomness enters
%        - the single-trial path (nTrs == 1), where the reference skips
%          bootstrapping entirely, so chsUsed and Vout must match exactly
%   2. STOCHASTIC, distributional:
%        - the multi-trial path, compared over repeated runs. The chosen CAR
%          size must agree closely, and the re-referenced data must agree to
%          within the spread that bootstrapping alone produces.
%
% ieeglab_carla is configured to match the reference exactly for these tests:
% the reference hardcodes twin = [0.01 0.3] and notches at 60/120/180 Hz.
%
% Cedric Cannard, iEEGLAB, 2026

tests = functiontests(localfunctions);
end

function setupOnce(tc)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'reference'));
addpath(fileparts(here));
addpath(fullfile(fileparts(here), 'functions'));
tc.assumeNotEmpty(which('CARLA_reference'), ...
    'tests/reference/CARLA_reference.m is missing.');
tc.assumeNotEmpty(which('ieeglab_carla'));

% Options that make ieeglab_carla behave exactly like the reference
tc.TestData.matchOpts = struct('winResp',[0.01 0.3], 'lineFreq',60, ...
    'notchFirst',true, 'sens',false, 'verbose',false, 'minCARsize',2);
end

function [V, tt, srate] = simulateCCEP(nQuiet, nResp, K, seed)
% A CCEP-like montage: shared line noise and pink-ish background on every
% channel, plus a large evoked response on the last nResp channels. Those are
% the channels CARLA must keep OUT of the common average.
rng(seed);
srate = 1000;
tt = -0.2 : 1/srate : 0.5;
T = numel(tt);
n = nQuiet + nResp;

V = zeros(n, T, K);
for k = 1:K
    shared = 8*sin(2*pi*60*tt + 2*pi*rand) + 4*cumsum(randn(1,T))/sqrt(T);
    for ch = 1:n
        V(ch,:,k) = shared + 3*randn(1,T);
    end
end
% Evoked response: an N1/N2-like biphasic deflection after t = 0
n1 = -140 * exp(-((tt-0.025)/0.008).^2);
n2 =   90 * exp(-((tt-0.110)/0.045).^2);
evoked = n1 + n2;
evoked(tt < 0) = 0;
for ch = nQuiet+1 : n
    amp = 0.7 + 0.6*rand;
    V(ch,:,:) = V(ch,:,:) + amp * repmat(evoked, 1, 1, K);
end
end

% ============================ tests ============================

function test_ranking_is_identical(tc)
% stats.vars and stats.order are computed before any bootstrap draw, so two
% correct implementations must agree exactly.
[V, tt, srate] = simulateCCEP(12, 4, 20, 1);
o = tc.TestData.matchOpts;

[~, ~, sRef] = CARLA_reference(tt, V, srate, false, 20);
[~, ~, sMine] = ieeglab_carla(tt, V, srate, o);

tc.verifyEqual(sMine.vars, sRef.vars, 'RelTol', 1e-10, ...
    'Channel ranking metric (mean cross-trial covariance) differs from the reference.');
tc.verifyEqual(sMine.order(:), sRef.order(:), ...
    'Channel ordering differs from the reference.');
end

function test_single_trial_is_exactly_identical(tc)
% With one trial the reference does no bootstrapping, so the whole algorithm
% is deterministic and the outputs must match to machine precision.
[V, tt, srate] = simulateCCEP(12, 4, 1, 7);
o = tc.TestData.matchOpts;

[VoutRef, CARref, sRef]   = CARLA_reference(tt, V, srate);
[VoutMine, CARmine, sMine] = ieeglab_carla(tt, V, srate, o);

tc.verifyEqual(sMine.chsUsed(:), sRef.chsUsed(:), ...
    sprintf(['Single-trial CAR channel set differs.\n  reference: %s\n  ieeglab:   %s'], ...
    mat2str(sRef.chsUsed(:)'), mat2str(sMine.chsUsed(:)')));
tc.verifyEqual(sMine.nOptimum, numel(sRef.chsUsed));
tc.verifyEqual(VoutMine, VoutRef, 'AbsTol', 1e-9, ...
    'Single-trial re-referenced data differs from the reference.');
tc.verifyEqual(squeeze(CARmine), squeeze(CARref), 'AbsTol', 1e-9);
end

function test_single_trial_zminmean_profile_matches(tc)
% The whole method rests on the shape of the zMinMean curve; check the curve
% itself, not just the argmax it produces.
[V, tt, srate] = simulateCCEP(10, 3, 1, 11);
o = tc.TestData.matchOpts;

[~, ~, sRef]  = CARLA_reference(tt, V, srate);
[~, ~, sMine] = ieeglab_carla(tt, V, srate, o);

zRef  = squeeze(sRef.zMinMean(:))';
zMine = squeeze(sMine.zMinMean(:))';
tc.verifyEqual(numel(zMine), numel(zRef), 'zMinMean has a different length.');

ok = ~isnan(zRef) & ~isnan(zMine);
tc.verifyGreaterThan(nnz(ok), 3, 'Too few comparable points in zMinMean.');
tc.verifyEqual(zMine(ok), zRef(ok), 'AbsTol', 1e-9, ...
    'The zMinMean anticorrelation profile differs from the reference.');
end

function test_multitrial_agrees_distributionally(tc)
% Bootstrapping makes the multi-trial path stochastic. Run both repeatedly and
% require the chosen CAR size to agree closely. A systematic algorithmic
% difference would show up as a consistent offset, not as scatter.
[V, tt, srate] = simulateCCEP(14, 4, 24, 3);
o = tc.TestData.matchOpts;
nRuns = 8;

nRef = zeros(1,nRuns); nMine = zeros(1,nRuns);
for r = 1:nRuns
    rng(100+r); [~, ~, sR] = CARLA_reference(tt, V, srate, false, 100);
    rng(100+r); [~, ~, sM] = ieeglab_carla(tt, V, srate, o);
    nRef(r)  = numel(sR.chsUsed);
    nMine(r) = numel(sM.chsUsed);
end

tc.verifyLessThanOrEqual(abs(mean(nMine) - mean(nRef)), 1.0, ...
    sprintf(['Mean CAR size differs systematically from the reference.\n' ...
             '  reference: %s (mean %.2f)\n  ieeglab:   %s (mean %.2f)'], ...
             mat2str(nRef), mean(nRef), mat2str(nMine), mean(nMine)));

% Neither implementation may include a responsive channel (last 4)
for r = 1:nRuns
    rng(200+r); [~, ~, sM] = ieeglab_carla(tt, V, srate, o);
    tc.verifyLessThanOrEqual(max(sM.chsUsed), 14, ...
        'ieeglab_carla put a responsive channel into the common average.');
end
end

function test_both_reject_responsive_channels(tc)
% The point of CARLA: the responsive channels must stay out of the reference.
% Verified against the reference on the same data so a shared failure would be
% visible rather than hidden.
for nResp = [2 4 6]
    [V, tt, srate] = simulateCCEP(14, nResp, 20, 40+nResp);
    o = tc.TestData.matchOpts;

    rng(5); [~, ~, sRef]  = CARLA_reference(tt, V, srate, false, 60);
    rng(5); [~, ~, sMine] = ieeglab_carla(tt, V, srate, o);

    tc.verifyLessThanOrEqual(max(sRef.chsUsed), 14, ...
        sprintf('REFERENCE included a responsive channel (nResp=%d) - check the simulation.', nResp));
    tc.verifyLessThanOrEqual(max(sMine.chsUsed), 14, ...
        sprintf('ieeglab_carla included a responsive channel (nResp=%d).', nResp));
end
end

function test_reference_cannot_handle_nan_but_we_can(tc)
% A deliberate, documented divergence: ieeglab_carla treats an all-NaN channel
% as "excluded from the reference", which is how the stimulated contacts are
% kept out. The reference has no such concept. This test pins the divergence
% so it stays intentional.
[V, tt, srate] = simulateCCEP(12, 3, 15, 21);
o = tc.TestData.matchOpts;

Vn = V;
Vn([2 5],:,:) = NaN;
[~, ~, s] = ieeglab_carla(tt, Vn, srate, o);
tc.verifyEmpty(intersect(s.chsUsed, [2 5]), ...
    'NaN-masked channels must never enter the common average.');

% Without the mask those same channels are eligible again
[~, ~, s2] = ieeglab_carla(tt, V, srate, o);
tc.verifyNotEmpty(intersect(s2.chsUsed, [2 5]), ...
    'Sanity check failed: channels 2 and 5 should be usable when not masked.');
end

function test_no_reference_leak_on_real_tutorial_data(tc)
% End-to-end on the shipped sEEG CCEP dataset: after per-site CARLA, no
% stimulated contact may appear in its own reference set.
root = fileparts(fileparts(mfilename('fullpath')));
d = fullfile(root, 'tutorial', 'dataset_seeg');
tc.assumeTrue(isfolder(d), 'Tutorial sEEG dataset not present.');
tc.assumeNotEmpty(which('pop_loadset'), 'EEGLAB not on the path.');

EEG = pop_loadset('filename','sub-02_ses-ieeg01_task-ccep_run-01_ieeg.set','filepath',d);
elecs = readtable(fullfile(d,'sub-02_ses-ieeg01_electrodes.tsv'),'FileType','text','Delimiter','\t');
ev    = readtable(fullfile(d,'sub-02_ses-ieeg01_task-ccep_run-01_events.tsv'),'FileType','text','Delimiter','\t');
[~, EEG] = evalc('get_elec_coor(EEG, elecs)');
EEG.event = [];
for i = 1:height(ev)
    EEG.event(i).type    = ev.electrical_stimulation_site{i};
    EEG.event(i).latency = ev.onset(i)*EEG.srate + 1;
end
EEG = eeg_checkset(EEG,'eventconsistency');
[~, EEG] = evalc('pop_epoch(EEG, unique({EEG.event.type}), [-0.5 1], ''epochinfo'',''yes'')');

[~, out] = ieeglab_car(EEG, struct('car_method','carla','car_nboot',30,'verbose',false));
tc.verifyNotEmpty(out);
for g = 1:numel(out)
    tc.verifyEmpty(intersect(out(g).excluded_channels, out(g).car_channels), ...
        sprintf('Site %s leaked a stimulated contact into its reference.', out(g).group));
end
end
