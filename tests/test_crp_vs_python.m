function tests = test_crp_vs_python
% Cross-check run_CRP.m against the Python port used for the 2026-09-20 worklog
% numbers. Compute-only, no figures, safe under MATLAB -batch.
%
%   results = runtests('tests/test_crp_vs_python.m')
tests = functiontests(localfunctions);
end

function setupOnce(tc)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'functions'));
tc.TestData.dir = fullfile(here, 'reference');
end

function test_matches_python_port(tc)
d = tc.TestData.dir;
V   = readmatrix(fullfile(d, 'crp_test_input.csv'));     % time x trials
t   = readmatrix(fullfile(d, 'crp_test_time.csv'));      % seconds
ref = jsondecode(fileread(fullfile(d, 'crp_python_reference.json')));

[parms, projs] = run_CRP(V, t(:)', struct('verbose', false, 't_step', ref.t_step));

fprintf('\n  tR            MATLAB %.4f s   Python %.4f s\n', parms.tR, ref.tR);
fprintf('  mean expl var MATLAB %.4f     Python %.4f\n', mean(parms.expl_var), ref.mean_expl);

tc.verifyEqual(parms.tR, ref.tR, 'AbsTol', 1/ref.srate, ...
    'response duration tau_R differs by more than one sample');
tc.verifyEqual(mean(parms.expl_var), ref.mean_expl, 'AbsTol', 1e-6, ...
    'mean explained variance differs');
tc.verifyEqual(parms.expl_var(:), ref.expl_per_trial(:), 'AbsTol', 1e-6, ...
    'per-trial explained variance differs');

% canonical shape is defined up to a sign
c = parms.C(1:10); c = c * sign(dot(c, ref.canonical_first10));
tc.verifyEqual(c(:), ref.canonical_first10(:), 'AbsTol', 1e-8, ...
    'canonical shape differs');

al = parms.al(:) * sign(dot(parms.C(1:10), ref.canonical_first10));
tc.verifyEqual(al, ref.alpha(:), 'AbsTol', 1e-6, 'alpha weights differ');

% projection profile should be a sane length and finite
tc.verifyTrue(all(isfinite(projs.mean_proj_profile)), 'projection profile has non-finite values');
end
