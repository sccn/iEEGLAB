function tests = test_ieeglab_reject_trials
% Tests for automatic outlier-trial rejection. Headless.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(mfilename('fullpath')));
tc.assumeNotEmpty(which('eeglab'), 'EEGLAB is not on the MATLAB path.');
if isempty(which('pop_loadset')), evalc('eeglab nogui'); end
addpath(root); addpath(fullfile(root, 'functions'));
d = fullfile(root, 'tutorial', 'dataset_seeg');
EEG = pop_loadset('filename','sub-02_ses-ieeg01_task-ccep_run-01_ieeg.set','filepath',d);
[~, EEG] = evalc('ieeglab_load(EEG, struct())');
opt = struct('apply_highpass',true,'highpass',0.5,'apply_notch',false,'apply_lowpass',false, ...
    'apply_ds',false,'apply_epoch',true,'epoch_window',[-500 1000],'apply_car',false, ...
    'apply_baseline',true,'baseline_period',[-500 -50],'plot',false,'verbose',false);
[~, tc.TestData.E] = evalc('ieeglab_preprocess(EEG, opt)');
end

function test_injected_artifacts_are_found(tc)
E = tc.TestData.E;
sites = ieeglab_epoch_sites(E);
tr = find(sites == "ROP1-ROP2");            % the site with the most trials
% Real data can contain real outlier trials, so compare against the flags on
% the unmodified data rather than assuming none: injecting artifacts must add
% exactly the injected trials and leave every other decision unchanged.
[~, ~, T0] = evalc('ieeglab_reject_trials(E, struct(''action'',''mark'',''verbose'',false))');
bad = setdiff(tr, find(T0.bad));
bad = bad([2 5 9]);
E.data(:, :, bad) = E.data(:, :, bad) * 40;  % whole-trial artifact on every channel
[~, ~, T] = evalc('ieeglab_reject_trials(E, struct(''action'',''mark'',''verbose'',false))');
tc.verifyTrue(all(T.bad(bad)), 'Every injected artifact trial must be flagged.');
others = setdiff(1:E.trials, bad);
tc.verifyEqual(T.bad(others), T0.bad(others), ...
    'Injecting artifacts into three trials must not change the decision on any other trial.');
end

function test_mark_uses_eeglab_rejmanual_and_remove_drops(tc)
E = tc.TestData.E;
sites = ieeglab_epoch_sites(E);
tr = find(sites == "ROP1-ROP2");
E.data(:, :, tr(3)) = E.data(:, :, tr(3)) * 40;
[~, Em] = evalc('ieeglab_reject_trials(E, struct(''action'',''mark'',''verbose'',false))');
tc.verifyEqual(Em.trials, E.trials, 'mark must not remove trials.');
tc.verifyTrue(logical(Em.reject.rejmanual(tr(3))), 'mark must set EEG.reject.rejmanual.');
[~, Er] = evalc('ieeglab_reject_trials(E, struct(''action'',''remove'',''verbose'',false))');
tc.verifyEqual(Er.trials, E.trials - nnz(Em.reject.rejmanual));
end

function test_large_response_on_few_channels_is_kept(tc)
% A genuinely strong response on two contacts is not an artifact.
E = tc.TestData.E;
sites = ieeglab_epoch_sites(E);
tr = find(sites == "ROP1-ROP2");
E.data(1:2, :, tr(4)) = E.data(1:2, :, tr(4)) * 40;
[~, ~, T] = evalc('ieeglab_reject_trials(E, struct(''action'',''mark'',''verbose'',false))');
tc.verifyFalse(T.bad(tr(4)), 'Two outlying channels out of 14 usable must not reject the trial.');
end

function test_preprocess_option_and_clinician_agreement(tc)
% The tutorial events file marks two stimulation events as artifacts
% (status_description "... fixed when referencing"). The detector, which never
% sees those notes, must flag both - and preprocess must remove exactly the
% flagged trials when reject_trials is on.
root = fileparts(fileparts(mfilename('fullpath')));
d = fullfile(root, 'tutorial', 'dataset_seeg');
EEG = pop_loadset('filename','sub-02_ses-ieeg01_task-ccep_run-01_ieeg.set','filepath',d);
[~, EEG] = evalc('ieeglab_load(EEG, struct())');
ev = EEG.ieeglab.opt.events;
annotated = find(contains(string(ev.status_description), "artif"));
tc.assumeNumElements(annotated, 2, 'Expected two clinician-annotated artifact events in the tutorial.');
opt = struct('apply_highpass',true,'highpass',0.5,'apply_notch',false,'apply_lowpass',false, ...
    'apply_ds',false,'apply_epoch',true,'epoch_window',[-500 1000],'apply_car',false, ...
    'apply_baseline',true,'baseline_period',[-500 -50],'plot',false,'verbose',false);
[~, E] = evalc('ieeglab_preprocess(EEG, opt)');
tr0 = arrayfun(@(r) local_trial_of_urevent(E, r), annotated);
tc.assertFalse(any(isnan(tr0)), 'Annotated events must be traceable to trials through urevent.');
[~, ~, T] = evalc('ieeglab_reject_trials(E, struct(''action'',''mark'',''verbose'',false))');
tc.verifyTrue(all(T.bad(tr0)), sprintf('Both clinician-annotated artifact trials (%s) should be flagged.', mat2str(tr0(:)')));
opt.reject_trials = true;
[~, E2] = evalc('ieeglab_preprocess(EEG, opt)');
tc.verifyEqual(E2.trials, E.trials - nnz(T.bad), 'reject_trials in preprocess must remove exactly the flagged trials.');
end

function k = local_trial_of_urevent(E, r)
% Trial whose time-locking event (latency 0) is original event r.
k = NaN;
for t = 1:E.trials
    ue = E.epoch(t).eventurevent; lat = E.epoch(t).eventlatency;
    if ~iscell(ue), ue = {ue}; lat = {lat}; end
    [~, j] = min(cellfun(@(x) abs(double(x)), lat));
    if isequal(double(ue{j}), r), k = t; return; end
end
end

function test_stimulated_contacts_do_not_drive_rejection(tc)
E = tc.TestData.E;
sites = ieeglab_epoch_sites(E);
tr = find(sites == "ROP1-ROP2");
lab = {E.chanlocs.labels};
st = find(ismember(lab, {'ROP1','ROP2'}));
E.data(st, :, tr(6)) = E.data(st, :, tr(6)) * 1000;   % absurd, but only on the stimulated pair
[~, ~, T] = evalc('ieeglab_reject_trials(E, struct(''action'',''mark'',''verbose'',false,''frac'',0.1))');
tc.verifyFalse(T.bad(tr(6)), 'The stimulated contacts must be ignored when judging a trial.');
end
