function [EEG, out] = ieeglab_car(EEG, opt)
% ieeglab_car() - Adjusted common average referencing for iEEG, with per-
%                 stimulation-site channel selection and stimulated-contact
%                 exclusion.
%
% Usage:
%   EEG        = ieeglab_car(EEG)            % options read from EEG.ieeglab.opt
%   [EEG, out] = ieeglab_car(EEG, opt)       % options passed explicitly (headless)
%
% Options (fields of opt, or of EEG.ieeglab.opt when opt is omitted):
%   .car_method   'carla'     - CAR by Least Anticorrelation (Huang et al. 2024). DEFAULT.
%                               Adaptive: the number of reference channels is chosen
%                               from the data.
%                 'varsubset' - fixed-percentile lowest-covariance subset within
%                               64-channel blocks (Valencia et al. 2023 / HAPwave).
%                               Non-adaptive; kept for backwards compatibility.
%                 'car'       - plain common average over all good channels.
%                 'none'      - no re-referencing.
%   .car_timewin  [1x2] response window in MILLISECONDS used for ranking.
%                 Default [10 300] to match the CARLA reference implementation.
%   .car_fraction percentile (0-1) of channels kept, 'varsubset' only. Default 0.25.
%   .car_nboot    bootstrap samples for CARLA. Default 100.
%   .car_sens     logical, use CARLA's sensitive cutoff. Default false.
%   .car_linefreq line-noise fundamental, 50 or 60 Hz. Default 60.
%   .bad_channels indices or labels of channels to keep out of the reference.
%   .car_persite  logical, run the selection independently per stimulation site.
%                 Default true. Set false to pool all trials (not recommended
%                 for CCEP data - responsive channels are site-specific).
%
% Outputs:
%   EEG - re-referenced. EEG.ref and EEG.ieeglab.car record what was done.
%   out - struct array, one entry per stimulation-site group, with fields
%         group, trials, car_channels, excluded_channels, n_sel, stats
%
% Notes on the stimulated contacts: for CCEP data the stimulated pair MUST be
% excluded from the common average, otherwise the stimulation artifact is
% subtracted into every channel. The pair is resolved per epoch from, in order
% of preference, EEG.epoch(i).eventtype, the BIDS events table in
% opt.events, and EEG.event. Event types of the form 'RA1-RA2' are parsed as
% a stimulated pair.
%
% Cedric Cannard, iEEGLAB, 2026. Method attribution in ieeglab_carla.m.

if nargin < 2 || isempty(opt)
    if isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'opt') && isstruct(EEG.ieeglab.opt)
        opt = EEG.ieeglab.opt;
    else
        opt = struct();
    end
end

% Pre-1.0 field names, honoured only when the current name is absent, so an
% alias stored on the dataset by an old dialog never overrides an explicit
% request. (apply_acar is not read here: whether to re-reference at all is the
% caller's decision - ieeglab_preprocess gates this function on apply_car.)
if isfield(opt,'acar_timewin') && ~isempty(opt.acar_timewin) && ~(isfield(opt,'car_timewin') && ~isempty(opt.car_timewin))
    opt.car_timewin = opt.acar_timewin;
end
if isfield(opt,'acar_fraction') && ~isempty(opt.acar_fraction) && ~(isfield(opt,'car_fraction') && ~isempty(opt.car_fraction))
    opt.car_fraction = opt.acar_fraction;
end

% ---------------- defaults ----------------
def = struct('car_method','carla', 'car_timewin',[10 300], 'car_fraction',0.25, ...
             'car_nboot',100, 'car_sens',false, 'car_linefreq',60, ...
             'car_persite',true, 'bad_channels',[], 'verbose',true);
fn = fieldnames(def);
for ii = 1:numel(fn)
    if ~isfield(opt, fn{ii}) || isempty(opt.(fn{ii}))
        opt.(fn{ii}) = def.(fn{ii});
    end
end

out = struct('group',{},'trials',{},'car_channels',{},'excluded_channels',{},'n_sel',{},'stats',{});
if strcmpi(opt.car_method,'none')
    return
end

% CARLA is defined for CCEP data only. It picks the reference by finding the
% channels least anticorrelated with the response to a KNOWN stimulated pair,
% so on a dataset with no stimulation site there is nothing for it to rank
% against and the result is not the published method. Fall back rather than
% produce a number with no interpretation.
if strcmpi(opt.car_method,'carla')
    dmode = ieeglab_detect_mode(EEG);
    if ~strcmp(dmode, 'ccep')
        warning('ieeglab_car:carlaNotCCEP', ...
            ['CARLA applies to CCEP (single-pulse stimulation) data, but this dataset ' ...
             'looks like "%s". Falling back to a plain common average reference.\n' ...
             'Pass car_method=''car'' explicitly to silence this, or ''carla'' on ' ...
             'CCEP data whose event types name the stimulated pair (e.g. ''RA1-RA2'').'], dmode);
        opt.car_method = 'car';
    end
end

% ---------------- data ----------------
X = EEG.data;
if ismatrix(X), X = reshape(X, size(X,1), size(X,2), 1); end
X = double(X);
[C, T, N] = size(X);

if ~isfield(EEG,'times') || numel(EEG.times) ~= T
    error('ieeglab_car:noTimes', ...
        'EEG.times has %d entries but the data has %d time points. Epoch the data first.', ...
        numel(EEG.times), T);
end
tt_s = double(EEG.times(:))' / 1000;              % CARLA works in seconds
win_s = double(opt.car_timewin(:))' / 1000;

labels = string({EEG.chanlocs.labels});

% ---------------- bad channels ----------------
badIdx = local_resolve_channels(opt.bad_channels, labels);
% Channels marked bad (e.g. from channels.tsv) but not removed stay out of the reference
if isfield(EEG.chanlocs, 'status')
    isBad = cellfun(@(x) ~isempty(x) && strcmpi(char(x), 'bad'), {EEG.chanlocs.status});
    badIdx = unique([badIdx(:); find(isBad(:))]);
end

% ---------------- stimulated contacts per epoch ----------------
% Shared parser: the same definition of a stimulation site as N1, CRP and
% the connectivity matrix, read from EEG.epoch (epoch-aligned by construction).
[siteName, stimExcl] = ieeglab_epoch_sites(EEG);
nResolved = sum(~cellfun(@isempty, stimExcl));
if opt.verbose
    if nResolved == 0
        fprintf('[CAR] No stimulation sites resolved - treating data as non-CCEP (single group, no contact exclusion).\n');
    else
        fprintf('[CAR] Stimulated contacts resolved for %d/%d epochs across %d site(s).\n', ...
            nResolved, N, numel(unique(siteName(siteName~=""))));
    end
end

% ---------------- group epochs ----------------
% Per-site grouping only when stimulation sites were resolved. On non-CCEP
% data ieeglab_epoch_sites returns the condition name as the "site"; grouping
% by it would choose a different reference per condition and confound every
% condition contrast.
persite = opt.car_persite && nResolved > 0;
if persite
    [uSites, ~, grpIdx] = unique(siteName);
else
    uSites = "all"; grpIdx = ones(1,N);
end
nGrp = numel(uSites);
if ~persite && nResolved > 0 && opt.verbose
    fprintf('[CAR] Pooled reference (car_persite=false): each trial''s stimulated contacts are left out of that trial''s average.\n');
end

% ---------------- diagnostics before ----------------
r_before = local_rank_proxy(X);

% ---------------- re-reference, group by group ----------------
for g = 1:nGrp
    tr = find(grpIdx == g);
    if isempty(tr), continue; end

    % Channels barred from this group's reference selection: bad channels, and
    % for a site group its stimulated pair. A pooled group excludes stimulated
    % contacts trial by trial when the reference is applied (excluding the union
    % of every site's pair would leave no channel on a typical montage).
    stimG = [];
    if persite, stimG = vertcat(stimExcl{tr}); end
    excl = unique([badIdx(:); stimG(:)]);
    excl = excl(excl >= 1 & excl <= C);

    if strcmpi(opt.car_method,'carla') && numel(tr) > 1 && numel(tr) < 5
        warning('ieeglab_car:fewTrials', ...
            'Site "%s" has only %d trials; CARLA''s bootstrap estimate of the optimal reference size will be unstable.', ...
            uSites(g), numel(tr));
    end

    switch lower(opt.car_method)
        case 'carla'
            Vg = X(:,:,tr);
            Vg(excl,:,:) = NaN;                       % NaN = excluded from the average
            copts = struct('winResp', win_s, 'nboot', opt.car_nboot, ...
                           'sens', opt.car_sens, 'lineFreq', opt.car_linefreq, ...
                           'verbose', false);
            [~, ~, st] = ieeglab_carla(tt_s, Vg, EEG.srate, copts);
            carCh = st.chsUsed(:);
            stg   = st;
            blocks = {struct('rows', 1:C, 'sel', carCh')};

        case 'varsubset'
            [carCh, stg] = local_varsubset(X(:,:,tr), tt_s, win_s, excl, C, opt.car_fraction);
            blocks = cellfun(@(b) struct('rows', b.block, 'sel', b.selected), stg.blocks, 'UniformOutput', false);

        case 'car'
            carCh = setdiff((1:C)', excl);
            stg   = struct();
            blocks = {struct('rows', 1:C, 'sel', carCh')};

        otherwise
            error('ieeglab_car:badMethod', ...
                'Unknown car_method "%s". Use carla, varsubset, car or none.', opt.car_method);
    end

    if isempty(carCh)
        warning('ieeglab_car:noRefChannels', ...
            'Site "%s": no channels left for the reference after exclusions; skipping this group.', uSites(g));
        continue
    end

    % Apply per trial, so the common average tracks trial-level noise. Each
    % block (one for CARLA/CAR, one per 64 channels for varsubset, as in the
    % legacy HAPwave code) is referenced to the mean of its own selection.
    nNoRef = 0;
    for b = 1:numel(blocks)
        rows = blocks{b}.rows; sel = blocks{b}.sel;
        if persite || nResolved == 0
            X(rows,:,tr) = X(rows,:,tr) - mean(X(sel,:,tr), 1, 'omitnan');
        else
            for k = tr(:)'
                selk = setdiff(sel, stimExcl{k});
                if isempty(selk), nNoRef = nNoRef + 1; continue; end
                X(rows,:,k) = X(rows,:,k) - mean(X(selk,:,k), 1, 'omitnan');
            end
        end
    end
    if nNoRef > 0
        warning('ieeglab_car:trialsNotReferenced', ...
            '%d trial-block(s) had no reference channel left once their stimulated contacts were excluded, and were not re-referenced.', nNoRef);
    end

    out(end+1) = struct('group', char(uSites(g)), 'trials', tr, ...
        'car_channels', carCh(:)', 'excluded_channels', excl(:)', ...
        'n_sel', numel(carCh), 'stats', stg); %#ok<AGROW>

    if opt.verbose
        fprintf('[CAR] site %-14s %3d trial(s), %3d/%d channels in reference, %d excluded (%s)\n', ...
            char(uSites(g)), numel(tr), numel(carCh), C, numel(excl), ...
            strjoin(cellstr(labels(excl)), ','));
    end
end

% ---------------- write back ----------------
% EEG.ref names a reference only when one was actually applied.
if isempty(out)
    warning('ieeglab_car:nothingApplied', ...
        'No group could be re-referenced; the data are unchanged and EEG.ref is left as "%s".', local_refname(EEG));
    return
end
EEG.data = X;
switch lower(opt.car_method)
    case 'carla',     EEG.ref = 'CARLA';
    case 'varsubset', EEG.ref = 'CAR (variance subset)';
    case 'car',       EEG.ref = 'CAR';
end
EEG.ieeglab.car = struct('method', opt.car_method, 'timewin_ms', opt.car_timewin, ...
    'persite', persite, 'n_groups', numel(out));

r_after = local_rank_proxy(X);
if r_after < r_before
    warning('ieeglab_car:rankDrop', ...
        'Effective rank proxy decreased %d -> %d. This is expected for common average referencing.', ...
        r_before, r_after);
end

% ---------------- citation ----------------
if opt.verbose
    fprintf('\nPlease cite the following for the referencing performed here:\n');
    switch lower(opt.car_method)
        case 'carla'
            fprintf(['  Huang H., et al. (2024). CARLA: Adjusted common average referencing for\n' ...
                     '  cortico-cortical evoked potential data. J Neurosci Methods, 407:110153.\n' ...
                     '  https://doi.org/10.1016/j.jneumeth.2024.110153\n']);
        case 'varsubset'
            fprintf(['  Ojeda Valencia G., et al. (2023). Signatures of electrical stimulation driven\n' ...
                     '  network interactions in the human limbic system. J Neurosci, 43(39):6697-6711.\n' ...
                     '  https://doi.org/10.1523/JNEUROSCI.2201-22.2023\n']);
        otherwise
            fprintf('  (plain common average reference - no method citation required)\n');
    end
end

end

% ========================== local helpers ==========================

function idx = local_resolve_channels(sel, labels)
% Accept indices, logicals, a char/string, or a cellstr of labels.
idx = [];
if isempty(sel), return; end
if islogical(sel)
    idx = find(sel(:));
elseif isnumeric(sel)
    idx = sel(:);
else
    if ischar(sel), sel = {sel}; end
    sel = string(sel(:));
    [tf, loc] = ismember(upper(strtrim(sel)), upper(labels));
    idx = loc(tf);
    if any(~tf)
        warning('ieeglab_car:unknownChannel', ...
            'These channel names were not found and are ignored: %s', ...
            strjoin(cellstr(sel(~tf)), ', '));
    end
end
idx = unique(idx(idx > 0));
end

function [carCh, st] = local_varsubset(X, tt_s, win_s, excl, C, frac)
% Legacy behaviour: within each block of 64 channels, take the channels below
% the frac-th percentile of cross-trial covariance. Kept so previously analysed
% datasets can be reproduced. Adapted from apply_ieeg_car / ccep_CAR64blocks_percent
% (Multimodal Neuroimaging Lab, Mayo Clinic).
tmask = tt_s > win_s(1) & tt_s < win_s(2);      % open interval, as apply_ieeg_car
if ~any(tmask), tmask = true(1, size(X,2)); end
carCh = [];
st = struct('blocks', {{}});
for b = 1:ceil(C/64)
    blk = ((b-1)*64+1) : min(b*64, C);
    pool = setdiff(blk(:), excl);
    if isempty(pool), continue; end
    Xi = reshape(X(pool, tmask, :), numel(pool), []);
    v  = local_nanvar(Xi);
    th = quantile(v, frac);
    sel = pool(v <= th);
    if isempty(sel), sel = pool; end
    carCh = [carCh; sel(:)]; %#ok<AGROW>
    st.blocks{end+1} = struct('block', blk, 'selected', sel(:)', 'threshold', th);
end
carCh = unique(carCh);
end

function s = local_refname(EEG)
s = '';
if isfield(EEG,'ref') && (ischar(EEG.ref) || isstring(EEG.ref)), s = char(EEG.ref); end
end

function r = local_rank_proxy(X3)
X2 = reshape(X3, size(X3,1), []);
X2 = X2 - mean(X2, 2, 'omitnan');
X2(~isfinite(X2)) = 0;
e = eig(cov(X2'));
r = sum(e > max(e)*1e-7);
end

function v = local_nanvar(X)
% Row variance ignoring NaN, without var(..., 'omitnan'): that option errors on
% some MATLAB installs (issue #9, "Invalid option.
% Option must be 'omitnan' or 'includenan'"), likely from a toolbox var.m
% shadowing MATLAB's. Same result as var(X, 0, 2, 'omitnan').
X = double(X);
ok = ~isnan(X);
n = sum(ok, 2);
X(~ok) = 0;
m = sum(X, 2) ./ n;
D = (X - m) .* ok;
v = sum(D.^2, 2) ./ (n - 1);
v(n < 2) = NaN;
end
