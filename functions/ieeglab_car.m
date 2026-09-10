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

% Back-compatibility with the previous field names
if isfield(opt,'acar_timewin')  && ~isempty(opt.acar_timewin),  opt.car_timewin  = opt.acar_timewin;  end
if isfield(opt,'acar_fraction') && ~isempty(opt.acar_fraction), opt.car_fraction = opt.acar_fraction; end
if isfield(opt,'apply_acar') && ~opt.apply_acar, opt.car_method = 'none'; end

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

% ---------------- stimulated contacts per epoch ----------------
[stimExcl, siteName] = local_stim_sites(EEG, labels, N, opt);
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
if opt.car_persite && any(siteName ~= "")
    [uSites, ~, grpIdx] = unique(siteName);
else
    uSites = "all"; grpIdx = ones(1,N);
end
nGrp = numel(uSites);

% ---------------- diagnostics before ----------------
r_before = local_rank_proxy(X);

% ---------------- re-reference, group by group ----------------
for g = 1:nGrp
    tr = find(grpIdx == g);
    if isempty(tr), continue; end

    % Channels barred from this group's reference: stimulated pair + bad channels
    excl = unique([badIdx(:); vertcat(stimExcl{tr})]);
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
            carCh = st.chsUsed;
            stg   = st;

        case 'varsubset'
            [carCh, stg] = local_varsubset(X(:,:,tr), tt_s, win_s, excl, C, opt.car_fraction);

        case 'car'
            carCh = setdiff((1:C)', excl);
            stg   = struct();

        otherwise
            error('ieeglab_car:badMethod', ...
                'Unknown car_method "%s". Use carla, varsubset, car or none.', opt.car_method);
    end

    if isempty(carCh)
        warning('ieeglab_car:noRefChannels', ...
            'Site "%s": no channels left for the reference after exclusions; skipping this group.', uSites(g));
        continue
    end

    % Apply per trial, so the common average tracks trial-level noise
    mu = mean(X(carCh,:,tr), 1, 'omitnan');           % [1 x T x numel(tr)]
    X(:,:,tr) = X(:,:,tr) - mu;

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
EEG.data = X;
switch lower(opt.car_method)
    case 'carla',     EEG.ref = 'CARLA';
    case 'varsubset', EEG.ref = 'CAR (variance subset)';
    case 'car',       EEG.ref = 'CAR';
end
EEG.ieeglab.car = struct('method', opt.car_method, 'timewin_ms', opt.car_timewin, ...
    'persite', opt.car_persite, 'n_groups', numel(out));

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

function [stimExcl, siteName] = local_stim_sites(EEG, labels, N, opt)
% Resolve the stimulated contact indices for every epoch.
%
% The previous implementation relied on a BIDS events table whose row count had
% to equal the epoch count. pop_epoch silently drops boundary epochs, so that
% test failed and exclusion became a silent no-op. We now read EEG.epoch first,
% which pop_epoch always populates and which is aligned to the data by
% construction.
stimExcl = repmat({zeros(0,1)}, 1, N);
siteName = strings(1, N);

% --- source 1: EEG.epoch (always epoch-aligned) ---
if isfield(EEG,'epoch') && ~isempty(EEG.epoch) && numel(EEG.epoch) == N
    for i = 1:N
        ty = local_epoch_type(EEG.epoch(i));
        if ty == "", continue; end
        stimExcl{i} = local_tokens_to_idx(local_split_labels(ty), labels);
        % Group by the canonical (order-independent) pair so that 'ROP2-ROP4'
        % and 'ROP4-ROP2' are recognised as the same stimulation site.
        siteName(i) = local_canonical_site(ty, stimExcl{i}, labels);
    end
end

% --- source 2: BIDS events table, when it happens to be epoch-aligned ---
if all(cellfun(@isempty, stimExcl)) && isfield(opt,'events') && istable(opt.events) ...
        && height(opt.events) == N
    ev = opt.events;
    cand = {'electrical_stimulation_site','electrodes_involved_onset','stim_electrodes'};
    for i = 1:N
        for f = 1:numel(cand)
            if ~ismember(cand{f}, ev.Properties.VariableNames), continue; end
            tok = local_split_labels(ev{i, cand{f}});
            idx = local_tokens_to_idx(tok, labels);
            if ~isempty(idx)
                stimExcl{i} = idx;
                siteName(i) = string(strjoin(cellstr(labels(idx)), '-'));
                break
            end
        end
    end
end

% --- source 3: continuous data with events but no epochs ---
if all(cellfun(@isempty, stimExcl)) && N == 1 && isfield(EEG,'event') && ~isempty(EEG.event)
    types = unique(string({EEG.event.type}));
    idx = local_tokens_to_idx(local_split_labels(strjoin(cellstr(types),' ')), labels);
    stimExcl{1} = idx;
end
end

function site = local_canonical_site(ty, idx, labels)
% Order-independent name for a stimulation site, so that a pair recorded as
% 'ROP2-ROP4' in some trials and 'ROP4-ROP2' in others forms one group.
if isempty(idx)
    site = ty;
else
    site = string(strjoin(sort(cellstr(labels(idx))), '-'));
end
end

function ty = local_epoch_type(ep)
% The event type at latency 0 for this epoch. pop_epoch stores eventtype as a
% cell when several events fall inside the epoch window.
ty = "";
if ~isfield(ep,'eventtype') || isempty(ep.eventtype), return; end
t = ep.eventtype;
l = [];
if isfield(ep,'eventlatency') && ~isempty(ep.eventlatency), l = ep.eventlatency; end
if iscell(t)
    if iscell(l) && numel(l) == numel(t)
        lv = cellfun(@(x) abs(double(x(1))), l);
        [~, k] = min(lv);                 % the event the epoch is locked to
    else
        k = 1;
    end
    t = t{k};
end
if isnumeric(t), t = num2str(t); end
ty = strtrim(string(t));
end

function tokens = local_split_labels(val)
% Split a stimulation-site string such as 'RA1-RA2' into contact labels.
% BIDS missing markers are treated as absent rather than as a label, which is
% what previously let an 'n/a' column mask the real stimulation-site column.
tokens = {};
if isempty(val), return; end
if istable(val), val = val{1,1}; end
if iscell(val), val = string(val); end
if iscategorical(val), val = string(val); end
if isnumeric(val)
    if all(isnan(val(:))), return; end
    tokens = cellstr(string(val(:)')); return
end
s = char(strjoin(string(val(:))', ' '));
parts = regexp(s, '[,;+\-\/\|\s]+', 'split');
parts = parts(~cellfun(@isempty, parts));
bad = ~cellfun(@isempty, regexpi(parts, '^(n/?a|nan|none|undefined|missing|\?)$', 'once'));
tokens = parts(~bad);
end

function idx = local_tokens_to_idx(tokens, labels)
idx = zeros(0,1);
if isempty(tokens), return; end
tok = string(tokens(:));
[tf, loc] = ismember(upper(strtrim(tok)), upper(labels));
idx = loc(tf);
if isempty(idx)
    tok2 = regexprep(tok, '[^\w]', '');           % strip punctuation and retry
    [tf2, loc2] = ismember(upper(tok2), upper(regexprep(labels,'[^\w]','')));
    idx = loc2(tf2);
end
idx = unique(idx(idx > 0));
idx = idx(:);
end

function [carCh, st] = local_varsubset(X, tt_s, win_s, excl, C, frac)
% Legacy behaviour: within each block of 64 channels, take the channels below
% the frac-th percentile of cross-trial covariance. Kept so previously analysed
% datasets can be reproduced. Adapted from apply_ieeg_car / ccep_CAR64blocks_percent
% (Multimodal Neuroimaging Lab, Mayo Clinic).
tmask = tt_s >= win_s(1) & tt_s <= win_s(2);
if ~any(tmask), tmask = true(1, size(X,2)); end
carCh = [];
st = struct('blocks', {{}});
for b = 1:ceil(C/64)
    blk = ((b-1)*64+1) : min(b*64, C);
    pool = setdiff(blk(:), excl);
    if isempty(pool), continue; end
    Xi = reshape(X(pool, tmask, :), numel(pool), []);
    v  = var(Xi, 0, 2, 'omitnan');
    th = quantile(v, frac);
    sel = pool(v <= th);
    if isempty(sel), sel = pool; end
    carCh = [carCh; sel(:)]; %#ok<AGROW>
    st.blocks{end+1} = struct('block', blk, 'selected', sel(:)', 'threshold', th);
end
carCh = unique(carCh);
end

function r = local_rank_proxy(X3)
X2 = reshape(X3, size(X3,1), []);
X2 = X2 - mean(X2, 2, 'omitnan');
X2(~isfinite(X2)) = 0;
e = eig(cov(X2'));
r = sum(e > max(e)*1e-7);
end
