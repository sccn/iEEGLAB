function [EEG, T, com] = ieeglab_bad_channels(EEG, opt)
% ieeglab_bad_channels() - Mark or remove bad channels, primarily from the
%                          clinician / annotator labels that ship with BIDS data.
%
% Usage:
%   EEG         = ieeglab_bad_channels(EEG)            % options from EEG.ieeglab.opt
%   [EEG, T]    = ieeglab_bad_channels(EEG, opt)       % headless
%
% Sources, combined (each optional):
%   1. BIDS *_channels.tsv          opt.channels_tsv  (path, or 'auto' to find the
%      sibling of the dataset file). A row with status == 'bad' marks the channel
%      bad and keeps status_description as the reason. Rows whose type is not an
%      intracranial type (e.g. ECG, EMG, TRIG, MISC, DC) are flagged 'non-iEEG'.
%   2. electrodes.tsv 'seizure_zone' opt.elec_tsv (path or table). Recorded per
%      channel as clinical_zone ('SOZ', 'IrritativeZone', ...). Only removed when
%      opt.exclude_soz / opt.exclude_irritative ask for it.
%   3. opt.bad_channels             explicit labels or indices from the user.
%   4. opt.auto_detect              flat or extreme-variance channels. Off by
%      default: the clinician labels are the primary source.
%
% Options:
%   .channels_tsv        '' | path | 'auto'          default 'auto'
%   .elec_tsv            '' | path | table           default ''
%   .bad_channels        labels or indices           default []
%   .drop_non_ieeg       treat non-iEEG types as bad default true
%   .exclude_soz         also treat SOZ contacts as bad          default false
%   .exclude_irritative  also treat irritative-zone contacts     default false
%   .auto_detect         flat / noisy detection      default false
%   .auto_z              robust z threshold on log-SD for 'noisy' default 5
%   .action              'remove' (default) | 'mark'
%                        'mark' keeps the data but records status='bad', which the
%                        re-referencing step reads and keeps out of the reference.
%   .drop_bad_stim_sites for CCEP data, also drop trials whose stimulated pair
%                        includes a bad contact. Default true, as in the original
%                        pipeline (is_good_pair). Uses exact token matching.
%   .honor_previous      keep marks from an earlier pass (e.g. ieeglab_load). Default true.
%   .keep_channels       labels forced GOOD, overriding every source (a user
%                        decision to keep a clinician-marked channel). Default [].
%   .verbose             default true
%
% Outputs:
%   EEG  - chanlocs gain status / status_description / clinical_zone fields;
%          EEG.ieeglab.bad_channels holds T; removed labels are kept in
%          EEG.ieeglab.removed_channels.
%   T    - table: label, status, reason, clinical_zone, source
%   com  - EEGLAB history string
%
% Cedric Cannard, iEEGLAB, 2026

com = '';
if nargin < 2 || isempty(opt)
    if isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'opt') && isstruct(EEG.ieeglab.opt)
        opt = EEG.ieeglab.opt;
    else
        opt = struct();
    end
end
def = struct('channels_tsv','auto', 'elec_tsv','', 'bad_channels',[], ...
    'drop_non_ieeg',true, 'exclude_soz',false, 'exclude_irritative',false, ...
    'auto_detect',false, 'auto_z',5, 'action','remove', ...
    'drop_bad_stim_sites',true, 'honor_previous',true, 'keep_channels',[], 'verbose',true);
f = fieldnames(def);
for i = 1:numel(f)
    if ~isfield(opt, f{i}) || (isempty(opt.(f{i})) && ~ischar(def.(f{i})))
        opt.(f{i}) = def.(f{i});
    end
end
if ~any(strcmpi(opt.action, {'remove','mark'}))
    error('ieeglab_bad_channels:badAction', 'action must be ''remove'' or ''mark''; got ''%s''.', opt.action);
end

C = EEG.nbchan;
labels = string({EEG.chanlocs.labels});
L = upper(strtrim(labels));
status = repmat("good", 1, C);
reason = strings(1, C);
source = strings(1, C);
zone   = strings(1, C);

% ---------- 0. marks from an earlier pass ----------
% ieeglab_load annotates with action='mark'; preprocessing then calls this again
% with action='remove'. Earlier marks are honoured so that second call acts on
% them even when it is not re-reading the same files.
if opt.honor_previous && isfield(EEG.chanlocs, 'status')
    for c = 1:C
        s0 = EEG.chanlocs(c).status;
        if ~isempty(s0) && strcmpi(string(s0), "bad")
            status(c) = "bad";
            source(c) = "previous";
            % EEGLAB's BIDS import can leave NaN or 'n/a' here (bids_loadfile reads n/a as NaN)
            d0 = "";
            if isfield(EEG.chanlocs,'status_description'), d0 = local_text(EEG.chanlocs(c).status_description); end
            if strlength(d0) > 0
                reason(c) = d0;
            else
                reason(c) = "marked bad earlier";
            end
        end
        if isfield(EEG.chanlocs,'clinical_zone') && strlength(local_text(EEG.chanlocs(c).clinical_zone)) > 0
            zone(c) = local_text(EEG.chanlocs(c).clinical_zone);
        end
    end
end

% ---------- 1. BIDS channels.tsv ----------
chTsv = local_resolve_channels_tsv(opt.channels_tsv, EEG);
if chTsv ~= ""
    Tc = readtable(chTsv, 'FileType','text', 'Delimiter','\t', 'TextType','string');
    vn = lower(Tc.Properties.VariableNames);
    iName = find(strcmp(vn,'name'), 1);
    if isempty(iName)
        warning('ieeglab_bad_channels:noNameColumn', ...
            '%s has no "name" column; ignoring it.', chTsv);
    else
        names = upper(strtrim(string(Tc{:, iName})));
        iStat = find(strcmp(vn,'status'), 1);
        iDesc = find(strcmp(vn,'status_description'), 1);
        iType = find(strcmp(vn,'type'), 1);
        nMatched = 0;
        for r = 1:height(Tc)
            c = find(L == names(r), 1);
            if isempty(c), continue; end
            nMatched = nMatched + 1;
            if ~isempty(iStat) && strcmpi(strtrim(string(Tc{r, iStat})), "bad") && source(c) ~= "previous"
                status(c) = "bad";
                source(c) = "channels.tsv";
                if ~isempty(iDesc)
                    d = strtrim(string(Tc{r, iDesc}));
                    if ~ismissing(d) && d ~= "" && lower(d) ~= "n/a", reason(c) = d; end
                end
                if reason(c) == "", reason(c) = "marked bad in channels.tsv"; end
            end
            if opt.drop_non_ieeg && ~isempty(iType)
                ty = upper(strtrim(string(Tc{r, iType})));
                if ~ismissing(ty) && ty ~= "" && ty ~= "N/A" && ~ismember(ty, ["SEEG","ECOG","DBS","EEG"]) && status(c) == "good"
                    status(c) = "bad";
                    source(c) = "channels.tsv";
                    reason(c) = "non-iEEG channel type " + ty;
                end
            end
        end
        if opt.verbose
            fprintf('[bad channels] %s: %d/%d channels matched, %d marked bad.\n', ...
                local_short(chTsv), nMatched, C, nnz(source == "channels.tsv"));
        end
    end
end

% ---------- 2. seizure_zone from electrodes.tsv ----------
Te = local_electrodes_table(opt.elec_tsv);
if ~isempty(Te)
    vn = lower(Te.Properties.VariableNames);
    iName = find(strcmp(vn,'name'), 1);
    iZone = find(strcmp(vn,'seizure_zone'), 1);
    if ~isempty(iName) && ~isempty(iZone)
        names = upper(strtrim(string(Te{:, iName})));
        zones = strtrim(string(Te{:, iZone}));
        for c = 1:C
            r = find(names == L(c), 1);
            if isempty(r), continue; end
            z = zones(r);
            if ~ismissing(z) && z ~= "" && lower(z) ~= "n/a" && zone(c) == "", zone(c) = z; end
        end
    end
end
isSOZ = contains(lower(zone), "soz") | contains(lower(zone), "seizure");
isIRR = contains(lower(zone), "irrit");
if opt.exclude_soz
    idx = isSOZ & status == "good";
    status(idx) = "bad"; source(idx) = "seizure_zone"; reason(idx) = "seizure onset zone (" + zone(idx) + ")";
end
if opt.exclude_irritative
    idx = isIRR & status == "good";
    status(idx) = "bad"; source(idx) = "seizure_zone"; reason(idx) = "irritative zone (" + zone(idx) + ")";
end

% ---------- 3. explicit list ----------
if ~isempty(opt.bad_channels)
    b = opt.bad_channels;
    if isnumeric(b) || islogical(b)
        if islogical(b), b = find(b); end
        idx = false(1, C); idx(b(b>=1 & b<=C)) = true;
    else
        if ischar(b), b = {b}; end
        want = upper(strtrim(string(b(:)')));
        idx = ismember(L, want);
        unknown = want(~ismember(want, L));
        if ~isempty(unknown)
            warning('ieeglab_bad_channels:unknownLabel', ...
                'These bad-channel labels are not in the dataset and are ignored: %s', strjoin(unknown, ', '));
        end
    end
    newly = idx & status == "good";
    status(newly) = "bad"; source(newly) = "user"; reason(newly) = "listed by user";
end

% ---------- 4. automatic detection ----------
if opt.auto_detect
    X = reshape(double(EEG.data), C, []);
    sd = std(X, 0, 2, 'omitnan')';
    med = median(sd(sd > 0), 'omitnan');
    flat = ~isfinite(sd) | sd <= 1e-3 * med;
    lsd = log(max(sd, eps));
    ok = ~flat;
    mu = median(lsd(ok)); mad1 = 1.4826 * median(abs(lsd(ok) - mu));
    z = (lsd - mu) / max(mad1, eps);
    noisy = ok & z > opt.auto_z;
    newly = flat & status == "good";
    status(newly) = "bad"; source(newly) = "auto"; reason(newly) = "flat signal";
    newly = noisy & status == "good";
    status(newly) = "bad"; source(newly) = "auto";
    reason(newly) = "extreme variance (robust z = " + compose('%.1f', z(newly)) + ")";
end

% ---------- user overrides ----------
if ~isempty(opt.keep_channels)
    kc = opt.keep_channels;
    if isnumeric(kc) || islogical(kc)
        if islogical(kc), kc = find(kc); end
        kc = cellstr(labels(kc(kc>=1 & kc<=C)));
    end
    ov = ismember(L, upper(strtrim(string(kc)))) & status == "bad";
    for k = find(ov)
        reason(k) = "mark overridden by user (was: " + reason(k) + ")";
    end
    status(ov) = "good"; source(ov) = "user override";
end

% ---------- record on chanlocs ----------
for c = 1:C
    EEG.chanlocs(c).status             = char(status(c));
    EEG.chanlocs(c).status_description = char(reason(c));
    EEG.chanlocs(c).clinical_zone      = char(zone(c));
end
T = table(labels(:), status(:), reason(:), zone(:), source(:), ...
    'VariableNames', {'label','status','reason','clinical_zone','source'});
EEG.ieeglab.bad_channels = T;
if ~isfield(EEG.ieeglab,'orig_labels') || isempty(EEG.ieeglab.orig_labels)
    EEG.ieeglab.orig_labels = cellstr(labels);
end

isBad = status == "bad";
if opt.verbose
    if any(isBad)
        fprintf('[bad channels] %d bad: %s\n', nnz(isBad), strjoin(labels(isBad) + " (" + reason(isBad) + ")", ', '));
    else
        fprintf('[bad channels] No bad channels.\n');
    end
    if any(isSOZ) && ~opt.exclude_soz
        fprintf('[bad channels] %d seizure-onset-zone contact(s) kept (exclude_soz is off): %s\n', ...
            nnz(isSOZ), strjoin(labels(isSOZ), ', '));
    end
end

% A dataset in which every channel is bad has nothing left to analyse; say so
% before any trial is dropped (dropping first gave EEGLAB's cryptic 'dataset
% is empty', even for action='mark').
if all(isBad)
    error('ieeglab_bad_channels:allBad', 'Every channel is marked bad; nothing would remain to analyse.');
end

% ---------- trials stimulating a bad contact (CCEP) ----------
badLabels = cellstr(labels(isBad));
if opt.drop_bad_stim_sites && any(isBad)
    if EEG.trials > 1
        sites = ieeglab_epoch_sites(EEG);
        dropTr = false(1, EEG.trials);
        [~, stimIdx] = ieeglab_epoch_sites(EEG);
        badIdx = find(isBad);
        for i = 1:EEG.trials
            dropTr(i) = numel(stimIdx{i}) >= 1 && sites(i) ~= "" && any(ismember(stimIdx{i}, badIdx)) ...
                && numel(ieeglab_site_tokens(sites(i), cellstr(labels))) >= 2;
        end
        if all(dropTr)
            error('ieeglab_bad_channels:allTrialsBad', ...
                'Every trial stimulates a bad contact, so none would remain.');
        end
        if any(dropTr)
            if opt.verbose
                fprintf('[bad channels] Dropping %d/%d trials whose stimulated pair includes a bad contact.\n', ...
                    nnz(dropTr), EEG.trials);
            end
            EEG = pop_select(EEG, 'notrial', find(dropTr));
        end
    elseif isfield(EEG,'event') && ~isempty(EEG.event)
        dropEv = ieeglab_events_using(EEG.event, badLabels);
        if any(dropEv)
            if opt.verbose
                fprintf('[bad channels] Dropping %d/%d stimulation events whose pair includes a bad contact.\n', ...
                    nnz(dropEv), numel(dropEv));
            end
            EEG.event(dropEv) = [];
            if isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'opt') && isfield(EEG.ieeglab.opt,'events') ...
                    && istable(EEG.ieeglab.opt.events) && height(EEG.ieeglab.opt.events) == numel(dropEv)
                EEG.ieeglab.opt.events(dropEv,:) = [];
            end
            EEG = eeg_checkset(EEG, 'eventconsistency');
        end
    end
end

% ---------- act ----------
if strcmpi(opt.action, 'remove') && any(isBad)
    keepLabels = setdiff(cellstr(labels), badLabels, 'stable');
    if isempty(keepLabels)
        error('ieeglab_bad_channels:allBad', 'Every channel is marked bad; nothing would remain.');
    end
    EEG = pop_select(EEG, 'nochannel', find(isBad));
    prev = {};
    if isfield(EEG.ieeglab,'removed_channels'), prev = EEG.ieeglab.removed_channels; end
    EEG.ieeglab.removed_channels = unique([prev(:); badLabels(:)], 'stable');
end
EEG = eeg_checkset(EEG);
used = struct('channels_tsv', char(string(opt.channels_tsv)), 'bad_channels', opt.bad_channels, ...
    'drop_non_ieeg', opt.drop_non_ieeg, 'exclude_soz', opt.exclude_soz, ...
    'exclude_irritative', opt.exclude_irritative, 'auto_detect', opt.auto_detect, 'auto_z', opt.auto_z, ...
    'action', char(opt.action), 'drop_bad_stim_sites', opt.drop_bad_stim_sites, ...
    'honor_previous', opt.honor_previous, 'keep_channels', opt.keep_channels);
if ischar(opt.elec_tsv) || isstring(opt.elec_tsv), used.elec_tsv = char(opt.elec_tsv); end
if iscell(used.bad_channels) || isstring(used.bad_channels), used.bad_channels = cellstr(string(used.bad_channels)); end
com = sprintf('EEG = ieeglab_bad_channels(EEG, %s);', ieeglab_literal(used));
end

% ======================= helpers =======================

function p = local_resolve_channels_tsv(spec, EEG)
p = "";
if isempty(spec), return; end
spec = char(spec);
if ~strcmpi(spec, 'auto')
    if exist(spec, 'file') == 2
        p = string(spec);
    else
        warning('ieeglab_bad_channels:missingFile', 'channels.tsv not found: %s', spec);
    end
    return
end
% 'auto': the BIDS sidecar of THIS dataset; never another run's
p = string(ieeglab_bids_sibling(EEG, 'channels'));
end

function Te = local_electrodes_table(spec)
Te = [];
if istable(spec), Te = spec; return; end
if isempty(spec), return; end
if exist(char(spec),'file') == 2
    Te = readtable(char(spec), 'FileType','text', 'Delimiter','\t', 'TextType','string');
end
end

function s = local_short(p)
[~, n, e] = fileparts(char(p)); s = [n e];
end

function s = local_text(x)
% A free-text field as a string; "" for empty, NaN, missing or n/a.
s = "";
if isempty(x) || (isnumeric(x) && all(isnan(x(:)))), return; end
s = strtrim(string(x));
if numel(s) ~= 1 || ismissing(s) || strcmpi(s, "n/a"), s = ""; end
end
