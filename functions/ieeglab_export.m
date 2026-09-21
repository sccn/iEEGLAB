function [files, com] = ieeglab_export(EEG, outdir, opt)
% ieeglab_export() - Save iEEGLAB results as plain files anyone can open.
%
% Usage:
%   files = ieeglab_export(EEG)              % asks for a folder
%   files = ieeglab_export(EEG, outdir)
%   files = ieeglab_export(EEG, outdir, opt)
%
% Writes, for whatever results exist on the dataset:
%   <prefix>_desc-channels_ieeglab.tsv       one row per contact: label, x, y, z,
%                                            status, status_description,
%                                            clinical_zone, in_degree, n_sites_tested
%   <prefix>_desc-n1_ieeglab.tsv             N1 results, one row per site x contact
%                                            (every tested pair, with p_adj / significant)
%   <prefix>_desc-crp_ieeglab.tsv            CRP results, one row per site x contact
%   <prefix>_desc-ccepresponse_ieeglab.tsv   connectivity matrix, rows = stimulation
%                                            sites, columns = contacts: 1 response,
%                                            0 tested without response, n/a not tested
%   <prefix>_desc-ccepamplitude_ieeglab.tsv  N1 amplitude / latency (or CRP tau_R /
%   <prefix>_desc-cceplatency_ieeglab.tsv    explained variance) of the SIGNIFICANT
%                                            responses; n/a elsewhere. The values of
%                                            every tested pair are in the n1/crp TSV.
%   <prefix>_desc-pipeline_ieeglab.json      provenance: every option used, plugin
%                                            version, reference, bad channels
%   <prefix>_ieeglab.mat                     all of the above as MATLAB variables
%
% Tables are tab-separated with a header row, and missing values are written
% as n/a - the BIDS derivatives convention - so they open directly in Excel,
% R, pandas or MATLAB. A connectivity matrix that no longer matches the N1/CRP
% results next to it (they were recomputed after it was built) is not
% exported; a warning says so.
%
% Options:
%   .formats    cellstr, any of {'tsv','json','mat'}. Default all three.
%   .prefix     filename stem. Default derived from EEG.filename, e.g.
%               sub-02_ses-01_task-ccep_run-01
%   .overwrite  default true. When false, every target is checked before
%               anything is written, so a refusal leaves the folder untouched.
%   .verbose    default true
%
% Cedric Cannard, iEEGLAB, 2026

files = {};
com = '';
if nargin < 3 || isempty(opt), opt = struct(); end
def = struct('formats', {{'tsv','json','mat'}}, 'prefix', '', 'overwrite', true, 'verbose', true);
f = fieldnames(def);
for i = 1:numel(f)
    if ~isfield(opt,f{i}) || isempty(opt.(f{i})), opt.(f{i}) = def.(f{i}); end
end
opt.formats = cellstr(string(opt.formats));
fmt = lower(strtrim(opt.formats));
unknown = setdiff(fmt, {'tsv','json','mat'});
if ~isempty(unknown)
    error('ieeglab_export:badFormat', 'Unknown format(s): %s. Use tsv, json and/or mat.', strjoin(unknown, ', '));
end

if nargin < 2 || isempty(outdir)
    start = pwd;
    if isfield(EEG,'filepath') && ~isempty(EEG.filepath) && isfolder(EEG.filepath), start = EEG.filepath; end
    outdir = uigetdir(start, 'Choose a folder for the iEEGLAB results');
    if isequal(outdir, 0), return; end
end
outdir = char(outdir);

prefix = char(opt.prefix);
if isempty(prefix), prefix = local_prefix(EEG); end
fn = @(desc, ext) fullfile(outdir, sprintf('%s_desc-%s_ieeglab.%s', prefix, desc, ext));

% ---------- collect ----------
R = struct();
haveM = isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'ccep_matrix') && ~isempty(EEG.ieeglab.ccep_matrix);
if haveM
    [ok, why] = ieeglab_matrix_current(EEG);
    if ~ok
        warning('ieeglab_export:staleMatrix', 'The connectivity matrix is not exported: %s.', why);
        haveM = false;
    end
end
R.channels = local_channel_table(EEG, haveM);
if isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'n1') && isfield(EEG.ieeglab.n1,'table') && ~isempty(EEG.ieeglab.n1.table)
    R.n1 = EEG.ieeglab.n1.table;
end
if isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'stats') && isfield(EEG.ieeglab.stats,'table') && ~isempty(EEG.ieeglab.stats.table)
    R.crp = EEG.ieeglab.stats.table;
end
if haveM
    R.ccep_matrix = EEG.ieeglab.ccep_matrix;
end
R.pipeline = local_provenance(EEG, haveM);

% ---------- plan every file first ----------
plan = struct('file', {}, 'kind', {}, 'data', {});
if any(strcmp(fmt,'tsv'))
    plan(end+1) = struct('file', fn('channels','tsv'), 'kind', 'table', 'data', {R.channels});
    if isfield(R,'n1'),  plan(end+1) = struct('file', fn('n1','tsv'),  'kind', 'table', 'data', {R.n1}); end
    if isfield(R,'crp'), plan(end+1) = struct('file', fn('crp','tsv'), 'kind', 'table', 'data', {R.crp}); end
    if haveM
        M = R.ccep_matrix;
        descs = struct('response','ccepresponse', 'amplitude_uv','ccepamplitude', 'latency_ms','cceplatency', ...
                       'tR_ms','cceptr', 'explained_var','ccepexplainedvar');
        for m = fieldnames(descs)'
            if ~isfield(M, m{1}), continue; end
            A = M.(m{1});
            if ~strcmp(m{1}, 'response'), A(M.response ~= 1) = NaN; end   % significant responses only
            plan(end+1) = struct('file', fn(descs.(m{1}),'tsv'), 'kind', 'matrix', ...
                'data', {struct('A', A, 'rows', {M.sites}, 'cols', {M.channels})}); %#ok<AGROW>
        end
    end
end
if any(strcmp(fmt,'json'))
    plan(end+1) = struct('file', fn('pipeline','json'), 'kind', 'json', 'data', {R.pipeline});
end
if any(strcmp(fmt,'mat'))
    plan(end+1) = struct('file', fullfile(outdir, sprintf('%s_ieeglab.mat', prefix)), 'kind', 'mat', 'data', {R});
end

if ~opt.overwrite
    exists = arrayfun(@(p) exist(p.file, 'file') == 2, plan);
    if any(exists)
        [~, n, e] = cellfun(@fileparts, {plan(exists).file}, 'UniformOutput', false);
        error('ieeglab_export:exists', ...
            'overwrite is false and %d target file(s) already exist in %s (nothing was written): %s', ...
            nnz(exists), outdir, strjoin(strcat(n, e), ', '));
    end
end
if ~isfolder(outdir)
    [ok, msg] = mkdir(outdir);
    if ~ok, error('ieeglab_export:mkdir', 'Cannot create %s: %s', outdir, msg); end
end

% ---------- write ----------
for k = 1:numel(plan)
    switch plan(k).kind
        case 'table',  local_write_table(plan(k).data, plan(k).file);
        case 'matrix', local_write_matrix(plan(k).data, plan(k).file);
        case 'json',   local_write_json(plan(k).data, plan(k).file);
        case 'mat',    S = plan(k).data; save(plan(k).file, '-struct', 'S', '-v7');
    end
    files{end+1} = plan(k).file; %#ok<AGROW>
end

if opt.verbose
    fprintf('[export] Wrote %d file(s) to %s\n', numel(files), outdir);
    for i = 1:numel(files)
        [~, n, e] = fileparts(files{i}); fprintf('           %s%s\n', n, e);
    end
end
com = sprintf('ieeglab_export(EEG, ''%s'', %s);', strrep(outdir, '''', ''''''), ...
    ieeglab_literal(opt, {'verbose'}));
end

% ======================= helpers =======================

function p = local_prefix(EEG)
p = '';
if isfield(EEG,'filename') && ~isempty(EEG.filename)
    [~, p] = fileparts(EEG.filename);
    p = regexprep(p, '_(ieeg|eeg)$', '');
end
if isempty(p) && isfield(EEG,'setname') && ~isempty(EEG.setname)
    p = regexprep(EEG.setname, '[^\w-]', '_');
end
if isempty(p), p = 'ieeglab'; end
end

function T = local_channel_table(EEG, haveM)
C = EEG.nbchan;
cl = EEG.chanlocs;
label = strings(C,1); x = nan(C,1); y = nan(C,1); z = nan(C,1);
status = repmat("n/a", C, 1); desc = strings(C,1); zone = strings(C,1);
for i = 1:C
    label(i) = string(cl(i).labels);
    x(i) = local_num(local_field(cl(i),'X')); y(i) = local_num(local_field(cl(i),'Y')); z(i) = local_num(local_field(cl(i),'Z'));
    status(i) = local_text(local_field(cl(i),'status'), "n/a");
    desc(i)   = local_text(local_field(cl(i),'status_description'), "");
    zone(i)   = local_text(local_field(cl(i),'clinical_zone'), "");
end
T = table(label, x, y, z, status, desc, zone, ...
    'VariableNames', {'name','x','y','z','status','status_description','clinical_zone'});
if haveM
    M = EEG.ieeglab.ccep_matrix;
    [tf, loc] = ismember(upper(cellstr(label)), upper(M.channels));
    inDeg = nan(C,1); nTested = nan(C,1);
    inDeg(tf) = M.in_degree(loc(tf));
    nt = sum(~isnan(M.response), 1)';
    nTested(tf) = nt(loc(tf));
    T.in_degree = inDeg;
    T.n_sites_tested = nTested;
end
end

function v = local_field(s, fld)
if isfield(s, fld), v = s.(fld); else, v = []; end
end

function v = local_num(x)
if isempty(x) || ~isnumeric(x), v = NaN; else, v = double(x(1)); end
end

function s = local_text(x, default)
if isempty(x) || ~(ischar(x) || isstring(x)), s = default; return; end
s = string(x); s = s(1);
if ismissing(s) || strtrim(s) == "", s = default; end
end

function local_write_table(T, f)
fid = fopen(f, 'w');
if fid < 0, error('ieeglab_export:open', 'Cannot write %s', f); end
c = onCleanup(@() fclose(fid));
vn = T.Properties.VariableNames;
fprintf(fid, '%s\n', strjoin(vn, '\t'));
cols = cell(1, numel(vn));
for j = 1:numel(vn), cols{j} = local_col2str(T.(vn{j})); end
for i = 1:height(T)
    row = strings(1, numel(vn));
    for j = 1:numel(vn), row(j) = cols{j}(i); end
    fprintf(fid, '%s\n', strjoin(row, sprintf('\t')));
end
end

function local_write_matrix(D, f)
fid = fopen(f, 'w');
if fid < 0, error('ieeglab_export:open', 'Cannot write %s', f); end
c = onCleanup(@() fclose(fid));
fprintf(fid, 'stimulation_site\t%s\n', strjoin(cellstr(D.cols(:)'), '\t'));
for i = 1:numel(D.rows)
    vals = local_col2str(D.A(i,:)');
    fprintf(fid, '%s\t%s\n', D.rows{i}, strjoin(vals', sprintf('\t')));
end
end

function s = local_col2str(x)
% Column -> string column, BIDS n/a for missing, no embedded tabs/newlines.
n = size(x, 1);
if islogical(x)
    s = strings(n,1); s(x) = "true"; s(~x) = "false";
elseif isnumeric(x)
    x = double(x(:,1));
    s = strings(n,1);
    ok = isfinite(x);
    s(ok) = compose('%.6g', x(ok));
    s(~ok) = "n/a";
else
    s = string(x);
    s = s(:,1);
    s(ismissing(s) | strtrim(s) == "") = "n/a";
    s = regexprep(s, '[\t\r\n]+', ' ');
end
end

function local_write_json(S, f)
try
    txt = jsonencode(S, 'PrettyPrint', true);
catch
    txt = jsonencode(S);
end
fid = fopen(f, 'w');
if fid < 0, error('ieeglab_export:open', 'Cannot write %s', f); end
c = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', txt);
end

function P = local_provenance(EEG, haveM)
P = struct();
P.generated_by = 'iEEGLAB';
P.ieeglab_version = local_version();
P.date = char(datetime('now', 'Format', 'yyyy-MM-dd''T''HH:mm:ss'));
P.matlab_version = version;
P.dataset = struct('filename', local_str(EEG, 'filename'), 'setname', local_str(EEG, 'setname'), ...
    'srate', EEG.srate, 'nbchan', EEG.nbchan, 'trials', EEG.trials, 'pnts', EEG.pnts);
if isfield(EEG,'times') && ~isempty(EEG.times), P.dataset.epoch_ms = [EEG.times(1) EEG.times(end)]; end
P.reference = local_str(EEG, 'ref');
if isfield(EEG,'ieeglab')
    I = EEG.ieeglab;
    if isfield(I,'opt'),  P.options = local_sanitize(I.opt, 0); end
    if isfield(I,'car'),  P.rereferencing = local_sanitize(I.car, 0); end
    if isfield(I,'blank'), P.stimulation_blanking = local_sanitize(I.blank, 0); end
    if isfield(I,'removed_channels'), P.removed_channels = I.removed_channels; end
    if isfield(I,'n1') && isfield(I.n1,'opt'), P.n1_options = local_sanitize(I.n1.opt, 0); end
    if isfield(I,'stats') && isfield(I.stats,'opt'), P.crp_options = local_sanitize(I.stats.opt, 0); end
    if haveM
        M = I.ccep_matrix;
        P.ccep_matrix = struct('source', M.source, 'n_sites', numel(M.sites), ...
            'n_channels', numel(M.channels), 'n_tested', M.n_tested, ...
            'n_significant', M.n_significant, 'density', M.density, ...
            'metric_matrices', 'significant responses only; n/a elsewhere');
    end
end
end

function v = local_version()
v = 'unknown';
w = which('eegplugin_ieeglab');
if isempty(w), return; end
t = regexp(fileread(w), 'vers\s*=\s*''([^'']+)''', 'tokens', 'once');
if ~isempty(t), v = t{1}; end
end

function s = local_str(S, f)
s = '';
if isfield(S, f) && ~isempty(S.(f)) && (ischar(S.(f)) || isstring(S.(f))), s = char(S.(f)); end
end

function out = local_sanitize(x, depth)
% Keep what JSON can represent sensibly; summarise the rest.
if depth > 4, out = '<nested>'; return; end
if istable(x)
    out = sprintf('<table %dx%d>', height(x), width(x));
elseif isstruct(x)
    if numel(x) > 1
        out = sprintf('<struct array 1x%d>', numel(x));
        return
    end
    out = struct();
    fns = fieldnames(x);
    for i = 1:numel(fns)
        out.(fns{i}) = local_sanitize(x.(fns{i}), depth + 1);
    end
elseif (isnumeric(x) || islogical(x)) && numel(x) > 64
    out = sprintf('<%s %s>', class(x), mat2str(size(x)));
elseif iscell(x)
    if all(cellfun(@(c) ischar(c) || isstring(c), x(:))) && numel(x) <= 256
        out = cellstr(x(:)');
    else
        out = sprintf('<cell %s>', mat2str(size(x)));
    end
elseif isa(x, 'function_handle')
    out = func2str(x);
else
    out = x;
end
end
