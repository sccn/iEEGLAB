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
%                                            clinical_zone, in_degree, out of n tested
%   <prefix>_desc-n1_ieeglab.tsv             N1 results, one row per site x contact
%   <prefix>_desc-crp_ieeglab.tsv            CRP results, one row per site x contact
%   <prefix>_desc-ccepresponse_ieeglab.tsv   connectivity matrix, rows = stimulation
%   <prefix>_desc-ccepamplitude_ieeglab.tsv    sites, columns = contacts
%   <prefix>_desc-cceplatency_ieeglab.tsv
%   <prefix>_desc-pipeline_ieeglab.json      provenance: every option used, plugin
%                                            version, reference, bad channels
%   <prefix>_ieeglab.mat                     all of the above as MATLAB variables
%
% Tables are tab-separated with a header row, and missing values are written
% as n/a - the BIDS derivatives convention - so they open directly in Excel,
% R, pandas or MATLAB.
%
% Options:
%   .formats    cellstr, any of {'tsv','json','mat'}. Default all three.
%   .prefix     filename stem. Default derived from EEG.filename, e.g.
%               sub-02_ses-01_task-ccep_run-01
%   .overwrite  default true
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
if ischar(opt.formats), opt.formats = {opt.formats}; end
fmt = lower(opt.formats);

if nargin < 2 || isempty(outdir)
    start = pwd;
    if isfield(EEG,'filepath') && ~isempty(EEG.filepath) && isfolder(EEG.filepath), start = EEG.filepath; end
    outdir = uigetdir(start, 'Choose a folder for the iEEGLAB results');
    if isequal(outdir, 0), return; end
end
if ~isfolder(outdir)
    [ok, msg] = mkdir(outdir);
    if ~ok, error('ieeglab_export:mkdir', 'Cannot create %s: %s', outdir, msg); end
end

prefix = char(opt.prefix);
if isempty(prefix), prefix = local_prefix(EEG); end
fn = @(desc, ext) fullfile(outdir, sprintf('%s_desc-%s_ieeglab.%s', prefix, desc, ext));

R = struct();

% ---------- channels ----------
R.channels = local_channel_table(EEG);

% ---------- results tables ----------
if isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'n1') && isfield(EEG.ieeglab.n1,'table') && ~isempty(EEG.ieeglab.n1.table)
    R.n1 = EEG.ieeglab.n1.table;
end
if isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'stats') && isfield(EEG.ieeglab.stats,'table') && ~isempty(EEG.ieeglab.stats.table)
    R.crp = EEG.ieeglab.stats.table;
end
if isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'ccep_matrix') && ~isempty(EEG.ieeglab.ccep_matrix)
    R.ccep_matrix = EEG.ieeglab.ccep_matrix;
end
R.pipeline = local_provenance(EEG);

% ---------- write ----------
if any(strcmp(fmt,'tsv'))
    files{end+1} = local_write_table(R.channels, fn('channels','tsv'), opt.overwrite);
    if isfield(R,'n1'),  files{end+1} = local_write_table(R.n1,  fn('n1','tsv'),  opt.overwrite); end
    if isfield(R,'crp'), files{end+1} = local_write_table(R.crp, fn('crp','tsv'), opt.overwrite); end
    if isfield(R,'ccep_matrix')
        M = R.ccep_matrix;
        for m = {'response','amplitude_uv','latency_ms','tR_ms','explained_var'}
            if isfield(M, m{1})
                desc = ['ccep' regexprep(m{1}, '_(uv|ms)$', '')];
                desc = regexprep(desc, '[^A-Za-z0-9]', '');
                files{end+1} = local_write_matrix(M.(m{1}), M.sites, M.channels, fn(desc,'tsv'), opt.overwrite); %#ok<AGROW>
            end
        end
    end
end
if any(strcmp(fmt,'json'))
    files{end+1} = local_write_json(R.pipeline, fn('pipeline','json'), opt.overwrite);
end
if any(strcmp(fmt,'mat'))
    matFile = fullfile(outdir, sprintf('%s_ieeglab.mat', prefix));
    local_check_overwrite(matFile, opt.overwrite);
    save(matFile, '-struct', 'R', '-v7');
    files{end+1} = matFile;
end
files = files(~cellfun(@isempty, files));

if opt.verbose
    fprintf('[export] Wrote %d file(s) to %s\n', numel(files), outdir);
    for i = 1:numel(files)
        [~, n, e] = fileparts(files{i}); fprintf('           %s%s\n', n, e);
    end
end
com = sprintf('ieeglab_export(EEG, ''%s'');', outdir);
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

function T = local_channel_table(EEG)
C = EEG.nbchan;
cl = EEG.chanlocs;
get = @(fld, i) local_field(cl(i), fld);
label = strings(C,1); x = nan(C,1); y = nan(C,1); z = nan(C,1);
status = repmat("good", C, 1); desc = strings(C,1); zone = strings(C,1);
for i = 1:C
    label(i) = string(cl(i).labels);
    x(i) = local_num(get('X',i)); y(i) = local_num(get('Y',i)); z(i) = local_num(get('Z',i));
    s = get('status',i); if ~isempty(s), status(i) = string(s); end
    desc(i) = string(get('status_description',i));
    zone(i) = string(get('clinical_zone',i));
end
T = table(label, x, y, z, status, desc, zone, ...
    'VariableNames', {'name','x','y','z','status','status_description','clinical_zone'});
if isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'ccep_matrix') && ~isempty(EEG.ieeglab.ccep_matrix)
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

function f = local_write_table(T, f, overwrite)
local_check_overwrite(f, overwrite);
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

function f = local_write_matrix(A, rows, cols, f, overwrite)
local_check_overwrite(f, overwrite);
fid = fopen(f, 'w');
if fid < 0, error('ieeglab_export:open', 'Cannot write %s', f); end
c = onCleanup(@() fclose(fid));
fprintf(fid, 'stimulation_site\t%s\n', strjoin(cellstr(cols(:)'), '\t'));
for i = 1:numel(rows)
    vals = local_col2str(A(i,:)');
    fprintf(fid, '%s\t%s\n', rows{i}, strjoin(vals', sprintf('\t')));
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

function f = local_write_json(S, f, overwrite)
local_check_overwrite(f, overwrite);
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

function local_check_overwrite(f, overwrite)
if ~overwrite && exist(f, 'file') == 2
    error('ieeglab_export:exists', '%s already exists and overwrite is false.', f);
end
end

function P = local_provenance(EEG)
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
    if isfield(I,'ccep_matrix') && ~isempty(I.ccep_matrix)
        M = I.ccep_matrix;
        P.ccep_matrix = struct('source', M.source, 'n_sites', numel(M.sites), ...
            'n_channels', numel(M.channels), 'n_tested', M.n_tested, ...
            'n_significant', M.n_significant, 'density', M.density);
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
