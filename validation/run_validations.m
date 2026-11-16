function run_validations(outDir, varargin)
% run_validations() - runs the method validations and writes, to one folder, the
%       files that validation/figures/figures_validations.py draws.
%
% Usage:
%   run_validations(outDir, 'fieldtrip', uci29Dir, 'crp', crpDir, ...
%       'erdetect', {ds004080Dir, erdetectDir}, 'ica', uci29Dir, ...
%       'carla', {ds004977Dir, carlaDir, mnlDir})
%
% A validation runs when its data are given; each function's help says what to
% download.
%   'fieldtrip'  validate_fieldtrip_uci29: fieldtrip_uci29_figure.mat and
%                fieldtrip_uci29_band_stats.tsv
%   'crp'        validate_crp_scripts: crp_vs_crp_scripts.tsv
%   'erdetect'   validate_erdetect_ds004080: validate_erdetect_ds004080_out.mat
%   'ica'        validate_icaref_michelmann: validate_icaref_michelmann_out.mat
%   'carla'      validate_carla_huang_ds004977: carla_fig5.json
% The re-referencing comparison on ds004696 (compare_references_ds004696, then
% matched_pipeline_ds004696 and split_half_effects_ds004696) writes its own
% output folder, which validation/figures/figures_rereferencing.py reads.
%
% Then, in a shell:
%   python figures_validations.py <outDir> <figure folder> [<SubjectUCI29 folder>]
%
% Cedric Cannard, iEEGLAB, 2026

g = struct('fieldtrip', [], 'crp', [], 'erdetect', [], 'ica', [], 'carla', []);
if mod(numel(varargin), 2), error('run_validations:args', 'Options come in name-value pairs.'); end
for k = 1:2:numel(varargin)
    f = lower(char(varargin{k}));
    if ~isfield(g, f), error('run_validations:option', 'Unknown validation ''%s''.', f); end
    g.(f) = varargin{k + 1};
end
if ~isfolder(outDir), mkdir(outDir); end
addpath(fileparts(mfilename('fullpath')));
vis = get(0, 'DefaultFigureVisible'); set(0, 'DefaultFigureVisible', 'off');
cleanup = onCleanup(@() set(0, 'DefaultFigureVisible', vis));
tsv = {'FileType', 'text', 'Delimiter', '\t'};

if ~isempty(g.fieldtrip)
    out = validate_fieldtrip_uci29(g.fieldtrip);
    T = out.stats; out.stats = []; out = rmfield(out, 'stats_info');
    save(fullfile(outDir, 'fieldtrip_uci29_figure.mat'), '-struct', 'out', '-v7');
    writetable(T, fullfile(outDir, 'fieldtrip_uci29_band_stats.tsv'), tsv{:});
end
if ~isempty(g.crp)
    T = validate_crp_scripts(g.crp);
    writetable(T, fullfile(outDir, 'crp_vs_crp_scripts.tsv'), tsv{:});
end
if ~isempty(g.erdetect)
    out = validate_erdetect_ds004080(g.erdetect{:}); %#ok<NASGU> saved by name
    save(fullfile(outDir, 'validate_erdetect_ds004080_out.mat'), 'out');
end
if ~isempty(g.ica)
    out = validate_icaref_michelmann(g.ica); %#ok<NASGU> saved by name
    save(fullfile(outDir, 'validate_icaref_michelmann_out.mat'), 'out');
end
if ~isempty(g.carla)
    out = validate_carla_huang_ds004977(g.carla{:});
    fid = fopen(fullfile(outDir, 'carla_fig5.json'), 'w');
    fwrite(fid, jsonencode(out.figure)); fclose(fid);
end
fprintf('Validation outputs in %s\n', outDir);
end
