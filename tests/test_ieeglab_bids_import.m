function tests = test_ieeglab_bids_import
% Loading after EEGLAB's BIDS import (EEG-BIDS), and the built-in GIfTI reader.
% Headless. Uses the native-rate tutorial datasets in tutorial/ieeglab_tutorial_*.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(mfilename('fullpath')));
tc.assumeNotEmpty(which('eeglab'), 'EEGLAB is not on the MATLAB path.');
if isempty(which('pop_loadset')), evalc('eeglab nogui'); end
addpath(root); addpath(fullfile(root, 'functions'));
tc.TestData.root = root;
end

function test_gifti_reader_matches_mesh(tc)
f = fullfile(tc.TestData.root, 'tutorial', 'dataset_seeg', 'pial.L.surf.gii');
g = ieeglab_read_gifti(f);
tc.verifySize(g.vertices, [155488 3]);
tc.verifyEqual(size(g.faces, 2), 3);
tc.verifyEqual(min(g.faces(:)), 1);
tc.verifyEqual(max(g.faces(:)), size(g.vertices, 1));
tc.verifyTrue(all(abs(g.vertices(:)) < 200), 'Vertices must be in millimetres.');
end

function test_gifti_reader_gzip_ascii_any_order(tc)
% Triangles listed before vertices, ASCII faces, zlib-compressed float vertices
V = single([0 0 0; 1 0 0; 0 1 0; 0 0 1]); F = int32([0 1 2; 0 1 3; 0 2 3; 1 2 3]);
bs = java.io.ByteArrayOutputStream(); z = java.util.zip.DeflaterOutputStream(bs);
z.write(typecast(reshape(V', 1, []), 'int8')); z.close();
v64 = matlab.net.base64encode(typecast(int8(bs.toByteArray()), 'uint8'));
xml = sprintf(['<?xml version="1.0"?><GIFTI NumberOfDataArrays="2">' ...
  '<DataArray Intent="NIFTI_INTENT_TRIANGLE" DataType="NIFTI_TYPE_INT32" ArrayIndexingOrder="RowMajorOrder" ' ...
  'Dimensionality="2" Dim0="4" Dim1="3" Encoding="ASCII" Endian="LittleEndian"><Data>%s</Data></DataArray>' ...
  '<DataArray Intent="NIFTI_INTENT_POINTSET" DataType="NIFTI_TYPE_FLOAT32" ArrayIndexingOrder="RowMajorOrder" ' ...
  'Dimensionality="2" Dim0="4" Dim1="3" Encoding="GZipBase64Binary" Endian="LittleEndian"><Data>%s</Data></DataArray>' ...
  '</GIFTI>'], sprintf('%d ', reshape(F', 1, [])), v64);
f = [tempname '.surf.gii'];
fid = fopen(f, 'w'); fwrite(fid, xml); fclose(fid);
c = onCleanup(@() delete(f));
g = ieeglab_read_gifti(f);
tc.verifyEqual(g.vertices, double(V));
tc.verifyEqual(g.faces, double(F) + 1);
end

function test_sidecars_found_from_derivatives_folder(tc)
% pop_importbids saves to <root>/derivatives/eeglab/sub-.../ieeg by default;
% the sidecars stay in <root>/sub-.../ieeg.
tut = fullfile(tc.TestData.root, 'tutorial', 'ieeglab_tutorial_seeg');
tc.assumeTrue(isfolder(tut), 'Native-rate tutorial dataset not present.');
src = tempname; mkdir(fullfile(src, 'sub-02', 'ses-ieeg01'));
c = onCleanup(@() rmdir(src, 's'));
copyfile(fullfile(tut, 'sub-02', 'ses-ieeg01', 'ieeg', '*.tsv'), fullfile(src, 'sub-02', 'ses-ieeg01', 'ieeg'));
EEG = struct('filepath', fullfile(src, 'derivatives', 'eeglab', 'sub-02', 'ses-ieeg01', 'ieeg'), ...
             'filename', 'sub-02_ses-ieeg01_task-ccep_run-01_ieeg.set');
mkdir(EEG.filepath);
p = ieeglab_bids_sibling(EEG, 'channels');
tc.verifyEqual(p, fullfile(src, 'sub-02', 'ses-ieeg01', 'ieeg', 'sub-02_ses-ieeg01_task-ccep_run-01_channels.tsv'));
% and from a source path recorded by the importer, wherever the output went
EEG = struct('filepath', tempdir, 'filename', 'x.set', 'BIDS', struct('sourcefile', ...
    fullfile(src, 'sub-02', 'ses-ieeg01', 'ieeg', 'sub-02_ses-ieeg01_task-ccep_run-01_ieeg.set')));
tc.verifyEqual(ieeglab_bids_sibling(EEG, 'events'), ...
    fullfile(src, 'sub-02', 'ses-ieeg01', 'ieeg', 'sub-02_ses-ieeg01_task-ccep_run-01_events.tsv'));
end

function test_load_after_eeglab_bids_import(tc)
% EEG-BIDS moves the site column into EEG.event.type and (up to at least
% 2026-05) assigns electrodes.tsv rows to channels by position, which on this
% subject gives most contacts the wrong label.
tc.assumeNotEmpty(which('pop_importbids'), 'EEG-BIDS plugin not installed.');
src = fullfile(tc.TestData.root, 'tutorial', 'ieeglab_tutorial_seeg');
tc.assumeTrue(isfolder(src), 'Native-rate tutorial dataset not present.');
work = tempname; copyfile(src, work); c = onCleanup(@() rmdir(work, 's'));
[~, ~, ALLEEG] = evalc(['pop_importbids(work, ''bidsevent'', ''on'', ''bidschanloc'', ''on'', ' ...
    '''eventtype'', ''electrical_stimulation_site'')']);
d = dir(fullfile(src, 'sub-*', 'ses-*', 'ieeg', '*_ieeg.set'));
O = pop_loadset('filename', d.name, 'filepath', d.folder, 'verbose', 'off');
[~, E] = evalc('ieeglab_load(ALLEEG(1), struct())');
[tf, loc] = ismember({E.chanlocs.labels}, {O.chanlocs.labels});
tc.verifyTrue(all(tf));
same = arrayfun(@(k) isequal(E.data(k, 1:2000), O.data(loc(k), 1:2000)), 1:E.nbchan);
tc.verifyTrue(all(same), 'Every channel label must name its own data after loading.');
bad = {E.chanlocs(cellfun(@(x) strcmpi(char(x), 'bad'), {E.chanlocs.status})).labels};
tc.verifyEqual(sort(bad), {'RA15','RB14','RB15'});
ty = unique(cellfun(@(x) char(string(x)), {E.event.type}, 'UniformOutput', false));
tc.verifyEqual(ty, {'RA2-RA3','RA3-RA4','RA4-RA5'});
end
