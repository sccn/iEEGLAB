function g = ieeglab_read_gifti(f)
% ieeglab_read_gifti() - Read a GIfTI surface (.surf.gii) without vistasoft or SPM.
%
% Usage:
%   g = ieeglab_read_gifti('pial.L.surf.gii');
%   % g.vertices [N x 3] double, g.faces [M x 3] double (1-based)
%
% Handles the three GIfTI encodings (ASCII, Base64Binary, GZipBase64Binary),
% both array orders and both endiannesses. The vertices are the float array
% with intent POINTSET (or, failing that, the first float N x 3 array) and the
% faces the integer array with intent TRIANGLE: files do not always store the
% vertices first.
%
% Replaces the @gifti class that iEEGLAB used to borrow from vistasoft
% (issue #8), so the 3D plots need no other toolbox.
%
% Cedric Cannard, iEEGLAB, 2026

txt = fileread(f);
blocks = regexp(txt, '<DataArray([^>]*)>(.*?)</DataArray>', 'tokens');
if isempty(blocks)
    error('ieeglab_read_gifti:noData', 'No DataArray in %s: not a GIfTI file?', f);
end
g = struct('vertices', [], 'faces', []);
for b = 1:numel(blocks)
    a = local_attributes(blocks{b}{1});
    d = regexp(blocks{b}{2}, '<Data>(.*?)</Data>', 'tokens', 'once');
    if isempty(d), continue; end
    x = local_decode(strtrim(d{1}), a, f);
    intent = upper(local_get(a, 'Intent', ''));
    if contains(intent, 'TRIANGLE')
        g.faces = double(x) + 1;
    elseif contains(intent, 'POINTSET') || (isempty(g.vertices) && isfloat(x) && size(x, 2) == 3)
        g.vertices = double(x);
    end
end
if isempty(g.vertices) || isempty(g.faces)
    error('ieeglab_read_gifti:incomplete', 'Found no vertex or no triangle array in %s.', f);
end
if max(g.faces(:)) > size(g.vertices, 1) || min(g.faces(:)) < 1
    error('ieeglab_read_gifti:badFaces', 'Triangle indices in %s do not match its %d vertices.', f, size(g.vertices, 1));
end
end

function a = local_attributes(s)
kv = regexp(s, '(\w+)\s*=\s*"([^"]*)"', 'tokens');
a = struct();
for k = 1:numel(kv), a.(kv{k}{1}) = kv{k}{2}; end
end

function v = local_get(a, name, default)
if isfield(a, name), v = a.(name); else, v = default; end
end

function x = local_decode(s, a, f)
types = struct('NIFTI_TYPE_UINT8', 'uint8', 'NIFTI_TYPE_INT8', 'int8', 'NIFTI_TYPE_INT16', 'int16', ...
    'NIFTI_TYPE_UINT16', 'uint16', 'NIFTI_TYPE_INT32', 'int32', 'NIFTI_TYPE_UINT32', 'uint32', ...
    'NIFTI_TYPE_FLOAT32', 'single', 'NIFTI_TYPE_FLOAT64', 'double');
dt = local_get(a, 'DataType', '');
if ~isfield(types, dt), error('ieeglab_read_gifti:type', 'Unsupported GIfTI DataType %s in %s.', dt, f); end
cls = types.(dt);
nd = str2double(local_get(a, 'Dimensionality', '1'));
dims = arrayfun(@(k) str2double(local_get(a, sprintf('Dim%d', k), '1')), 0:nd-1);
if isscalar(dims), dims(2) = 1; end
switch local_get(a, 'Encoding', 'ASCII')
    case 'ASCII'
        x = cast(sscanf(s, '%f'), cls);
    case {'Base64Binary', 'GZipBase64Binary'}
        bytes = matlab.net.base64decode(s);
        if strcmp(a.Encoding, 'GZipBase64Binary'), bytes = local_inflate(bytes); end
        x = typecast(uint8(bytes(:))', cls);
        [~, ~, native] = computer;
        if strcmpi(local_get(a, 'Endian', 'LittleEndian'), 'BigEndian') ~= strcmp(native, 'B')
            x = swapbytes(x);
        end
    otherwise
        error('ieeglab_read_gifti:encoding', 'Unsupported GIfTI encoding %s in %s (external files are not read).', a.Encoding, f);
end
if numel(x) ~= prod(dims)
    error('ieeglab_read_gifti:size', 'A data array of %s holds %d values, expected %d.', f, numel(x), prod(dims));
end
if strcmp(local_get(a, 'ArrayIndexingOrder', 'RowMajorOrder'), 'RowMajorOrder')
    x = permute(reshape(x, fliplr(dims)), numel(dims):-1:1);
else
    x = reshape(x, dims);
end
end

function out = local_inflate(bytes)
% zlib stream (GIfTI's GZipBase64Binary), through the JVM that ships with MATLAB
% (Java arrays reach MATLAB as copies, so the stream is copied on the Java side.)
in  = java.util.zip.InflaterInputStream(java.io.ByteArrayInputStream(typecast(uint8(bytes(:)), 'int8')));
buf = java.io.ByteArrayOutputStream();
com.mathworks.mlwidgets.io.InterruptibleStreamCopier.getInterruptibleStreamCopier().copyStream(in, buf);
in.close();
out = typecast(int8(buf.toByteArray()), 'uint8');
end
