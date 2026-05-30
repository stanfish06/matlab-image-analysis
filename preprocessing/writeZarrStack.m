function writeZarrStack(stack, zarrPath, opts)
%WRITEZARRSTACK Write a Y-X-Z image stack as uncompressed OME-Zarr v2.
%   WRITEZARRSTACK(STACK, ZARRPATH) writes STACK to ZARRPATH. STACK is
%   interpreted as Y-by-X or Y-by-X-by-Z and stored as a single-scale OME-Zarr
%   array with axes t-c-z-y-x and singleton t/c dimensions.
%
%   WRITEZARRSTACK(..., OPTS) accepts:
%       chunks      1x5 chunk vector in t-c-z-y-x order.
%       axes        currently 'tczyx' (default).
%       scale       physical scale in t-c-z-y-x order (default ones).
%       translation physical translation in t-c-z-y-x order (default zeros).

if nargin < 3 || isempty(opts)
    opts = struct();
end
if ~isfield(opts, 'axes') || isempty(opts.axes)
    opts.axes = 'tczyx';
end
if ~strcmp(opts.axes, 'tczyx')
    error('writeZarrStack:UnsupportedAxes', 'Only t-c-z-y-x OME-Zarr export is currently supported.');
end

if ndims(stack) > 3
    error('writeZarrStack:UnsupportedStackShape', 'Expected a Y-X or Y-X-Z stack.');
end
if ~isa(stack, 'uint8') && ~isa(stack, 'uint16') && ~isa(stack, 'single') && ~isa(stack, 'double')
    error('writeZarrStack:UnsupportedDType', 'Supported stack classes are uint8, uint16, single, and double.');
end

stackSize = size(stack);
if numel(stackSize) < 3
    stackSize(3) = 1;
end
ySize = stackSize(1);
xSize = stackSize(2);
zSize = stackSize(3);
shape = [1 1 zSize ySize xSize];

if ~isfield(opts, 'chunks') || isempty(opts.chunks)
    opts.chunks = [1 1 min(zSize, 16) min(ySize, 512) min(xSize, 512)];
end
chunks = opts.chunks;
if numel(chunks) ~= 5
    error('writeZarrStack:InvalidChunks', 'opts.chunks must be a 1x5 vector in t-c-z-y-x order.');
end
chunks = min(chunks(:)', shape);
chunks = max(chunks, 1);

if ~isfield(opts, 'scale') || isempty(opts.scale)
    opts.scale = ones(1,5);
end
if ~isfield(opts, 'translation') || isempty(opts.translation)
    opts.translation = zeros(1,5);
end

if exist(zarrPath, 'dir')
    rmdir(zarrPath, 's');
end
mkdir(zarrPath);
arrayPath = fullfile(zarrPath, '0');
mkdir(arrayPath);

writeJson(fullfile(zarrPath, '.zgroup'), struct('zarr_format', 2));
writeOmeZarrAttrs(fullfile(zarrPath, '.zattrs'), opts);
writeZarray(fullfile(arrayPath, '.zarray'), shape, chunks, zarrDType(stack), cast(0, class(stack)));
writeJson(fullfile(arrayPath, '.zattrs'), struct());

zChunkStarts = 1:chunks(3):zSize;
yChunkStarts = 1:chunks(4):ySize;
xChunkStarts = 1:chunks(5):xSize;

for zci = 1:numel(zChunkStarts)
    zRange = zChunkStarts(zci):min(zChunkStarts(zci) + chunks(3) - 1, zSize);
    for yci = 1:numel(yChunkStarts)
        yRange = yChunkStarts(yci):min(yChunkStarts(yci) + chunks(4) - 1, ySize);
        for xci = 1:numel(xChunkStarts)
            xRange = xChunkStarts(xci):min(xChunkStarts(xci) + chunks(5) - 1, xSize);
            chunk = zeros(chunks(4), chunks(5), chunks(3), class(stack));
            chunk(1:numel(yRange), 1:numel(xRange), 1:numel(zRange)) = stack(yRange, xRange, zRange);
            chunk = permute(chunk, [2 1 3]); % C-order bytes for z-y-x
            chunkName = sprintf('0.0.%d.%d.%d', zci-1, yci-1, xci-1);
            writeChunk(fullfile(arrayPath, chunkName), chunk);
        end
    end
end

end

function writeOmeZarrAttrs(path, opts)
fid = fopen(path, 'w');
if fid == -1
    error('writeZarrStack:FileOpenFailed', 'Could not write %s.', path);
end
cleanupObj = onCleanup(@() fclose(fid));

fprintf(fid, '{\n');
fprintf(fid, '  "multiscales": [\n');
fprintf(fid, '    {\n');
fprintf(fid, '      "version": "0.4",\n');
fprintf(fid, '      "name": "image",\n');
fprintf(fid, '      "axes": [\n');
fprintf(fid, '        {"name": "t", "type": "time"},\n');
fprintf(fid, '        {"name": "c", "type": "channel"},\n');
fprintf(fid, '        {"name": "z", "type": "space"},\n');
fprintf(fid, '        {"name": "y", "type": "space"},\n');
fprintf(fid, '        {"name": "x", "type": "space"}\n');
fprintf(fid, '      ],\n');
fprintf(fid, '      "datasets": [\n');
fprintf(fid, '        {\n');
fprintf(fid, '          "path": "0",\n');
fprintf(fid, '          "coordinateTransformations": [\n');
fprintf(fid, '            {"type": "scale", "scale": %s}', jsonencode(opts.scale));
if any(opts.translation ~= 0)
    fprintf(fid, ',\n');
    fprintf(fid, '            {"type": "translation", "translation": %s}', jsonencode(opts.translation));
end
fprintf(fid, '\n');
fprintf(fid, '          ]\n');
fprintf(fid, '        }\n');
fprintf(fid, '      ]\n');
fprintf(fid, '    }\n');
fprintf(fid, '  ]\n');
fprintf(fid, '}\n');
end

function dtype = zarrDType(stack)
switch class(stack)
    case 'uint8'
        dtype = '|u1';
    case 'uint16'
        dtype = '<u2';
    case 'single'
        dtype = '<f4';
    case 'double'
        dtype = '<f8';
end
end

function writeZarray(path, shape, chunks, dtype, fillValue)
fid = fopen(path, 'w');
if fid == -1
    error('writeZarrStack:FileOpenFailed', 'Could not write %s.', path);
end
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, '{\n');
fprintf(fid, '  "zarr_format": 2,\n');
fprintf(fid, '  "shape": %s,\n', jsonencode(shape));
fprintf(fid, '  "chunks": %s,\n', jsonencode(chunks));
fprintf(fid, '  "dtype": "%s",\n', dtype);
fprintf(fid, '  "compressor": null,\n');
fprintf(fid, '  "fill_value": %s,\n', jsonencode(fillValue));
fprintf(fid, '  "order": "C",\n');
fprintf(fid, '  "filters": null\n');
fprintf(fid, '}\n');
end

function writeJson(path, value)
fid = fopen(path, 'w');
if fid == -1
    error('writeZarrStack:FileOpenFailed', 'Could not write %s.', path);
end
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', jsonencode(value, 'PrettyPrint', true));
end

function writeChunk(path, chunk)
fid = fopen(path, 'w', 'ieee-le');
if fid == -1
    error('writeZarrStack:FileOpenFailed', 'Could not write %s.', path);
end
cleanupObj = onCleanup(@() fclose(fid));
fwrite(fid, chunk(:), class(chunk));
end
