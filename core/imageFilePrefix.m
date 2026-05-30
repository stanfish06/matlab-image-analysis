function prefix = imageFilePrefix(meta, posIndex, channel, time, dataDir)
%IMAGEFILEPREFIX Filename prefix (no extension) for one position's processed image.
%   PREFIX = IMAGEFILEPREFIX(META, POSINDEX, CHANNEL, TIME) returns the base
%   filename (without extension) of the preprocessed image for position
%   POSINDEX (1-based), channel CHANNEL (0-based) and time point TIME (0-based).
%   Append '.tif', '_cp_masks.tif', '_FinalSegmentation.tif' or '_masks.mat' as
%   needed.
%
%   It transparently handles both preprocessing outputs, so the my-pipeline
%   steps run on stitched OR non-stitched data without editing hard-coded
%   'stitched_' filenames:
%     - stitched data (stitchImageMontages): META.fileNames is empty and
%       META.filenameFormat is 'stitched_p%.4d_w%.4d_t%.4d.tif', giving
%       'stitched_p####_w####_t####'.
%     - non-stitched z-stacks (writeZstacks): META.fileNames holds the raw file
%       names; the prefix keeps each raw file's base name, e.g.
%       'myexp_p####_w####_t####' (or 'myexp_w####_t####' if there is no
%       per-file position index).
%
%   POSINDEX maps to the on-disk index the same way QUANTIFY does (stitched:
%   POSINDEX-1; non-stitched: the index parsed from the raw file name), so this
%   matches the mask/positions files QUANTIFY reads and writes.
%
%   PREFIX = IMAGEFILEPREFIX(META, POSINDEX, CHANNEL, TIME, DATADIR) uses DATADIR
%   to resolve raw file names (default: current directory).
%
%   See also QUANTIFY, WRITEZSTACKS, STITCHIMAGEMONTAGES, PARSEFILENAME.

    if nargin < 5 || isempty(dataDir)
        dataDir = pwd;
    end

    if isempty(meta.fileNames)
        % stitched (or fixed-format) data: build from the filename format
        [~, bare] = fileparts(meta.filenameFormat);
        id = posIndex - 1;                       % on-disk index is 0-based
        prefix = sprintf(bare, id, channel, time);
    else
        % non-stitched z-stacks: keep the raw file's base name
        [barefname, id] = parseFilename(meta.fileNames{posIndex}, dataDir);
        if isempty(id)
            prefix = sprintf([barefname '_w%.4d_t%.4d'], channel, time);
        else
            prefix = sprintf([barefname '_p%.4d_w%.4d_t%.4d'], id, channel, time);
        end
    end
end
