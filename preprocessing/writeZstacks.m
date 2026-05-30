function writeZstacks(dataDir,opts)
%function to write tiff stacks of z slices for segmentation
%%% TODO - handle all positions being saved in a single multi-series nd2
%%% file that is not used for stitching %%%

if isfield(opts,'writeDir')
    %can specify a separate folder to which to write the z slices if
    %desired
    writeDir = opts.writeDir;
    if ~exist(writeDir,'dir'), mkdir(writeDir); end
else
    %by default write the z stacks in the data directory
    writeDir = dataDir;
end

if isfield(opts,'exportZarr')
    exportZarr = opts.exportZarr;
else
    exportZarr = false;
end
if exportZarr
    if isfield(opts,'zarrDir') && ~isempty(opts.zarrDir)
        zarrDir = opts.zarrDir;
    else
        zarrDir = fullfile(writeDir, 'Zarr');
    end
    if ~exist(zarrDir,'dir'), mkdir(zarrDir); end
    if isfield(opts,'zarrOptions')
        zarrOptions = opts.zarrOptions;
    else
        zarrOptions = struct();
    end
else
    zarrOptions = struct();
end

if isfield(opts, 'writeZsubStacks')
    writeZsubStacks = opts.writeZsubStacks;
    maxZslices = opts.maxZslices;
else
    writeZsubStacks = false;
    maxZslices = 0;
end

if isfield(opts,'filepattern') && ~isempty(opts.filepattern)
    %manually specify the image file extension
    filepattern = opts.filepattern;
else
    %look for files from a list of microscope image file formats
    filepatterns = {'*FusionStitcher*.ims','*.ims','*.nd2','*.lif','*.vsi'};
    idx = 1;
    criterion = true;
    while criterion && idx <= length(filepatterns)
        filepattern = filepatterns{idx};
        listing = dir(fullfile(dataDir,filepattern));
        if ~isempty(listing)
            criterion = false;
        end
        idx = idx + 1;
    end
    disp(strcat("writing stacks for files matching pattern ",filepattern))
end

if isfield(opts,'nucChannel')
    nucChannel = opts.nucChannel;
else
    nucChannel = 0;
end

if ~isfield(opts,'tmax') || isempty(opts.tmax)
    tmax = Inf;
else
    tmax = opts.tmax;
end

%input additional channels for which to write z slices, indexed from 0
if ~isfield(opts,'writeChannels')
    writeChannels = [];
else
    writeChannels = opts.writeChannels;
end

%always write z slices for the nuclear channel
writeChannels = union(nucChannel,writeChannels);
nChannels = numel(writeChannels);

filelist = dir(fullfile(dataDir,filepattern));
%do we also want to find image files in subfolders like for Jo's paper?
%that would be implemented as follows:
%filelist = dir(fullfile(dataDir,['**/',filepattern]));
nss = ones(1,numel(filelist));
for fi = 1:numel(filelist)

    %name of the z slice file should follow
    fname = filelist(fi).name;
    [~,~,ext] = fileparts(fname);
    [barefname, pidx] = parseFilename(fname,dataDir);

    disp(['writing z slices for ' fname]);
    tic
    r = bfGetReader(fullfile(dataDir, fname));
    r.setSeries(0);

    nZslices = r.getSizeZ();
    nt = r.getSizeT();

    if strcmp(ext,'.nd2') || strcmp(ext,'.lif')
        %nd2 and lif files may have multiple fields of view saved as
        %different series of one file
        ns = r.getSeriesCount;
    else
        %ims files may have multiple series for the same field of view that
        %just consist of downsampled versions of the same data
        ns = 1;
    end
    nss(fi) = ns;

    %account for nd2 and lif files with multiple series here
    for si = 1:ns
        if isempty(pidx) && ns > 1
            id = si;
        else
            id = pidx;
        end

        r.setSeries(si-1);

        if isempty(id)
            outfnamePattern = [barefname, '_w%.4d_t%.4d.tif'];
        else
            outfnamePattern = [barefname, sprintf('_p%.4d',id), '_w%.4d_t%.4d.tif'];
        end

        if writeZsubStacks
            if maxZslices == 1
                nSubStacks = nZslices;
                outfnamePattern = [barefname, '_w%.4d_t%.4d_z%.4d.tif'];
            else
                nSubStacks = ceil(nZslices / maxZslices);
                outfnamePattern = [barefname, '_w%.4d_t%.4d_z%.4d-%.4d.tif'];
            end
        else
            nSubStacks = 1;
        end
        maxX = 4000;
        maxY = 4000;
        nXblocks = ceil(r.getSizeX() / maxX);
        nYblocks = ceil(r.getSizeY() / maxY);
        readSizeX = ceil(r.getSizeX() / nXblocks);
        readSizeY = ceil(r.getSizeY() / nYblocks);
        for ti = 1:min(tmax,nt)
            for cii = 1:nChannels
                ci = writeChannels(cii);
                for bi = 1:nSubStacks
                    if writeZsubStacks
                        zl = (bi - 1) * maxZslices;
                        zr = min((bi * maxZslices) - 1, nZslices - 1);
                    else
                        zl = 0;
                        zr = nZslices - 1;
                    end
                    if exportZarr
                        stackForZarr = zeros(r.getSizeY(), r.getSizeX(), zr-zl+1, 'uint16');
                    end
                    for zi = zl:zr
                        fprintf([sprintf('process plane %.4d', zi), newline]);
                        if zi == zl
                            mode = 'overwrite';
                        else
                            mode = 'append';
                        end
                        
                        xo = 1;
                        yo = 1;
                        im_zi = zeros(r.getSizeY(), r.getSizeX(), 'uint16');
                        for xi = 1:nXblocks
                            readSizeX_adj = readSizeX - max(((xo + readSizeX) - r.getSizeX()), 0);
                            for yi = 1:nYblocks
                                fprintf('.');
                                readSizeY_adj = readSizeY - max(((yo + readSizeY) - r.getSizeY()), 0);
                                im_zi_sub = bfGetPlane(r, r.getIndex(zi,ci,ti-1)+1, xo, yo, readSizeX_adj, readSizeY_adj);
                                im_zi((yo):(yo+readSizeY_adj-1),(xo):(xo+readSizeX_adj-1)) = im_zi_sub;
                                yo = yo + readSizeY_adj;
                            end
                            yo = 1;
                            xo = xo + readSizeX_adj;
                            fprintf(newline);
                        end
                        %should we check for blank/mostly blank z slices and skip?
                        if writeZsubStacks
                            if maxZslices == 1
                                writename = fullfile(writeDir,sprintf(outfnamePattern,ci,ti-1,zl));
                            else
                                writename = fullfile(writeDir,sprintf(outfnamePattern,ci,ti-1,zl,zr));
                            end
                        else
                            writename = fullfile(writeDir,sprintf(outfnamePattern,ci,ti-1));
                        end
                        imwrite(uint16(im_zi),writename,'WriteMode',mode)
                        if exportZarr
                            stackForZarr(:,:,zi-zl+1) = uint16(im_zi);
                        end
                    end
                    if exportZarr
                        [~, zarrBareName] = fileparts(writename);
                        zarrname = fullfile(zarrDir, [zarrBareName, '.zarr']);
                        writeZarrStack(stackForZarr, zarrname, addImageScale(zarrOptions, dataDir));
                    end
                end
            end
        end
    end
    toc
end

load(fullfile(dataDir,'meta.mat'),'meta')
meta.fileNames = {filelist.name};
%sum of the number of series in each file should give total positions
meta.nPositions = sum(nss);
save(fullfile(dataDir,'meta.mat'),'meta')


end

function zarrOptions = addImageScale(zarrOptions, dataDir)
if ~isfield(zarrOptions, 'scale') || isempty(zarrOptions.scale)
    metaFile = fullfile(dataDir, 'meta.mat');
    if exist(metaFile, 'file')
        S = load(metaFile, 'meta');
        zres = 1;
        hasZres = (isobject(S.meta) && isprop(S.meta, 'zres')) || (isstruct(S.meta) && isfield(S.meta, 'zres'));
        if hasZres && ~isempty(S.meta.zres) && ~isnan(S.meta.zres)
            zres = S.meta.zres;
        end
        zarrOptions.scale = [1 1 zres S.meta.yres S.meta.xres];
    end
end
end
