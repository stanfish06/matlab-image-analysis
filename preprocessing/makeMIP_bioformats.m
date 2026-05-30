function MIP = makeMIP_bioformats(...
    filename, barefname, position, channel, outputdir, saveidx, tmax, zrange)
%MAKEMIP_BIOFORMATS Maximum-intensity projection of one channel via Bio-Formats.
%   MIP = MAKEMIP_BIOFORMATS(FILENAME, BAREFNAME, POSITION, CHANNEL, OUTPUTDIR,
%   SAVEIDX, TMAX, ZRANGE) reads FILENAME with Bio-Formats and computes the
%   z-MIP of channel CHANNEL (0-based) for each time point up to TMAX over
%   z-slices ZRANGE. Output is written to OUTPUTDIR (default <dataDir>/MIP);
%   if SAVEIDX is true the per-pixel argmax-z index image is also saved.
%   Loops over series for .lif/.nd2 files.

    r = bfGetReader(filename);
    dataDir = fileparts(filename);

    % some input checking
    if ~exist('saveidx','var')
        saveidx = false;
    end
    if ~exist('outputdir','var') || isempty(outputdir)
        outputdir = fullfile(dataDir,'MIP');
    end
    if ~exist(outputdir,'dir')
        mkdir(outputdir);
    end
    
    % we only want to loop over series for lif files
    % because vsi files have the same picture at different resolution as
    % the series
    if strcmp(filename(end-3:end),'.lif') || strcmp(filename(end-3:end),'.nd2')
        N = r.getSeriesCount;
    else
        N = 1;
    end
    
    for si = 1:N
        
        disp(['reading series ' num2str(si)]);
        r.setSeries(si-1);
        
        ci = channel;
        if ci+1 > r.getSizeC()
            error('channel exceeds number of channels');
        end
        if ~exist('tmax','var') || isempty(tmax)
            tmax = r.getSizeT();
        else
            tmax = min(r.getSizeT(), tmax);
        end
        if ~exist('zrange','var') || isempty(zrange)
            zrange = 1:r.getSizeZ();
        end

        MIP = zeros([r.getSizeY() r.getSizeX() tmax], 'uint16');
        MIPidx = zeros([r.getSizeY() r.getSizeX() tmax], 'uint8');
        maxX = 4000;
        maxY = 4000;
        nXblocks = ceil(r.getSizeX() / maxX);
        nYblocks = ceil(r.getSizeY() / maxY);
        readSizeX = ceil(r.getSizeX() / nXblocks);
        readSizeY = ceil(r.getSizeY() / nYblocks);       
        fprintf([sprintf('%.4d planes to process', r.getSizeZ()), newline]);
        for ti = 1:tmax
            if numel(zrange) > 1
                im_mip = zeros([r.getSizeY() r.getSizeX()]);
                im_mip_idx = ones([r.getSizeY() r.getSizeX()]);
                for zi = 1:numel(zrange)
                    fprintf([sprintf('process plane %.4d', zi), newline]);
                    im_zi = zeros([r.getSizeY() r.getSizeX()]);
                    xo = 1;
                    yo = 1;
                    for xi = 1:nXblocks
                        readSizeX_adj = readSizeX - max(((xo + readSizeX) - r.getSizeX()), 0);
                        for yi = 1:nYblocks
                            fprintf('.');
                            readSizeY_adj = readSizeY - max(((yo + readSizeY) - r.getSizeY()), 0);
                            im_zi_sub = bfGetPlane(r, r.getIndex(zrange(zi)-1,ci,ti-1)+1, xo, yo, readSizeX_adj, readSizeY_adj);
                            im_zi((yo):(yo+readSizeY_adj-1),(xo):(xo+readSizeX_adj-1)) = im_zi_sub;
                            yo = yo + readSizeY_adj;
                        end
                        yo = 1;
                        xo = xo + readSizeX_adj;
                        fprintf(newline);
                    end
                    if zi == 1
                        im_mip = im_zi;
                    else
                        [im_mip, im_mip_idx_lin] = max(cat(3, im_zi, im_mip),[],3,"linear");
                        im_mip_idx_tmp = cat(3, ones([r.getSizeY() r.getSizeX()]) * zi, im_mip_idx);
                        img_mip_idx = im_mip_idx_tmp(im_mip_idx_lin);
                    end
                end
                MIP(:,:,ti) = im_mip;
                MIPidx(:,:,ti) = img_mip_idx;
            else
                MIP(:,:,ti) = bfGetPlane(r, r.getIndex(0,ci,ti-1)+1);
            end
            fprintf('.');
            if mod(ti,60)==0
                fprintf('\n');
            end
        end
        fprintf('\n');

        % save result
        %-------------

        % MIP
        
        if r.getSeriesCount > 1 && (strcmp(filename(end-3:end),'.lif') || strcmp(filename(end-3:end),'.nd2'))
            pi = si-1;
            warning('multiple series lif, assuming all positions are in single file');
        else
            pi = position;
        end

        if ~isempty(pi)
            fname = fullfile(outputdir, sprintf([barefname '_MIP_p%.4d_w%.4d.tif'],pi,ci));
            idxfname = fullfile(outputdir, sprintf([barefname '_MIPidx_p%.4d_w%.4d.tif'],pi,ci));
        else
            fname = fullfile(outputdir, sprintf([barefname '_MIP_w%.4d.tif'],ci));
            idxfname = fullfile(outputdir, sprintf([barefname '_MIPidx_w%.4d.tif'],ci));
        end
        if exist(fname,'file')
            delete(fname);
        end

        MIP = squeeze(MIP);
        disp('saving MIP');
        imwrite(MIP(:,:,1), fname);
        fprintf('.');
        for i = 2:size(MIP,3)
            fprintf('.');
            if mod(i,60)==0
                fprintf('\n');
            end
            imwrite(MIP(:,:,i), fname,'WriteMode','Append');
        end
        fprintf('\n');

        if saveidx
            disp('saving MIPidx');
            MIPidx = squeeze(MIPidx);
            fprintf('.');
            imwrite(MIPidx(:,:,1), idxfname);
            for i = 2:size(MIPidx,3)
                fprintf('.');
                if mod(i,60)==0
                    fprintf('\n');
                end
                imwrite(MIPidx(:,:,i), idxfname,'WriteMode','Append');
            end
            fprintf('\n');
        end
    end
	r.close();
end
