% ANALYZE_LIVE_STAN  Live time-lapse quantification/analysis driver (Stan's) [template].
%   Loads meta.mat and reads out per-cell intensities over positions, time and
%   channels (optionally on MIPs) into a time-indexed cellData, for downstream
%   signaling plots (see plot_marker_dynamics). Edit bare/imageType/useMIP near
%   the top; run from the data directory.
clear; close all; clc

scriptPath = pwd;
dataDir = scriptPath;

%% setup
load(fullfile(dataDir,'meta.mat'),'meta')
channelLabel = meta.channelLabel;
nucChannel = meta.nucChannel;
npos = meta.nPositions;
ntime = meta.nTime;
nchannels = meta.nChannels;
bare = '20260204_BMP_actd_chx_sox2-h2b_p%.4d_w0001_t%.4d'

% define imageType as micropatterned 'MP' or 'disordered'
imageType = 'disordered';
% radii of micropatterns in micron for each condition (ignored if not 'MP')
radii = 700*ones(1,meta.nWells)/2;

% use maximum intensity projection (MIP) across z before computing intensity values
useMIP = true; 

segDir = dataDir;

%% iterate over positions, time points, and channels, and read out intensities

if strcmp(imageType,'MP')
    positions(meta.nPositions) = Colony();
elseif strcmp(imageType,'disordered')
    positions(meta.nPositions) = Position();
end

tic
for pi = 1:meta.nPositions
    condi = find(pi >= meta.conditionStartPos,1,'last'); % condition index
          
    fprintf('position %d of %d\n',pi,npos)

    if strcmp(imageType,'MP')
        positions(pi) = Colony(meta, pi);
        positions(pi).setRadius(radii(condi), meta.xres);
        positions(pi).well = condi;
    elseif strcmp(imageType,'disordered')
        positions(pi) = Position(meta, pi);
    end
    
    % shared properties between micropattern and disordered
    positions(pi).dataChannels = 1:positions(pi).nChannels;
    clear CData
    Cdata(ntime) = struct; %#ok<SAGROW>
    for ti = 1:ntime
        prefix = sprintf(bare,pi-1,ti-1);
        segname = [prefix '_masks.mat'];
        load(fullfile(segDir,segname),'masks','cellData','bgmask');

        % Compute XY positions, Z, and nuclear area from mask indices
        % Need image size to convert linear indices to subscripts
        % Load one image to get dimensions
        img_for_size = squeeze(readStack([dataDir '/' sprintf(bare, pi-1, ti-1) '.tif']));
        img_size = size(img_for_size);
        if numel(img_size) == 2
            img_size = [img_size 1];  % Add z dimension if 2D
        end
        cellData = struct();
        ncells = length(masks);
        XY = zeros(ncells, 2);
        Z = zeros(ncells, 1);
        nucArea = zeros(ncells, 1);
        for ci_mask = 1:ncells
            nucmask_idx = masks(ci_mask).nucmask;
            if ~isempty(nucmask_idx)
                % Convert linear indices to subscripts [row, col, z]
                [row, col, z] = ind2sub(img_size, nucmask_idx);
                % XY is [col, row] to match image coordinates (x=col, y=row)
                XY(ci_mask, 1) = mean(col);  % X coordinate (column)
                XY(ci_mask, 2) = mean(row);  % Y coordinate (row)
                Z(ci_mask) = mean(z);        % Z coordinate
                nucArea(ci_mask) = length(nucmask_idx);  % Volume in voxels
            end
        end
        cellData.XY = XY;
        cellData.Z = Z;
        cellData.nucArea = nucArea;

        fields = fieldnames(cellData);
        for fi = 1:length(fields)
            CData(ti).(fields{fi}) = cellData.(fields{fi}); %#ok<SAGROW>
        end
        ncells = length(masks);
        positions(pi).ncells(ti) = ncells;
        CData(ti).nucLevel = NaN(ncells, nchannels); 
        CData(ti).cytLevel = NaN(ncells, nchannels);
        CData(ti).NCratio = NaN(ncells, nchannels); 
        CData(ti).background = NaN(1,nchannels);
        for ci = 1:nchannels
            img = positions(pi).loadImage(dataDir,ci-1,ti);

            if useMIP && ndims(img) == 3
                % MIP across z
                img_mip = max(img, [], 3);

                % Convert 3D masks to 2D
                masks_2d = masks;
                for mi = 1:length(masks)
                    % Convert 3D linear indices to 2D
                    [row, col, ~] = ind2sub(img_size, masks(mi).nucmask);
                    masks_2d(mi).nucmask = unique(sub2ind(img_size(1:2), row, col));
                    [row, col, ~] = ind2sub(img_size, masks(mi).cytmask);
                    masks_2d(mi).cytmask = unique(sub2ind(img_size(1:2), row, col));
                end

                % Convert 3D bgmask to 2D (background if background in all z)
                bgmask_2d = all(bgmask, 3);

                [nL,cL,ncR,bg] = readIntensityValues(img_mip, masks_2d, bgmask_2d);
            else
                [nL,cL,ncR,bg] = readIntensityValues(img,masks,bgmask);
            end
            CData(ti).nucLevel(:,ci) = nL;
            CData(ti).cytLevel(:,ci) = cL;
            CData(ti).NCratio(:,ci) = ncR;
            CData(ti).background(ci) = bg;
        end
        fprintf('.')
        if mod(ti,50) == 0
            fprintf('\n')
        end
    end
    fprintf('\n')
    positions(pi).cellData = CData;

    if strcmp(imageType,'MP')
        % margin around nominal radius to exclude cells from, in micron
        % 96h colonies extend beyond 350 um
        margin = 50; 
        positions(pi).setCenter(margin);
    end
end
toc

save(fullfile(dataDir,'positions.mat'),'positions')
