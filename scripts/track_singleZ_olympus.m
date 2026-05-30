% TRACK_SINGLEZ_OLYMPUS  Single-z live tracking pipeline for Olympus data [template].
%   Loads/patches metadata (Olympus resolution and time must be set manually),
%   sets segmentation + cleanup options, then segments, extracts per-frame
%   cellData (extractDataSingleFrame_stan) and tracks cells. Edit meta.xres/yres,
%   nTime and the opts blocks near the top.
clear; close all;
%% setup
open track_MIP.m
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
load("meta.mat")
% not sure what is going on but olympus resolution cannot be determined automatically
% also you need to specify time
meta.xres = 1.3;
meta.yres = 1.3;
meta.nTime = 1;
save(fullfile(scriptPath, 'meta_update'), 'meta');
clear meta;
%% opts
load('./meta_update.mat');
opts = struct(...
    'cytoSize',             4,...
    'cytoMargin',           -1,...
    'nucShrinkage',         2,...
    'cytoplasmicLevels',    false);

opts.cleanupOptions = struct('separateFused', false,...
    'clearBorder',true,...
    'minAreaStd', 1,...
    'minSolidity',0.95,...
    'minArea', 0,...
    'openSize', 1,...
    'fillholes', true);

% dividing parent and early daughter cells are quite small
minarea = 0/(meta.xres^2); %area threshold for getting rid of junk
maxarea = 300/(meta.xres^2); %area threshold for discarding fused clumps of nuclei
% convex decomposition options
% when nuclei are not perfectly convex, set tau1 and tau2 lower
tau1 = 0.6; %relative concavity threshold (dimensionless)
tau2 = 0.4/meta.xres; %absolute concavity threshold (in pixels)
tau3 = 1.9/meta.xres; %larger, overriding absolute concavity threshold

opts.decompopts = struct(...
    'flag',         false,...
    'tau3',        tau3,...
    'useMinArea',   true,...
    'minArea',      minarea,...
    'tau1',         tau1,...
    'tau2',         tau2,...
    'maxArea',      maxarea,...
    'ignoreholes',  false);
%% clean up segmentation
if ~exist('opts', 'var')
    error('The variable "opts" does not exist. Please run the opts block first.');
end
% 1 represents regular cell
% 2 represents daughter cell
% 3 represents dividing parent cell
% 4 represents junk
fgChannels = [1, 2, 3];
divChannels = [2, 3];

overlaydir = fullfile(scriptPath, 'SegOverlays');
if ~exist(overlaydir,'dir')
    mkdir(overlaydir);
end
list = dir(fullfile(scriptPath, "processed", "/*Object Predictions.h5"));
seg_names  = fullfile(scriptPath, "processed", {list.name});
list = dir(fullfile(scriptPath, "processed", ['\*', sprintf('_w%.4d.tif', meta.nucChannel)]));
img_names  = fullfile(scriptPath, "processed", {list.name});
npos = meta.nPositions;
%%
cmap_ilastik = [1, 1, 0; 0, 0, 1; 1, 0, 0; 0, 1, 1;];
nucChannel = meta.nucChannel;
pattern = 'clean_seg_p%.4d_w%.4d.tif';
for ii = 1:npos
    fprintf('position %d of %d\n',ii,npos)
    fname = seg_names{ii};
    allseg = ilastikRead(fname);
    nt = size(allseg, 3);
    fname_img = img_names{ii};
    newsegs = zeros(size(allseg),'uint8');
    parfor ti = 1:nt
        fprintf('.');
        img = imread(fname_img, ti);
        I = imadjust(img,stretchlim(img),[],0.8);
        seg = allseg(:, :, ti);
        B = labeloverlay(I, seg, 'Colormap', cmap_ilastik, 'Transparency', 0.8);
        [~,name,~] = fileparts(fname_img);
        savename = [name, sprintf('_t%.4d',ti-1), '_SegOverlay_Object_classification.jpg'];
        savename = fullfile(overlaydir,savename);
        savename_img = [name, sprintf('_t%.4d',ti-1), '.png'];
        savename_img = fullfile(overlaydir,savename_img);
        imwrite(B,savename);
        imwrite(I,savename_img);
        tic
        % further clean up the segmentation
        newseg = newNuclearCleanup(seg, fgChannels, divChannels, opts);
        newseg = uint8(newseg);
        newsegs(:,:,ti) = newseg;
        toc
        B = visualize_nuclei_v2(newseg>0,I);
        savename = [name, sprintf('_t%.4d',ti-1), '_SegOverlay_cleaned.jpg'];
        savename = fullfile(overlaydir,savename);
        imwrite(B,savename);
    end
    writename = fullfile(scriptPath,"processed",sprintf(pattern, ii-1, nucChannel));
    t = Tiff(writename, 'w8');
    tagstruct.ImageLength = size(newsegs, 1);
    tagstruct.ImageWidth = size(newsegs, 2);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 8;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Compression = Tiff.Compression.None;
    for ti = 1:nt
        t.setTag(tagstruct);
        t.write(newsegs(:,:,ti));
        if ti < meta.nTime
            t.writeDirectory();
        end
    end
    t.close();
end
%% create position object for the final segmentation
if ~exist('opts', 'var')
    error('The variable "opts" does not exist. Please run the opts block first.');
end
positions(meta.nPositions) = Position;
clear cellData;
list = dir(fullfile(scriptPath, "processed", 'clean_seg_*.tif'));
names = {list.name};
opts.tMax = meta.nTime;
opts.dataChannels = 0;
list = dir(fullfile(scriptPath, "processed", ['stitched', '*', sprintf('_w%.4d.tif', meta.nucChannel)]));
names_nuc = {list.name};
segnames = fullfile(scriptPath, "processed", names);
for ci = 1:length(meta.conditionStartPos)
    for cpi = 1:meta.posPerCondition(ci)
        pidx = meta.conditionStartPos(ci)+cpi-1;
        fprintf('Processing position %d\n',pidx)
        segs = cell(1,1,meta.nTime);
        for ti = 1:meta.nTime
            segs{ti} = imread(segnames{pidx},ti) == 1;
        end
        segs = cell2mat(segs);
        opts.nuclearSegmentation = segs;
        positions(pidx) = Position(meta, pidx);
        % here just extract nuclear data, but could be adapted to multiple channels
        fname = fullfile(scriptPath, "processed", names_nuc{pidx});
        im_nuc = readStack(fname);
        nchannels = 1;
        for ti = 1:meta.nTime
            fprintf('Processing time point %d\n',ti)
            imgs = cell(1,nchannels);
            for cii = 1:nchannels
                %xyczt
                imgs{cii} = squeeze(im_nuc(:,:,opts.dataChannels(cii)+1,:,ti));
            end
            nucseg = segs(:,:,ti);
            [~, CData] = extractDataSingleFrame_stan(imgs,nucseg,opts);
            fields = fieldnames(CData);
            for fi = 1:length(fields)
                cellData(ti).(fields{fi}) = CData.(fields{fi});
            end
            if any(isnan(cellData(ti).background)) && ti > 1
                cellData(ti).background = cellData(ti-1).background;
            end
        end
        positions(pidx).cellData = cellData;
        positions(pidx).ncells = cellfun(@(x) size(x,1),{cellData.XY})';
        clear cellData
        save(fullfile(pwd, 'positions'), 'positions');
    end
end
