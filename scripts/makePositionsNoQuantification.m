% MAKEPOSITIONSNOQUANTIFICATION  Create Position/Colony objects without quantifying [template].
%   Builds the Position (or Colony, for micropatterns) array and geometry from
%   meta.mat without reading out intensities, saving positions.mat. Useful when
%   only colony geometry/segmentation is needed. Edit imageType/radii near the top.
clear; close all; clc

scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
dataDir = scriptPath;

%% setup
load(fullfile(dataDir,'meta.mat'),'meta') %load metadata

% define imageType as micropatterned 'MP' or 'disordered' - for this
% application, not sure if we really need to include disordered
imageType = 'MP';
% radii of micropatterns in micron for each condition (ignored if not 'MP')
radii = 700*ones(1,meta.nWells)/2; 

npos = meta.nPositions;
filenames = meta.fileNames;
bare = meta.filenameFormat;

if isempty(filenames)
    stitchflag = true;
else
    stitchflag = false;
end

if ~stitchflag && (length(filenames) == 1) && npos > 1
    [~,~,ext] = fileparts(filenames{1});
    if strcmp(ext,'.nd2') || strcmp(ext,'.lif')
        filenames = repmat(filenames,1,npos);
    end
end

nucChannel = meta.nucChannel;

%% initialize positions
if strcmp(imageType,'MP')
    positions(meta.nPositions) = Colony();
elseif strcmp(imageType,'disordered')
    positions(meta.nPositions) = Position();
end

conditionStartPos = meta.conditionStartPos;

for pi = 1:meta.nPositions

    condi = find(pi >= conditionStartPos,1,'last'); % condition index
          
    fprintf('colony %d of %d\n',pi,npos)

    if strcmp(imageType,'MP')
        positions(pi) = Colony(meta, pi);
        positions(pi).setRadius(radii(condi), meta.xres);
        positions(pi).well = condi;
    elseif strcmp(imageType,'Disordered')
        positions(pi) = Position(meta, pi);
    end

    % shared properties between micropattern and disordered
    positions(pi).dataChannels = 1:positions(pi).nChannels;
    if stitchflag
        positions(pi).filenameFormat = bare;
    else
        positions(pi).filenameFormat = filenames{pi};
    end
end

save(fullfile(dataDir,'positions.mat'),'positions')

%% set centers
close all

r = round(10/meta.xres);
SE = strel("disk",r);
figure
for pi = 1:length(positions)
    condi = find(pi >= conditionStartPos,1,'last'); % condition index
    
    clf
    img = positions(pi).loadImage(dataDir,nucChannel,1);
    im = max(img,[],3);
    
    %functions for image binarization don't work well when we have
    %zero-padding so fill in zero regions from non-zero regions before
    %doing adaptive thresholding
    idxs = find(im==0); nzidxs = find(im>0); test = im;
    test(idxs) = test(nzidxs(1:length(idxs)));
    %binarize the image and clean up very small foreground and background
    %islands
    seg = imbinarize(im,adaptthresh(test));
    seg = imopen(imclose(seg,strel('disk',1)),strel('disk',3));
    %image close with a structure element with a radius of ~10 um
    J = imclose(seg,SE);
    
    %find the centroid of the largest object in the closed segmentation
    props = regionprops(J);
    [~,I] = max([props.Area]);
    center = props(I).Centroid;
    
    %show the image with the centroid and a circle outlining the colony
    imshow(imadjust(im,stitchedlim(im)))
    hold on
    scatter(center(1),center(2),100,'filled')
    viscircles(center,radii(condi)/meta.xres);
    hold off
    title(sprintf('colony %d',pi))
    cleanSubplot
    drawnow
    
    positions(pi).center = center;
end

save(fullfile(dataDir,'positions.mat'),'positions')


















