% PREPROCESS_IMAGES_NONOVERLAP  Preprocessing for non-overlapping (unstitched) data [template].
%   Like PREPROCESS_IMAGES but for data acquired without montage overlap (no
%   stitching): sets metadata (nucChannel, conditions, channelLabel,
%   posPerCondition) manually and runs MIP/z-stack preprocessing. Edit the
%   metadata block near the top.
clear;
%% setup
scriptPath = pwd;
dataDir = scriptPath;

% -----------metadata------------------------
%nuclear marker channel (e.g., DAPI, HOECHST; indexed from 0)
nucChannel = 1;
%description of conditions in order (should be alphabetical order that the
%system finds image files but usually that is the same as the order of
%imaging)
conditions = {'BMP100','BMP100-LDN250', 'BMP100-LDN100', 'BMP100-LDN50'};
%names of channels, e.g. channelLabel = {'DAPI','SMAD23','pAKT','SOX17'};
channelLabel = {'SOX2','H2B'};
%if the number of positions per condition varies between conditions, set it
%manually; e.g., posPerCondition = [4 4 4 4 5]. If left empty ([]) then
%posPerCondition divides the number of positions by the number of
%conditions and assumes it is the same for each condition
posPerCondition = [4 4 4 4];

manualMeta = struct();
manualMeta.nucChannel = nucChannel;
manualMeta.channelLabel = channelLabel;
manualMeta.conditions = conditions;
manualMeta.posPerCondition = posPerCondition;

meta = Metadata(dataDir, manualMeta);
save(fullfile(dataDir,'meta.mat'),'meta');
%%
barefname = 'Image_???.vsi';
filelist = cellfun(@(x){fullfile(dataDir,x)}, sort({dir(fullfile(dataDir,barefname)).name}));
load('correctionMats_olympus_20x_lwd.mat','Gps','D')
ffchannels = {'GFP','RFP','Cy5'};
for f = filelist
    [~, baseName, ~] = fileparts(f{1});
    img = readStack(f{1});
    for t = 1:size(img, 5)
        img_slice = squeeze(img(:,:,:,:,t));
        for i = 1:length(meta.channelNames)
            ch_name = meta.channelNames{i};
            G = Gps{strcmp(ffchannels',ch_name)};
            mG = mean(G(:));
            G = (G - mG) * 0.2 + mG;
            img_slice_corrected = flatFieldCorrection({squeeze(img_slice(:,:,i))},G,D);
            img_slice(:,:,i) = img_slice_corrected{1};
        end
        img(:,:,:,:,t) = img_slice;
    end
    ImageJ;
    img_fiji = copytoImagePlus(img, 'XYCZT');
    img_fiji.show();
    ij.IJ.saveAs("tiff", [baseName '.tif']);
    img_fiji.close();
end
ij.IJ.run("Quit","");