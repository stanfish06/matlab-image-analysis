% REMOVE_BLACK_MARGIN  Crop black borders from stitched images and re-save [template].
%   Loads meta.mat and, for each position, removes the zero-valued border left by
%   stitching and writes the trimmed image to a 'processed' subfolder. Run from
%   the data directory.
clear; close all;
%% setup
open remove_black_margin.m
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
dataDir = scriptPath;
load(fullfile(dataDir,'meta.mat'));
bare_file = meta.filenameFormat;
npos =  meta.nPositions;
nt = 1;
nuc_ch = meta.nucChannel;
xsize = meta.xSize;
ysize = meta.ySize;
%% save each position as tiff
if ~exist(fullfile(dataDir, "processed"), 'dir')
    mkdir(fullfile(dataDir, "processed"))
end
for i = 1:npos
    boundary_mask = zeros(ysize, xsize);
    img_stack = repmat({zeros(ysize, xsize)}, 1, nt);
    for t = 1:nt
        fname = sprintf(bare_file, i-1, nuc_ch, t-1);
        img = imread(fullfile(dataDir, fname));
        img_stack{t} = img;
        boundary_mask = boundary_mask | (img == 0);
    end
    bbox = LargestRectangle(~boundary_mask, 5, 1, 0, 89.9999, 0);
    left_idx = max(ceil(max(bbox(2,1),bbox(5,1))), 1);
    right_idx = min(floor(min(bbox(3,1),bbox(4,1))), xsize);
    top_idx = max(ceil(max(bbox(2,2),bbox(3,2))), 1);
    bot_idx = min(floor(min(bbox(4,2),bbox(5,2))), ysize);
    if nt > 1
        img_stack_new = uint16(zeros(bot_idx-top_idx+1, right_idx-left_idx+1, nt));
        for t = 1:nt
            img_stack_new(:,:,t) = img_stack{t}(top_idx:bot_idx, left_idx:right_idx);
        end
        options.overwrite = true;
        saveastiff(img_stack_new, regexprep(fname, "_t\d{4}", ''), options);
    else
        img_stack_new = img_stack{1}(top_idx:bot_idx, left_idx:right_idx);
        imwrite(squeeze(img_stack_new), regexprep(fname, "_t\d{4}", ''));
    end
    movefile(fullfile(dataDir, regexprep(fname, "_t\d{4}", '')), fullfile(dataDir, "processed", regexprep(fname, "_t\d{4}", '')));
end
