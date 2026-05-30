% EXTRACT_COLONIES  Split quantified positions into individual colonies [template].
%   Loads meta.mat + positions.mat, detects colonies within each position, and
%   writes a per-colony Position array. Template variant using 'Image_%.3d_w%.4d'
%   file naming; edit the file-format strings and colony-picking parameters near
%   the top. Run from the data directory.
clear; close all; clc

scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
baseDir = scriptPath;
%% load in meta and position
meta = load(fullfile(baseDir, 'meta.mat')).meta;
positions = load(fullfile(baseDir, 'positions.mat')).positions;
% file format for the mask and img
bare_mask = "Image_%.3d_w%.4d_t0000_masks.mat";
bare_img = "Image_%.3d_w%.4d_t0000.tif";
npos = meta.nPositions;
% load mask
disp('loading nuclear and cytoplasmic masks')
tic
for i = 1:npos
    pos_idx = i; % note: if position index starts from 0, change i to i-1
    segname = sprintf(bare_mask,pos_idx,meta.nucChannel);
    imgname = sprintf(bare_img,pos_idx,meta.nucChannel);

    mask = load(fullfile(baseDir,segname));
    img_zstack = readStack(fullfile(baseDir, imgname));
    masks{i} = mask.masks;
    imgs{i} = img_zstack;
    cellDatas{i} = mask.cellData;
    bgmasks{i} = mask.bgmask;
end
toc
%% extract colonies
% number of colonies to extract per position
n_top_cc = 7;
% minimum area (in pixel) of each colony
min_area = 1000;
% max area (in pixel) of each colony
max_area = 60000;
% max distance (in pixel) to the center of the image
% this is useful because colonies near the boundary are often incomplete
max_dist = 1200;
% if dont want to exclude boundary colonies, set this to false
threshold_by_distance = true;
% this is an important parameter as it controls how much the mask will be
% dilated. If it is too high, two nearby colonies might be merged, and they
% will be assigned to a single colony. If it is too low, a colony might be
% incomplete.
dilate_radius = 10;

% this stores the montage images for the result of colony picking
% one can use this to check if colony picking is fine or not
colony_img_montage_per_colony_z = {};
colony_img_seg_montage_per_position_mip = {};
colony_img_montage_per_colony_mip = {};
% mask overlay transparency
mask_transparency = 0.25;

% this stores the updated position object per colony
positions_by_colony = {};
nPositions = 0;
nPositions_per_conditions = meta.posPerCondition * 0;
nPositions_per_conditions_original = meta.posPerCondition;
pos_count = 1;
condi_count = 1;
for i = 1:npos
    fprintf(['start position: ', num2str(i), '\n']);
    nuc_mask = {masks{i}.nucmask};
    cyto_mask = {masks{i}.cytmask};
    img = imgs{i};
    % if stitched, images from different positions might have different sizes, so using the image size stored in meta
    % can sometimes be wrong
    xsize = size(img, 2);
    ysize = size(img, 1);
    zsize = size(img, 5);
    img_center = [ysize/2 xsize/2];
    position = positions(i);
    % reconstruct binary nuclear mask
    % this will be used to find colonies
    % In brief words, we first get mip of the binary nuclear mask because colonies only differ in terms of their XY positiosn.
    % Then, we dilate the binary mip nuclear mask in order to connect cells
    % belonging to the same colony. Finally, we find the connected
    % componenets in the dilated binary mask.
    nuc_mask_bin = zeros([ysize xsize zsize]);
    for j = 1:length(nuc_mask)
        nuc_mask_bin(nuc_mask{j}) = 1;
    end
    % combine masks aross z
    nuc_mask_bin = max(nuc_mask_bin,[],3)>0;
    % dilate mask
    nuc_mask_bin = imdilate(nuc_mask_bin,strel('disk',dilate_radius));
    % fine connected components
    rc = regionprops(nuc_mask_bin, 'Area', 'Centroid', 'boundingbox');
    cc_centers = cell2mat({rc.Centroid}');
    cc_area = cell2mat({rc.Area});
    cc_dist = sqrt(sum((cc_centers - img_center).^2,2))';
    % remove colonies that are too small or near the boundary
    if threshold_by_distance
        rc = rc(cc_area > min_area & cc_area < max_area & cc_dist < max_dist);
    else
        rc = rc(cc_area > min_area & cc_area < max_area);
    end
    % get bounding boxes of the remaining colony
    cc_bbox = cell2mat({rc.BoundingBox}');
    cc_area = cell2mat({rc.Area});
    % sort by area
    [~,idx] = sort(cc_area, 'descend');
    % keep top n largest colonies
    nseg = min(n_top_cc,height(rc));
    idx = idx(1:nseg);
    % get bounding boxes of the top n colonies
    bbox = cc_bbox(idx,:);

    colony_mask = zeros([ysize xsize zsize]);
    colony_img_per_colony_z = {};
    colony_img_per_colony_mip = {};
    for j = 1:nseg
        x = bbox(j,1);
        y = bbox(j,2);
        w = bbox(j,3);
        h = bbox(j,4);
        colony_mask(ceil(y):floor(y+h), ceil(x):floor(x+w), :) = j;
    end

    colony_lab = zeros([ysize xsize zsize]);
    cell_lab = zeros([ysize xsize zsize]);
    next_cell_idx = 1;
    % position objects of colonies from the same fov
    positions_cc = {};
    % clean mask
    for j = 1:nseg
        position_cc = Position();
        position_old = struct(position);
        fieldNames = fieldnames(position_old);
        % copy original data to a new position object
        for k = 1:numel(fieldNames)
            try
                position_cc.(fieldNames{k}) = position_old.(fieldNames{k});
            catch
            end
        end

        colony_cell_idx = [];
        for k = 1:length(nuc_mask)
            colony_idx = max(colony_mask(nuc_mask{k}));
            if colony_idx == j
                if isempty(cyto_mask{k}) || isempty(nuc_mask{k})
                    continue
                end
                cell_lab(nuc_mask{k}) = next_cell_idx;
                next_cell_idx = next_cell_idx + 1;
                colony_lab(nuc_mask{k}) = j;
                colony_cell_idx(end + 1) = k;
            end
        end

        ncells_original = position_cc.ncells;
        position_cc.ncells = length(colony_cell_idx);
        cellData_colony = position_cc.cellData;
        % subset cellData
        fields = fieldnames(cellData_colony);
        for k = 1:length(fields)
            fieldName = fields{k};
            fieldValue = cellData_colony.(fieldName);
            % Check if the field is numeric and has the matching number of rows
            if isnumeric(fieldValue) && size(fieldValue, 1) == ncells_original
                % Subset the field using the indexArray
                cellData_colony.(fieldName) = fieldValue(colony_cell_idx, :);
            end
        end
        position_cc.cellData = cellData_colony;
        positions_cc{end + 1} = position_cc;
    end

    colony_lab = max(uint16(colony_lab), [], 3);
    img_display = imadjust(max(img, [], 5));
    colony_img_seg_montage_per_position_mip{end + 1} = visualize_nuclei_v2(colony_lab,img_display,mask_transparency);

    overlay_montage_per_position = {};
    for j = 1:nseg
        x = bbox(j,1);
        y = bbox(j,2);
        w = bbox(j,3);
        h = bbox(j,4);
        colony_img_per_colony_mip = [colony_img_per_colony_mip {colony_img_seg_montage_per_position_mip{end}(ceil(y):floor(y+h), ceil(x):floor(x+w), :)}];
        img_stack_colony = squeeze(img(ceil(y):floor(y+h), ceil(x):floor(x+w), :, :, :));
        img_stack_colony = imadjustn(img_stack_colony);
        overlay_montage = {colony_img_per_colony_mip{j}};
        for k = 1:zsize
            img_slice = img_stack_colony(:,:,k);
            cell_lab_colony = cell_lab(ceil(y):floor(y+h), ceil(x):floor(x+w), k);
            mask_overlay = visualize_nuclei_v2(cell_lab_colony>0,img_slice,mask_transparency);
            overlay_montage{end + 1} = mask_overlay;
        end
        overlay_montage_per_position = [overlay_montage_per_position overlay_montage];
    end
    colony_img_per_colony_z = [colony_img_per_colony_z {imtile(overlay_montage_per_position,"GridSize",[nseg zsize+1])}];
    colony_img_montage_per_colony_mip = [colony_img_montage_per_colony_mip {colony_img_per_colony_mip}];
    colony_img_montage_per_colony_z = [colony_img_montage_per_colony_z {colony_img_per_colony_z}];

    positions_by_colony = [positions_by_colony positions_cc];
    nPositions = nPositions + nseg;
    if pos_count <= nPositions_per_conditions_original(condi_count)
        pos_count = pos_count + 1;
        nPositions_per_conditions(condi_count) = nPositions_per_conditions(condi_count) + nseg;
    else
        condi_count = condi_count + 1;
        pos_count = 2;
        nPositions_per_conditions(condi_count) = nPositions_per_conditions(condi_count) + nseg;
    end
end
%% sanity check for colony picking
if ~exist(fullfile(baseDir, 'colony_picking'), 'dir')
    mkdir(fullfile(baseDir, 'colony_picking'));
end
condi_ct = 1;
ct = 1;
height_target = 5000;
for i = 1:npos
    if ct > nPositions_per_conditions_original(condi_ct)
        condi_ct = condi_ct + 1;
        ct = 2;
        pos_idx = 1;
    else
        pos_idx = ct;
        ct = ct + 1;
    end
    condi = meta.conditions(condi_ct);
    f = figure;
    img1 = imadjust(squeeze(max(imgs{i}, [], 5)));
    img2 = colony_img_seg_montage_per_position_mip{i};
    img3 = colony_img_montage_per_colony_z{i}{1};
    sz1 = size(img1);
    sz2 = size(img2);
    sz3 = size(img3);
    % make sure images have the same height
    img1 = imresize(img1, [height_target height_target / sz1(1) * sz1(2)]);
    img2 = imresize(img2, [height_target height_target / sz2(1) * sz2(2)]);
    img3 = imresize(img3, [height_target height_target / sz3(1) * sz3(2)]);
    montage({img1, img2, img3}, 'Size', [1 3], 'BorderSize', 5, 'BackgroundColor', 'w');
    title([condi{1} ': ' sprintf('pos\\_%.2d', pos_idx)], "FontSize", 15, "FontWeight", 'bold');
    exportgraphics(gca, fullfile(baseDir, 'colony_picking', [condi{1} '_' sprintf('position_%.4d.png', i)]), 'Resolution', 300);
    close(f);
end
%% if things look ok, update meta and save new meta and positions
meta.nPositions = nPositions;
meta.posPerCondition = nPositions_per_conditions;
positions = positions_by_colony;
save("positions_colony.mat", 'positions');
save("meta_colony.mat", "meta");