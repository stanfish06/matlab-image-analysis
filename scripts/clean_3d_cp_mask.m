%%
load("./meta.mat")
xres = meta.xres;
yres = meta.yres;
zres = 1.87 / 2;
npos = meta.nPositions;

nplanes_resize = 44; % make sure it matches the number of z splices (rescaled) when you are doing segmentation in cellpose, check cellpose output
nz = meta.nZslices;
rescale_ratio = nplanes_resize / nz;
% bare_fname = "stitched_p%.4d_w0000_t0000_cp_masks.tif";
% bare_fname_qc_img = "stitched_p%.4d_w0000_t0000.tif";
% bare_fname_save = "stitched_p%.4d_w0000_t0000_FinalSegmentation.tif";
bare_fname = "stitched_p%.4d_w0000_t0000_cp_masks.tif";
bare_fname_qc_img = "stitched_p%.4d_w0000_t0000.tif";
bare_fname_save = "stitched_p%.4d_w0000_t0000_FinalSegmentation.tif";
% background intensity for the qc image
% it is the easiest to just figure this out manually, like via fiji or something
% quantile based method will be affected by the number of positive objects in the image
% the purpose of this is to remove masks that segment super dark objects, which happens when you took many z slices and you only have a few number of cells in fov
% check the downstream scatter plots, if you see many points near origin.
% you should be able to estimate background intensity by taking complement of the cell mask
% TODO: can you filter out small pieces more efficiently?
background_intensity = 600; % pixel based filtering, useful to trim ambiant dark masks
junk_intensity = 5000; % not useful, using this will negativelly affect mean-based junk removal. For dead cells, it is more effective to remove as a whole
background_intensity_mean = 800; % should be greater then pixel bg thresh
junk_intensity_mean = 6000;
% when the sample is super sparse and dirth, relative volume threshold becomes the segmentation will contain large protion of junks. In that case, use absolute volume threshold.
do_abs_vol_thresh = true;
min_volume = 400;
split_large = false;
debug = true; % save MIP overlay of cleaned masks on QC image for visual inspection
debug_dir = fullfile(pwd, 'debug_overlay');
if debug && ~exist(debug_dir, 'dir')
    mkdir(debug_dir);
end

%%
parfor p = 1:npos
    % make sure to change this (e.g. p or p - 1) so that it matches the position indices of the files
    pos_idx = p - 1;
    lab_img = tiffreadVolume(sprintf(bare_fname, pos_idx));
    qc_img = tiffreadVolume(sprintf(bare_fname_qc_img, pos_idx));
    lab_img(qc_img < background_intensity) = 0;
    lab_img = imclose(lab_img,[1 1 1;1 1 1;1 1 1]);
    % TODO: imclose may still lead some broken pieces. You may need to keep only the largest cc

    sz = size(lab_img);
    lab_img_resolved = zeros(sz);
    rp = regionprops3(lab_img, "Volume", "BoundingBox", "VoxelIdxList");
    v_mean = mean(rp.Volume);
    next_id = 1;
    for i = 1:height(rp)
        if do_abs_vol_thresh
            v_ratio = rp(i,:).Volume;
        else
            v_ratio = rp(i,:).Volume / v_mean;
        end

        if (v_ratio < 0.1 && ~do_abs_vol_thresh) || (v_ratio < min_volume && do_abs_vol_thresh)
            % too small, remove
            lab_img_resolved(rp(i,:).VoxelIdxList{1}) = 0;
        elseif split_large && ((v_ratio > 1.5 && ~do_abs_vol_thresh) || (v_ratio > 3 * min_volume && do_abs_vol_thresh))
            % reduce this if you want to split more
            % do watershed to split undersegmented objects
            bbox = floor(rp(i,:).BoundingBox);
            bbox(bbox == 0) = 1;
            img_bbox = lab_img(bbox(2):min(bbox(2)+bbox(5),sz(1)),bbox(1):min(bbox(1)+bbox(4),sz(2)),bbox(3):min(bbox(3)+bbox(6),sz(3)));
            bw_bbox = img_bbox == i;

            % resize the bounding box to the original size
            % this is important as we use distance transform for watershed, so you need to make sure that z distance is not underestimated
            sz_bbox = size(bw_bbox);
            bw_bbox = imresize3(bw_bbox, [sz_bbox(1) sz_bbox(2) round(sz_bbox(3)*rescale_ratio)]);

            D = bwdist(~bw_bbox);
            regionalMaxima = imregionalmax(D);
            markers = bwlabeln(regionalMaxima);
            rp_bbox = regionprops3(markers, "VoxelIdxList");

            if height(rp_bbox) <= 1
                % we cannot split in this case
                lab_img_resolved(rp(i,:).VoxelIdxList{1}) = next_id;
                next_id = next_id + 1;
                continue;
            else
                % split as much as we can
                [~, idx] = sort(-cellfun(@(x)mean(D(x)),rp_bbox.VoxelIdxList));
                 % you may also manually increase this if you want to split more
                nseg = min(height(rp_bbox), ceil(v_ratio));
                keep_ids = idx(1:nseg);
                markers(~ismember(markers, keep_ids)) = 0;
                D = bwdist(markers ~= 0);
                D(~bw_bbox) = Inf;
                DL = uint16(watershed(D));
                DL(~bw_bbox) = 0;
                % down sample to the original size
                DL = DL(:,:,round(linspace(1, size(bw_bbox, 3), sz_bbox(3))));
                cell_idx = unique(DL(:))';
                cell_idx_count = 1;
                for k = cell_idx(cell_idx > 0)
                    DL(DL == k) = cell_idx_count;
                    cell_idx_count = cell_idx_count + 1;
                end

                bg_mask_bbox = DL == 0;
                DL = DL + next_id - 1;
                DL(bg_mask_bbox) = 0;
                nseg = sum(unique(DL(:)) > 0);
                bbox_resolved_current = lab_img_resolved(bbox(2):min(bbox(2)+bbox(5),sz(1)),bbox(1):min(bbox(1)+bbox(4),sz(2)),bbox(3):min(bbox(3)+bbox(6),sz(3)));
                DL(bg_mask_bbox) = bbox_resolved_current(bg_mask_bbox);
                lab_img_resolved(bbox(2):min(bbox(2)+bbox(5),sz(1)),bbox(1):min(bbox(1)+bbox(4),sz(2)),bbox(3):min(bbox(3)+bbox(6),sz(3))) = DL;
                next_id = next_id + nseg;
            end
        else
            mean_qc_intense = mean(qc_img(rp(i,:).VoxelIdxList{1}));
            if mean_qc_intense > junk_intensity_mean || mean_qc_intense < background_intensity_mean
                lab_img_resolved(rp(i,:).VoxelIdxList{1}) = 0   
            else
                lab_img_resolved(rp(i,:).VoxelIdxList{1}) = next_id;
                next_id = next_id + 1;
            end
        end
    end
    fname_save = sprintf(bare_fname_save, pos_idx);
    saveastiff(uint16(lab_img_resolved), char(fname_save),struct('overwrite',true));

    if debug
        qc_mip = max(qc_img, [], 3);
        qc_mip_norm = im2uint8(imadjust(mat2gray(qc_mip)));
        lab_mip = max(lab_img_resolved, [], 3);
        nlabels = max(lab_mip(:));
        ncolors = 64;
        base_cmap = hsv(ncolors);
        cmap = repmat(base_cmap, ceil(nlabels / ncolors), 1);
        cmap = cmap(1:max(nlabels,1), :);
        rng(0); cmap = cmap(randperm(size(cmap, 1)), :);
        overlay = labeloverlay(qc_mip_norm, lab_mip, 'Colormap', cmap, 'Transparency', 0.5);
        imwrite(overlay, fullfile(debug_dir, sprintf('debug_p%.4d.png', pos_idx)));
    end
end

%%
dataDir = pwd;
load(fullfile(dataDir,'meta.mat'));
nucChannel = meta.nucChannel;

filelist = dir(fullfile(dataDir,['*',sprintf('_w%.4d_t0000.tif',nucChannel)]));
npos = length(filelist);
names = {filelist.name};
barenames = cell(size(names));
for ii = 1:npos
    [~,barename,~] = fileparts(names{ii});
    barenames{ii} = barename;
end

allpmasks = {};
allbgmasks = {};
allcellData = {};

parfor pidx = 1:npos
    prefix = barenames{pidx};
    disp(prefix)
    % change this to cellpose as this code does not look good
    segname = fullfile(dataDir,[prefix,'_FinalSegmentation.tif']);
    cp_mask = tiffreadVolume(segname);
    rp = regionprops3(cp_mask, "VoxelIdxList");
    masks = struct('nucmask', cell(1, height(rp)), 'cytmask', cell(1, height(rp)));

    for i = 1:height(rp)
        masks(i).nucmask = rp.VoxelIdxList{i};
        masks(i).cytmask = rp.VoxelIdxList{i};
    end

    allpmasks{pidx} = masks;
    allbgmasks{pidx} = cp_mask == 0;
end

for pidx = 1:npos
    prefix = barenames{pidx};
    masks = allpmasks{pidx};
    bgmask = allbgmasks{pidx};
    save(fullfile(dataDir,[prefix,'_masks.mat']),'masks','bgmask')
end
