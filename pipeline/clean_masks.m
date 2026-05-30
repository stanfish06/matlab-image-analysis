function clean_masks(baseDir)
%CLEAN_MASKS Filter Cellpose masks and export MATLAB segmentation masks.

if nargin < 1 || isempty(baseDir)
    baseDir = pwd;
end
originalDir = pwd;
cleanupObj = onCleanup(@() cd(originalDir));
cd(baseDir);

%%
load("./meta.mat")
xres = meta.xres;
yres = meta.yres;
zres = 1.87 / 2;
npos = meta.nPositions;

nplanes_resize = 44; % make sure it matches the number of z splices (rescaled) when you are doing segmentation in cellpose, check cellpose output
nz = meta.nZslices;
rescale_ratio = nplanes_resize / nz;
nuc_ch = meta.nucChannel;
% Per-position filenames are derived from meta via imageFilePrefix so this runs
% on stitched AND non-stitched (writeZstacks) data without editing hard-coded
% 'stitched_' names. Images are read at time point 0.
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

% diagnostics: preallocate arrays for per-position stats
diag_n_total = zeros(1, npos);
diag_n_removed_small = zeros(1, npos);
diag_n_removed_dark = zeros(1, npos);
diag_n_removed_bright = zeros(1, npos);
diag_n_kept = zeros(1, npos);
diag_n_split = zeros(1, npos);

%%
dataDir = pwd;
logmsg('clean_masks: %d positions (min_volume=%d, bg_thresh=%d, abs_vol=%d, split_large=%d)', ...
    npos, min_volume, background_intensity, do_abs_vol_thresh, split_large);
parfor p = 1:npos
    % position index (0-based) used for the debug-overlay filename label
    pos_idx = p - 1;
    logmsg('clean_masks: position %d/%d start', p, npos);
    prefix = imageFilePrefix(meta, p, nuc_ch, 0, dataDir);
    lab_img = tiffreadVolume([prefix '_cp_masks.tif']);
    qc_img = tiffreadVolume([prefix '.tif']);
    lab_img(qc_img < background_intensity) = 0;
    lab_img = imclose(lab_img,[1 1 1;1 1 1;1 1 1]);
    % TODO: imclose may still lead some broken pieces. You may need to keep only the largest cc

    sz = size(lab_img);
    lab_img_resolved = zeros(sz);
    rp = regionprops3(lab_img, "Volume", "BoundingBox", "VoxelIdxList");
    v_mean = mean(rp.Volume);
    next_id = 1;

    % per-position counters
    cnt_total = height(rp);
    cnt_small = 0;
    cnt_dark = 0;
    cnt_bright = 0;
    cnt_kept = 0;
    cnt_split = 0;

    for i = 1:height(rp)
        if do_abs_vol_thresh
            v_ratio = rp(i,:).Volume;
        else
            v_ratio = rp(i,:).Volume / v_mean;
        end

        if (v_ratio < 0.1 && ~do_abs_vol_thresh) || (v_ratio < min_volume && do_abs_vol_thresh)
            % too small, remove
            lab_img_resolved(rp(i,:).VoxelIdxList{1}) = 0;
            cnt_small = cnt_small + 1;
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
                cnt_kept = cnt_kept + 1;
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
                cnt_split = cnt_split + 1;
                cnt_kept = cnt_kept + nseg;
            end
        else
            mean_qc_intense = mean(qc_img(rp(i,:).VoxelIdxList{1}));
            if mean_qc_intense > junk_intensity_mean || mean_qc_intense < background_intensity_mean
                lab_img_resolved(rp(i,:).VoxelIdxList{1}) = 0;
                if mean_qc_intense < background_intensity_mean
                    cnt_dark = cnt_dark + 1;
                else
                    cnt_bright = cnt_bright + 1;
                end
            else
                lab_img_resolved(rp(i,:).VoxelIdxList{1}) = next_id;
                next_id = next_id + 1;
                cnt_kept = cnt_kept + 1;
            end
        end
    end
    fname_save = [prefix '_FinalSegmentation.tif'];
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

    logmsg('clean_masks: position %d/%d done (%d/%d masks kept, %d small, %d dark, %d bright)', ...
        p, npos, cnt_kept, cnt_total, cnt_small, cnt_dark, cnt_bright);

    % store diagnostics
    diag_n_total(p) = cnt_total;
    diag_n_removed_small(p) = cnt_small;
    diag_n_removed_dark(p) = cnt_dark;
    diag_n_removed_bright(p) = cnt_bright;
    diag_n_kept(p) = cnt_kept;
    diag_n_split(p) = cnt_split;
end

%% print cleanup diagnostics
fprintf('\n--- Mask Cleanup Summary ---\n');
fprintf('%-6s %8s %8s %8s %8s %8s %8s %8s\n', ...
    'Pos', 'Total', 'Small', 'Dark', 'Bright', 'Split', 'Kept', 'Kept%');
for p = 1:npos
    if diag_n_total(p) > 0
        pct = 100 * diag_n_kept(p) / diag_n_total(p);
    else
        pct = 0;
    end
    fprintf('p%.4d %8d %8d %8d %8d %8d %8d %7.1f%%\n', ...
        p-1, diag_n_total(p), diag_n_removed_small(p), ...
        diag_n_removed_dark(p), diag_n_removed_bright(p), ...
        diag_n_split(p), diag_n_kept(p), pct);
end
fprintf('%-6s %8d %8d %8d %8d %8d %8d %7.1f%%\n', ...
    'TOTAL', sum(diag_n_total), sum(diag_n_removed_small), ...
    sum(diag_n_removed_dark), sum(diag_n_removed_bright), ...
    sum(diag_n_split), sum(diag_n_kept), ...
    100 * sum(diag_n_kept) / max(sum(diag_n_total), 1));

% save diagnostics as CSV
diag_table = table((0:npos-1)', diag_n_total', diag_n_removed_small', ...
    diag_n_removed_dark', diag_n_removed_bright', diag_n_split', diag_n_kept', ...
    'VariableNames', {'position', 'total', 'removed_small', 'removed_dark', ...
    'removed_bright', 'split', 'kept'});
writetable(diag_table, 'cleanup_stats.csv');
fprintf('Diagnostics saved to cleanup_stats.csv\n');

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

allpmasks = cell(1, npos);
allbgmasks = cell(1, npos);

parfor pidx = 1:npos
    prefix = barenames{pidx};
    logmsg('clean_masks: saving masks for position %d/%d (%s)', pidx, npos, prefix);
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

end
