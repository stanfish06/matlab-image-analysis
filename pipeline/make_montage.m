function make_montage(baseDir)
%MAKE_MONTAGE Generate per-channel and merged image montage summaries.

if nargin < 1 || isempty(baseDir)
    baseDir = pwd;
end
originalDir = pwd;
cleanupObj = onCleanup(@() cd(originalDir));
cd(baseDir);

%% load metadata and determine number of positions automatically
load('meta.mat');
channels = meta.channelLabel;
npos = meta.nPositions;
nch = length(channels);

% per-position filenames derived from meta (imageFilePrefix): works for stitched
% and non-stitched data
img_chs = cell(1, nch); % img_chs{ch}{pos} = MIP for that channel/position
for ch = 1:nch
    img_chs{ch} = cell(1, npos);
end
bboxes = cell(1, npos);

all_ch_imgs = cell(1, npos); % temporary: all_ch_imgs{i}{ch} = cropped MIP
parfor i = 1:npos
    ch_imgs = cell(1, nch);
    for ch = 1:nch
        ch_imgs{ch} = max(tiffreadVolume([imageFilePrefix(meta, i, ch-1, 0, baseDir) '.tif']), [], 3);
    end
    % crop to content using first channel mask
    img_mask = ch_imgs{1} > 0;
    rp = regionprops(img_mask, 'BoundingBox');
    sz = size(ch_imgs{1});
    bbox = floor(rp(1).BoundingBox);
    bbox(bbox == 0) = 1;
    bboxes{i} = bbox;
    for ch = 1:nch
        ch_imgs{ch} = ch_imgs{ch}(bbox(2):min(bbox(2)+bbox(4),sz(1)), bbox(1):min(bbox(1)+bbox(3),sz(2)));
    end
    all_ch_imgs{i} = ch_imgs;
end
% unpack into img_chs{ch}{pos} layout
for i = 1:npos
    for ch = 1:nch
        img_chs{ch}{i} = all_ch_imgs{i}{ch};
    end
end
clear all_ch_imgs;

%% auto contrast: estimate per-channel bounds from percentiles across all positions
% set auto_contrast = false to use manual ch_bounds instead
auto_contrast = true;
% percentiles for contrast stretching (adjust if images look too bright/dark)
pct_low = 1;    % percentile for black point
pct_high = 99.5; % percentile for white point

if auto_contrast
    ch_bounds = zeros(nch, 2);
    for ch = 1:nch
        % sample pixels from all positions to estimate range
        all_pixels = [];
        for i = 1:npos
            img_ch = img_chs{ch}{i};
            % subsample to keep memory reasonable
            pixels = img_ch(img_ch > 0); % exclude zero-padded background
            if numel(pixels) > 50000
                pixels = pixels(randperm(numel(pixels), 50000));
            end
            all_pixels = [all_pixels; double(pixels(:))];
        end
        p = prctile(all_pixels, [pct_low pct_high]);
        max_val = double(intmax('uint16'));
        ch_bounds(ch, :) = p / max_val;
    end
    fprintf('Auto contrast bounds (normalized):\n');
    for ch = 1:nch
        fprintf('  %s: [%.4f, %.4f]\n', channels{ch}, ch_bounds(ch,1), ch_bounds(ch,2));
    end
else
    % manual contrast bounds (normalized to [0,1] of uint16 range)
    % example ENDO (DAPI,TFAP2C,EOMES,SOX17):
    % ch_bounds = [0, 0.005; 0.01, 0.07; 0.008, 0.035; 0.01, 0.06];
    ch_bounds = repmat([0, 0.1], nch, 1);
end

%%
montageDir = 'mip_montage';
if ~exist(montageDir, 'dir')
    mkdir(montageDir);
end

%%
% to set manually: colors = {[0.5 0.5 0.5], [1 0 0], [0 1 0], [0 0 1]}; scales = [1, 1, 0.9, 1.0];
default_colors = {[0.5 0.5 0.5], [1 0 0], [0 1 0], [0 0 1], [1 1 0], [1 0 1], [0 1 1]};
colors = default_colors(1:nch);
scales = ones(1, nch);
textbox = [0.1 0.18 0.8 0.1];
scatter_csv = fullfile('scatter', 'scatterPercent.csv');
plot_percent = exist(scatter_csv, 'file') == 2;
if plot_percent
    dat = readtable(scatter_csv);
else
    fprintf('scatter/scatterPercent.csv not found, skipping percent overlay\n');
end

%% Population overlay settings
do_population_overlay = false;
overlay_on_single_channels = true;

if exist('positions_colony.mat', 'file')
    position_file = 'positions_colony.mat';
else
    position_file = 'positions.mat';
end
if do_population_overlay
    load(position_file, 'positions');
end

% example ENDO: channelThresholds = [100, 200, 150, 250];
channelThresholds = repelem(100, nch);

% example ENDO populations:
% populations = {
%     'EOMES+SOX17+', [2, 4], [1, 1], [0 1 1], 'disk', 5;
%     'EOMES+SOX17-', [2, 4], [1, -1], [1 0 0], 'diamond', 5;
%     'TFAP2C+SOX17+', [3, 4], [1, 1], [1 0 0], 'square', 10;
% };
populations = {};

createMarkerKernel = @(shape, radius) getnhood(strel(shape, radius));

%%
for i = 1:length(meta.conditions)
    condi = meta.conditions{i};
    start_i = meta.conditionStartPos(i);
    if i < length(meta.conditions)
        end_i = meta.conditionStartPos(i + 1) - 1;
    else
        end_i = npos;
    end
    n_panels = end_i - start_i + 1;

    well_montage = cell(1, nch); % per-channel montage panels
    for ch = 1:nch
        well_montage{ch} = cell(1, n_panels);
    end
    well_montage_merged = cell(1, n_panels);

    for j = 1:n_panels
        pos_idx = start_i + j - 1;
        % stack all channels
        img = zeros([size(img_chs{1}{pos_idx}) nch], 'uint16');
        for ch = 1:nch
            img(:,:,ch) = img_chs{ch}{pos_idx};
        end

        % apply contrast adjustment per channel
        for k = 1:nch
            if k == 1
                img(:,:,k) = medfilt2(imadjust(img(:,:,k), stretchlim(img(:,:,k))'+ch_bounds(k,:), [0 1]), [3 3]);
            else
                img(:,:,k) = medfilt2(imadjust(img(:,:,k), ch_bounds(k,:), [0 1]), [3 3]);
            end
        end

        for ch = 1:nch
            well_montage{ch}{j} = img(:,:,ch);
        end

        % create merged RGB image
        well_montage_merged{j} = uint16(zeros([size(img, [1 2]) 3]));
        for ci = 1:3
            for k = 1:nch
                well_montage_merged{j}(:,:,ci) = well_montage_merged{j}(:,:,ci) + img(:,:,k) * scales(k) * colors{k}(ci);
            end
        end

        % Inject population markers directly into the merged image pixels
        if do_population_overlay && ~isempty(populations)
            bbox = bboxes{pos_idx};
            img_size = size(well_montage_merged{j});

            if iscell(positions)
                pos_data = positions{pos_idx};
            else
                pos_data = positions(pos_idx);
            end

            if ~isempty(pos_data.cellData) && isfield(pos_data.cellData, 'XY') && ~isempty(pos_data.cellData.XY)
                for pop_idx = 1:size(populations, 1)
                    pop_channels = populations{pop_idx, 2};
                    pop_signs = populations{pop_idx, 3};
                    pop_color = populations{pop_idx, 4};
                    pop_shape = populations{pop_idx, 5};
                    pop_radius = populations{pop_idx, 6};

                    marker_kernel = createMarkerKernel(pop_shape, pop_radius);
                    kernel_size = size(marker_kernel);
                    half_k = floor(kernel_size / 2);

                    nucLevel = pos_data.cellData.nucLevel;
                    bg = pos_data.cellData.background;
                    cell_idx = true(size(nucLevel, 1), 1);
                    for k_ch = 1:length(pop_channels)
                        ch = pop_channels(k_ch);
                        level = nucLevel(:, ch) - bg(ch);
                        if pop_signs(k_ch) > 0
                            cell_idx = cell_idx & (level > channelThresholds(ch));
                        else
                            cell_idx = cell_idx & (level <= channelThresholds(ch));
                        end
                    end

                    xy = pos_data.cellData.XY(cell_idx, :);
                    if ~isempty(xy) && size(xy, 1) > 0
                        xy_local = zeros(size(xy));
                        xy_local(:,1) = xy(:,1) - bbox(1) + 1;
                        xy_local(:,2) = xy(:,2) - bbox(2) + 1;
                        max_val = double(intmax('uint16'));
                        marker_color_uint16 = uint16(pop_color * max_val);

                        for cell_i = 1:size(xy_local, 1)
                            cx = round(xy_local(cell_i, 1));
                            cy = round(xy_local(cell_i, 2));
                            row_start = max(1, cy - half_k(1));
                            row_end = min(img_size(1), cy + half_k(1));
                            col_start = max(1, cx - half_k(2));
                            col_end = min(img_size(2), cx + half_k(2));
                            k_row_start = row_start - (cy - half_k(1)) + 1;
                            k_row_end = kernel_size(1) - ((cy + half_k(1)) - row_end);
                            k_col_start = col_start - (cx - half_k(2)) + 1;
                            k_col_end = kernel_size(2) - ((cx + half_k(2)) - col_end);
                            kernel_region = marker_kernel(k_row_start:k_row_end, k_col_start:k_col_end);
                            for rgb_ch = 1:3
                                img_region = well_montage_merged{j}(row_start:row_end, col_start:col_end, rgb_ch);
                                img_region(kernel_region) = marker_color_uint16(rgb_ch);
                                well_montage_merged{j}(row_start:row_end, col_start:col_end, rgb_ch) = img_region;
                            end
                        end
                    end
                end
            end
        end

        % Overlay markers on single channel images
        if do_population_overlay && overlay_on_single_channels && ~isempty(populations)
            for ch = 1:nch
                well_montage{ch}{j} = repmat(well_montage{ch}{j}, [1 1 3]);
            end
            single_ch_images = cell(1, nch);
            for ch = 1:nch
                single_ch_images{ch} = well_montage{ch}{j};
            end

            bbox = bboxes{pos_idx};
            img_size = size(single_ch_images{1});

            if iscell(positions)
                pos_data = positions{pos_idx};
            else
                pos_data = positions(pos_idx);
            end

            if ~isempty(pos_data.cellData) && isfield(pos_data.cellData, 'XY') && ~isempty(pos_data.cellData.XY)
                for pop_idx = 1:size(populations, 1)
                    pop_channels = populations{pop_idx, 2};
                    pop_signs = populations{pop_idx, 3};
                    pop_color = populations{pop_idx, 4};
                    pop_shape = populations{pop_idx, 5};
                    pop_radius = populations{pop_idx, 6};

                    marker_kernel = createMarkerKernel(pop_shape, pop_radius);
                    kernel_size = size(marker_kernel);
                    half_k = floor(kernel_size / 2);

                    nucLevel = pos_data.cellData.nucLevel;
                    bg = pos_data.cellData.background;
                    cell_idx = true(size(nucLevel, 1), 1);
                    for k_ch = 1:length(pop_channels)
                        ch = pop_channels(k_ch);
                        level = nucLevel(:, ch) - bg(ch);
                        if pop_signs(k_ch) > 0
                            cell_idx = cell_idx & (level > channelThresholds(ch));
                        else
                            cell_idx = cell_idx & (level <= channelThresholds(ch));
                        end
                    end

                    xy = pos_data.cellData.XY(cell_idx, :);
                    if ~isempty(xy) && size(xy, 1) > 0
                        xy_local = zeros(size(xy));
                        xy_local(:,1) = xy(:,1) - bbox(1) + 1;
                        xy_local(:,2) = xy(:,2) - bbox(2) + 1;
                        max_val = double(intmax('uint16'));
                        marker_color_uint16 = uint16(pop_color * max_val);

                        for cell_i = 1:size(xy_local, 1)
                            cx = round(xy_local(cell_i, 1));
                            cy = round(xy_local(cell_i, 2));
                            row_start = max(1, cy - half_k(1));
                            row_end = min(img_size(1), cy + half_k(1));
                            col_start = max(1, cx - half_k(2));
                            col_end = min(img_size(2), cx + half_k(2));
                            k_row_start = row_start - (cy - half_k(1)) + 1;
                            k_row_end = kernel_size(1) - ((cy + half_k(1)) - row_end);
                            k_col_start = col_start - (cx - half_k(2)) + 1;
                            k_col_end = kernel_size(2) - ((cx + half_k(2)) - col_end);
                            kernel_region = marker_kernel(k_row_start:k_row_end, k_col_start:k_col_end);
                            for ch_idx = 1:nch
                                for rgb_ch = 1:3
                                    img_region = single_ch_images{ch_idx}(row_start:row_end, col_start:col_end, rgb_ch);
                                    img_region(kernel_region) = marker_color_uint16(rgb_ch);
                                    single_ch_images{ch_idx}(row_start:row_end, col_start:col_end, rgb_ch) = img_region;
                                end
                            end
                        end
                    end
                end
            end

            for ch = 1:nch
                well_montage{ch}{j} = single_ch_images{ch};
            end
        end
    end

    ngrid = ceil(sqrt(n_panels));

    well_montage_tiled = cell(1, nch);
    for ch = 1:nch
        well_montage_tiled{ch} = imtile(well_montage{ch}, 'GridSize', [1, ngrid], 'BackgroundColor', 'w', 'BorderSize', 10);
    end
    well_montage_merged = imtile(well_montage_merged, 'GridSize', [1, ngrid], 'BackgroundColor', 'w', 'BorderSize', 10);

    title_display = condi;
    fname = strrep(condi, '*', 'star');
    fs = 20;

    % save per-channel images
    for ch = 1:nch
        f = figure;
        imshow(well_montage_tiled{ch});
        title(title_display, 'FontSize', fs);
        f.Position = f.Position + [0 -200 0 250];
        annotation('textbox', textbox, ...
            'String', channels{ch}, ...
            'Color', 'black', ...
            'FontSize', fs, ...
            'FontWeight', 'bold', ...
            'EdgeColor', 'none')
        pause(1);
        exportgraphics(gca, fullfile(montageDir, sprintf('%s_ch%d.png', fname, ch)), 'Resolution', 300);
        close(f);
    end

    % save merged image
    channelLength = cellfun(@(x)length(x), channels);
    f = figure;
    imshow(well_montage_merged);
    % example ENDO percent overlay:
    % if plot_percent
    %     dat_sub_1 = dat(strcmp(dat.channel_x, 'SOX17') & strcmp(dat.channel_y, 'TFAP2C') & dat.condition_idx == i, :);
    %     dat_sub_2 = dat(strcmp(dat.channel_x, 'SOX17') & strcmp(dat.channel_y, 'EOMES') & dat.condition_idx == i, :);
    %     title({title_display, ...
    %         ['\color[rgb]{0.247, 1.0, 1.0}TFAP2C+SOX17+\color{black}=' num2str(round(dat_sub_1(1,:).pp, 1)), '%', ...
    %          ' ', '\color{blue}SOX17+\color{black}=' num2str(round(dat_sub_2(1,:).pp + dat_sub_2(1,:).np, 1)), '%'], ...
    %         ['\color[rgb]{0.9, 0, 1}EOMES+SOX17+\color{black}=' num2str(round(dat_sub_2(1,:).pp, 1)), '%', ...
    %          ' ', '\color{red}EOMES+\color{black}=' num2str(round(dat_sub_2(1,:).pp + dat_sub_2(1,:).pn, 1)), '%'], ...
    %         ['\color{black}Ncell\color{black}=' num2str(round(dat_sub_1(1,:).mean_ncells, 1)) ...
    %          '(' num2str(round(dat_sub_1(1,:).sd_ncells, 1)) ')']}, 'FontSize', fs);
    % end
    title({title_display}, 'FontSize', fs);

    f.Position = f.Position + [0 -200 0 250];
    annotation('textbox', textbox, ...
        'String', channels{1}, ...
        'Color', colors{1}, ...
        'FontSize', fs, ...
        'FontWeight', 'bold', ...
        'EdgeColor', 'none')
    % example ENDO order: label_order = [4 2 3];
    label_order = 2:nch;
    for j = 2:length(channels)
        annotation('textbox', textbox + [sum(channelLength(label_order(1:j-1))) * 0.035 + 0.0025 * (j - 1) 0 0 0], ...
            'String', channels{label_order(j-1)}, ...
            'Color', colors{label_order(j-1)}, ...
            'FontSize', fs, ...
            'FontWeight', 'bold', ...
            'EdgeColor', 'none')
    end
    pause(1);
    exportgraphics(gca, fullfile(montageDir, [fname '_merged.png']), 'Resolution', 300);
    close(f);
end

end
