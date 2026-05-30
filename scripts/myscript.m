% MYSCRIPT  Build cropped per-channel MIP montages from stitched images [scratch script].
%   Reads stitched images for positions 0-17, computes per-channel MIPs, crops
%   each to its content bounding box, and tiles them into montages saved under
%   'mip_montage'. Ad-hoc scratch script; edit the file format and position range.
bare_fname = 'stitched_p%.4d_w%.4d_t0000.tif';
img_ch1 = {};
img_ch2 = {};
img_ch3 = {};
img_ch4 = {};
bboxes = {};  % Store bounding boxes for coordinate transformation
parfor i = 0:17
    img_ch1{i+1} = max(tiffreadVolume(sprintf(bare_fname,i,0)),[],3);
    img_ch2{i+1} = max(tiffreadVolume(sprintf(bare_fname,i,1)),[],3);
    img_ch3{i+1} = max(tiffreadVolume(sprintf(bare_fname,i,2)),[],3);
    img_ch4{i+1} = max(tiffreadVolume(sprintf(bare_fname,i,3)),[],3);
    img_ch1_mask = img_ch1{i+1} > 0;
    rp = regionprops(img_ch1_mask,'BoundingBox');
    sz = size(img_ch1{i+1});
    bbox = floor(rp(1).BoundingBox);
    bbox(bbox == 0) = 1;
    bboxes{i+1} = bbox;  % Store bbox for later coordinate transformation
    img_ch1{i+1} = img_ch1{i+1}(bbox(2):min(bbox(2)+bbox(4),sz(1)),bbox(1):min(bbox(1)+bbox(3),sz(2)));
    img_ch2{i+1} = img_ch2{i+1}(bbox(2):min(bbox(2)+bbox(4),sz(1)),bbox(1):min(bbox(1)+bbox(3),sz(2)));
    img_ch3{i+1} = img_ch3{i+1}(bbox(2):min(bbox(2)+bbox(4),sz(1)),bbox(1):min(bbox(1)+bbox(3),sz(2)));
    img_ch4{i+1} = img_ch4{i+1}(bbox(2):min(bbox(2)+bbox(4),sz(1)),bbox(1):min(bbox(1)+bbox(3),sz(2)));
end
%%
montageDir = 'mip_montage';
if ~exist(montageDir, 'dir')
    mkdir(montageDir);
end
%%
load('meta.mat');
channels = meta.channelLabel;
colors = {[0.5 0.5 0.5], [1 0 0], [0 1 0], [0 0 1]};
% colors = {[0.5 0.5 0.5], [1 0 0], [0 1 1], [1 0 0]};
scales = [0.5, 1, 1.0, 1.0];
textbox = [0.1 0.18 0.8 0.1];
plot_percent = true;
if plot_percent
    dat = readtable("scatterPercent.csv");
end
%%
ch_bounds = [
    0.005, 0.005;
    0.01, 0.05;
    0.01, 0.04;
    0.01, 0.07;
    ];
%%
%% Population overlay settings
do_population_overlay = true;

% Load position data - choose the appropriate file for your data
% Use 'positions.mat' for full-image positions
% Use 'positions_colony.mat' for colony-extracted positions [BUG]
position_file = 'positions.mat';  % Change this as needed
if do_population_overlay
    load(position_file, 'positions');
end

% Channel thresholds for positive/negative classification
% Format: [DAPI, EOMES, TFAP2C, SOX17] - adjust based on your data
% These should match the thresholds used in makeScatter.m
channelThresholds = [0, 250, 350, 450];

% Define populations to overlay
% Each row: {name, [channel indices], [+1 for positive, -1 for negative], color, marker_shape, marker_size_pixels}
% Channel indices: 1=DAPI, 2=EOMES, 3=TFAP2C, 4=SOX17
% Signs: +1 means above threshold (positive), -1 means below threshold (negative)
% Marker shapes: 'disk', 'square', 'diamond', 'octagon' (uses strel)
% Marker size: radius in pixels
populations = {
    % Double positive/negative for EOMES vs SOX17
    %'EOMES+SOX17+', [2, 4], [1, 1], [0.9 0 1], 'disk', 5; % Magenta disks
    %'EOMES+SOX17-', [2, 4], [1, -1], [1 0 0], 'diamond', 5; % Red diamonds
    %'EOMES-SOX17+', [2, 4], [-1, 1], [0 0 1], 'square', 4; % Blue squares
    %'EOMES-SOX17-', [2, 4], [-1, -1], [0.5 0.5 0.5], 'disk', 2; % Gray dots (often many cells)

    % Double positive/negative for TFAP2C vs SOX17
    'TFAP2C+SOX17+', [3, 4], [1, 1], [1 0 0], 'disk', 5; % Cyan disks
    %'TFAP2C+SOX17-', [3, 4], [1, -1], [0 1 0], 'square', 4; % Green squares
    %'TFAP2C-SOX17+', [3, 4], [-1, 1], [0 0 1], 'diamond', 5; % Blue diamonds

    % Triple positive/negative (uncomment as needed)
    %'TriplePos', [2,3,4], [1,1,1], [1 1 0], 'octagon', 6; % Yellow octagons
    %'TripleNeg', [2,3,4], [-1,-1,-1], [0 0 0], 'disk', 2; % Black dots
};

% Helper function to create marker kernel using strel
% Returns a logical matrix for the marker shape
createMarkerKernel = @(shape, radius) getnhood(strel(shape, radius));
%%
%1:length(meta.conditions)
for i = 1:length(meta.conditions)
    condi = meta.conditions{i};
    start_i = meta.conditionStartPos(i);
    if i < length(meta.conditions)
        end_i = meta.conditionStartPos(i + 1) - 1;
    else
        end_i = length(img_ch1);
    end
    well_montage_ch1 = {};
    well_montage_ch2 = {};
    well_montage_ch3 = {};
    well_montage_ch4 = {};
    well_montage_merged = {};
    for j = 1:(end_i - start_i + 1)
        pos_idx = start_i + j - 1;
        img = cat(3, img_ch1{pos_idx}, img_ch2{pos_idx}, img_ch3{pos_idx}, img_ch4{pos_idx});
        % img = img_all{pos_idx};
        img = max(img, [], 4);
        for k = 1:4
            % img(:,:,k) = imtophat(img(:,:,k), strel('disk', 200));
            if k == 1
                img(:,:,k) = medfilt2(imadjust(img(:,:,k),stretchlim(img(:,:,k))'+ch_bounds(k,:),[0 1]), [3 3]);
            else
                img(:,:,k) = medfilt2(imadjust(img(:,:,k),ch_bounds(k,:),[0 1]), [3 3]);
            end
        end
        well_montage_ch1{j} = img(:,:,1);
        well_montage_ch2{j} = img(:,:,2);
        well_montage_ch3{j} = img(:,:,3);
        well_montage_ch4{j} = img(:,:,4);
    well_montage_merged{j} = uint16(zeros([size(img, [1 2]) 3]));
    for ci = 1:3
        for k = 1:4
        well_montage_merged{j}(:,:,ci) = well_montage_merged{j}(:,:,ci) + img(:,:,k) * scales(k) * colors{k}(ci);
        end
    end

        % Inject population markers directly into the merged image pixels
        if do_population_overlay && ~isempty(populations)
            % Get bounding box for this position
            bbox = bboxes{pos_idx};
            img_size = size(well_montage_merged{j});

            % Get position data (handle both cell array and array of Position objects)
            if iscell(positions)
                pos_data = positions{pos_idx};
            else
                pos_data = positions(pos_idx);
            end

            % Skip if no cells in this position
            if ~isempty(pos_data.cellData) && isfield(pos_data.cellData, 'XY') && ~isempty(pos_data.cellData.XY)
                % Overlay each population
                for pop_idx = 1:size(populations, 1)
                    pop_channels = populations{pop_idx, 2};
                    pop_signs = populations{pop_idx, 3};
                    pop_color = populations{pop_idx, 4};
                    pop_shape = populations{pop_idx, 5};
                    pop_radius = populations{pop_idx, 6};

                    % Create marker kernel using strel
                    marker_kernel = createMarkerKernel(pop_shape, pop_radius);
                    kernel_size = size(marker_kernel);
                    half_k = floor(kernel_size / 2);

                    % Get cells matching this population
                    nucLevel = pos_data.cellData.nucLevel;
                    bg = pos_data.cellData.background;

                    % Start with all cells
                    cell_idx = true(size(nucLevel, 1), 1);

                    % Apply each channel criterion
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
                        % Transform coordinates: subtract bbox offset
                        xy_local = zeros(size(xy));
                        xy_local(:,1) = xy(:,1) - bbox(1) + 1;  % X (column)
                        xy_local(:,2) = xy(:,2) - bbox(2) + 1;  % Y (row)

                        % Convert color to uint16 (scale to max intensity)
                        max_val = double(intmax('uint16'));
                        marker_color_uint16 = uint16(pop_color * max_val);

                        % Stamp marker at each cell location
                        for cell_i = 1:size(xy_local, 1)
                            cx = round(xy_local(cell_i, 1));  % column (x)
                            cy = round(xy_local(cell_i, 2));  % row (y)

                            % Calculate stamp bounds with boundary checking
                            row_start = max(1, cy - half_k(1));
                            row_end = min(img_size(1), cy + half_k(1));
                            col_start = max(1, cx - half_k(2));
                            col_end = min(img_size(2), cx + half_k(2));

                            % Calculate corresponding kernel indices
                            k_row_start = row_start - (cy - half_k(1)) + 1;
                            k_row_end = kernel_size(1) - ((cy + half_k(1)) - row_end);
                            k_col_start = col_start - (cx - half_k(2)) + 1;
                            k_col_end = kernel_size(2) - ((cx + half_k(2)) - col_end);

                            % Get the kernel region
                            kernel_region = marker_kernel(k_row_start:k_row_end, k_col_start:k_col_end);

                            % Stamp the marker color where kernel is true
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
    end
    n_panels = length(well_montage_ch1);
    ngrid = ceil(sqrt(n_panels));
    disp(ngrid);

    well_montage_ch1 = imtile(well_montage_ch1, 'GridSize', [1, ngrid], 'BackgroundColor', 'w', 'BorderSize', 10);
    well_montage_ch2 = imtile(well_montage_ch2, 'GridSize', [1, ngrid], 'BackgroundColor', 'w', 'BorderSize', 10);
    well_montage_ch3 = imtile(well_montage_ch3, 'GridSize', [1, ngrid], 'BackgroundColor', 'w', 'BorderSize', 10);
    well_montage_ch4 = imtile(well_montage_ch4, 'GridSize', [1, ngrid], 'BackgroundColor', 'w', 'BorderSize', 10);
    well_montage_merged = imtile(well_montage_merged, 'GridSize', [1, ngrid], 'BackgroundColor', 'w', 'BorderSize', 10);

    % title_display = strrep(condi," ",newline);
    title_display = condi;
    fname = strrep(condi, '*', 'star');

    fs = 20;
    f = figure;
    imshow(well_montage_ch1);
    title(title_display, 'FontSize', fs);
    f.Position = f.Position +  [0 -200 0 250];
    annotation('textbox', textbox, ...
        'String', channels{1}, ...
        'Color', 'black', ...
        'FontSize', fs, ...
        'FontWeight', 'bold', ...
        'EdgeColor', 'none')
    pause(1);
    exportgraphics(gca, fullfile(montageDir, [fname '_ch1.png']), 'Resolution', 300);
    close(f);

    f = figure;
    imshow(well_montage_ch2);
    title(title_display, 'FontSize', fs);
    f.Position = f.Position +  [0 -200 0 250];
    annotation('textbox', textbox, ...
        'String', channels{2}, ...
        'Color', 'black', ...
        'FontSize', fs, ...
        'FontWeight', 'bold', ...
        'EdgeColor', 'none')
    pause(1);
    exportgraphics(gca, fullfile(montageDir, [fname '_ch2.png']), 'Resolution', 300);
    close(f);

    f = figure;
    imshow(well_montage_ch3);
    title(title_display, 'FontSize', fs);
    f.Position = f.Position +  [0 -200 0 250];
    annotation('textbox', textbox, ...
        'String', channels{3}, ...
        'Color', 'black', ...
        'FontSize', fs, ...
        'FontWeight', 'bold', ...
        'EdgeColor', 'none')
    pause(1);
    exportgraphics(gca, fullfile(montageDir, [fname '_ch3.png']), 'Resolution', 300);
    close(f);

    f = figure;
    imshow(well_montage_ch4);
    title(title_display, 'FontSize', fs);
    f.Position = f.Position +  [0 -200 0 250];
    annotation('textbox', textbox, ...
        'String', channels{4}, ...
        'Color', 'black', ...
        'FontSize', fs, ...
        'FontWeight', 'bold', ...
        'EdgeColor', 'none')
    pause(1);
    exportgraphics(gca, fullfile(montageDir, [fname '_ch4.png']), 'Resolution', 300);
    close(f);

    channelLength = cellfun(@(x)length(x), channels);
    f = figure;
    imshow(well_montage_merged);
    if plot_percent
        dat_sub_1 = dat(strcmp(dat.channel_x, 'SOX17') & strcmp(dat.channel_y, 'TFAP2C') & dat.condition_idx == i, :);
        dat_sub_2 = dat(strcmp(dat.channel_x, 'SOX17') & strcmp(dat.channel_y, 'EOMES') & dat.condition_idx == i, :);
        title({title_display, ['\color[rgb]{0.247, 1.0, 1.0}TFAP2C+SOX17+\color{black}=' num2str(round(dat_sub_1(1,:).pp, 1)), '%', ' ', '\color{blue}SOX17+\color{black}=' num2str(round(dat_sub_2(1,:).pp + dat_sub_2(1,:).np, 1)), '%'], ['\color[rgb]{0.9, 0, 1}EOMES+SOX17+\color{black}=' num2str(round(dat_sub_2(1,:).pp, 1)), '%', ' ', '\color{red}EOMES+\color{black}=' num2str(round(dat_sub_2(1,:).pp + dat_sub_2(1,:).pn, 1)), '%'], ['\color{black}Ncell\color{black}=' num2str(round(dat_sub_1(1,:).mean_ncells, 1)) '(' num2str(round(dat_sub_1(1,:).sd_ncells, 1)) ')']}, 'FontSize', fs);
    else
        title({title_display}, 'FontSize', fs);
    end

    f.Position = f.Position +  [0 -200 0 250];
    annotation('textbox', textbox, ...
        'String', channels{1}, ...
        'Color', colors{1}, ...
        'FontSize', fs, ...
        'FontWeight', 'bold', ...
        'EdgeColor', 'none')
    label_order = [4 2 3];
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
