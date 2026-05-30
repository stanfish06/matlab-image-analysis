function [ target_indices, source_indices, total_cost] = newhungarianlinker(source_info,...
    target_info, imsize, max_distance)
    %The linker constructs a cost function for assigning each point in the
    %source set to a point in the target set, or to nothing. 
    %Source info and target info contain columns for x, y, area, and
    %intensity of each point in the source and target frames:
    source = source_info(:,1:2);
    target = target_info(:,1:2);
    source_a = source_info(:,3);
    target_a = target_info(:,3);
    source_i = source_info(:,4);
    target_i = target_info(:,4);
    
    n_source_points = size(source, 1);
    n_target_points = size(target, 1);

    %find the source and target areas in the lower 20%
    ssarea = sortedQuantileValue(source_a, 0.2, 'floor');
    starea = sortedQuantileValue(target_a, 0.2, 'floor');
    %find the source and target intensities in the upper 80%
    bsintensity = sortedQuantileValue(source_i, 0.8, 'ceil');
    btintensity = sortedQuantileValue(target_i, 0.8, 'ceil');

    % Build frame-to-frame cost matrix. Rows are source points and columns
    % are target points.
    square_dist = (source(:,1) - target(:,1)').^2 + (source(:,2) - target(:,2)').^2;
    diff_i = 2*abs(source_i - target_i')./(source_i + target_i');

    % If the source cell is small and bright, reduce the cost of linking to
    % other small bright cells (can be refined further to account for
    % expected direction of dividing cell movement, added criteria to
    % determine whether a cell is likely to divide (shape)).
    divs = ones(n_source_points, n_target_points);
    source_div_like = source_a < ssarea & source_i > bsintensity;
    target_div_like = target_i' > btintensity & target_a' < starea;
    divs(source_div_like, target_div_like) = 0.5;

    %A1 = square_dist.*(1 + diff_a).*(1 + diff_i);
    A1 = square_dist.*(1 + diff_i).*divs;

    % Deal with maximal linking distance: we simply mark these links as already
    % treated, so that they can never generate a link.
    A1 ( A1 > max_distance * max_distance ) = Inf;
    
    %Determine the cost to be ~>110% of maximum cell-cell linking cost
    linking_costs = A1(~isinf(A1)); %extract the finite costs
    if ~isempty(linking_costs)
        thresh = 1.1*max(linking_costs,[],'all'); %use max value
    else
        thresh = 1; %arbitrary value, linking costs are all Inf
    end

    %Define the weight for disappearance of each point from the source based on
    %distance from edges of frame
    edge_buffer = 50; %set some buffer so costs at the edge are less dramatically low
    wd = 3; %make a weight to adjust the cost for disappearance
    dx = min(source(:,1), imsize(1) - source(:,1)) + edge_buffer;
    dy = min(source(:,2), imsize(2) - source(:,2)) + edge_buffer;
    ds = wd*(dx.*dy./(dx + dy)).^2;
    ds(ds > max_distance^2) = thresh;

    %Define the cost for appearance of each point in the target based on
    %distance from edges of frame
    wb = 3; %make a weight to adjust the cost for appearance
    dx = min(target(:,1), imsize(1) - target(:,1)) + edge_buffer;
    dy = min(target(:,2), imsize(2) - target(:,2)) + edge_buffer;
    bs = wb*(dx.*dy./(dx + dy)).^2;
    bs(bs > max_distance^2) = thresh;

    A2 = Inf(n_source_points);
    A2(eye(n_source_points) == 1) = ds;

    A3 = Inf(n_target_points);
    A3(eye(n_target_points) == 1) = bs;

    A4 = transpose(A1);

    CM = [A1, A2; A3, A4];

    % Find the optimal assignment using Yi Cao's optimization algorithm
    [CM_indices, total_cost] = lapjv(CM); %munkres(CM);

    target_indices = CM_indices(1:n_source_points);
    source_indices = CM_indices((n_source_points+1):(n_source_points+n_target_points)) - n_target_points;

    target_indices(target_indices > n_target_points) = -1;
    source_indices(source_indices < 1) = -1;

end

function value = sortedQuantileValue(values, fraction, roundingMode)
    sortedValues = sort(values);
    if strcmp(roundingMode, 'floor')
        idx = floor(fraction*length(sortedValues));
    else
        idx = ceil(fraction*length(sortedValues));
    end
    idx = max(1, min(length(sortedValues), idx));
    value = sortedValues(idx);
end
