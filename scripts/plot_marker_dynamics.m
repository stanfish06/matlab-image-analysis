%% Plot Marker Dynamics
% Simple script to plot nuclear/cytoplasmic marker intensity over time
% Requires: positions array and meta from analyze_live_stan.m

%% User parameters
timeInterval = '15 min';  % manually set time interval
globalTimeUnit = 'hour';   % display/fit unit: 'min' or 'hour'
treatmentTime = 1;        % frame number of treatment (1 = no pre-treatment)
minNCells = 10;           % minimum cells per timepoint
SNRcutoff = 1;            % signal-to-noise cutoff (1 = no filtering)
outputDir = '.';          % directory to save plots

% Mode per channel:
%   'N' = nuclear (raw), 'C' = cytoplasmic (raw)
%   'N-bg' = nuclear - background, 'C-bg' = cytoplasmic - background
%   'N:C' = ratio, 'C:N' = inverse ratio
% Must match order of meta.channelLabel
channelModes = {'N-bg', 'N-bg'};  % e.g., {'N:C', 'N-bg'} for SMAD + H2B

% Normalization: 'none', 'baseline' (divide by baseline),
% 'baseline_subtract' (subtract baseline)
normalization = 'baseline';
baselineFrames = [1 5];  % range of frames for baseline

% Drift correction settings.
% Drift is fit from one reference condition/channel/mode and then applied
% to all conditions/channels by division.
enableDriftCorrection = true;
driftFitCondition = 'mTe';   % condition name or index
driftFitChannel = 'SOX2';    % channel label or index
driftFitMode = 'N-bg';       % mode for reference fit
driftFitModel = 'linear';    % currently supports: 'linear'
driftSlopeScale = 1;        % multiply fitted slope (set 1 for default fit)
driftMinPoints = 5;          % minimum valid points required for fit
driftEpsilon = 1e-6;         % denominator safety lower bound

% Error bars: 'none', 'std', 'sem' (standard error of mean)
showError = 'sem';

% Axis limits: set to [] for auto, or [min max] to customize.
% Units follow globalTimeUnit.
xLimits = [0 5];  % e.g., [0 480] if globalTimeUnit='min', [0 8] if 'hour'
yLimits = [0.5, 1.3];  % or cell array per channel: {[0 1000], [0 500]}

% Legend placement:
%   'inside'        -> use legendInsideLocation
%   'rightoutside'  -> legend outside right side
%   'bottomoutside' -> legend outside bottom
legendPlacement = 'bottomoutside';
legendInsideLocation = 'best';
legendNumColumns = 1;  % [] for auto, or positive integer

% Channel-plot layout and labels
% Figure size is [width height] in pixels.
% Plot aspect is x:y ratio in the plotting box (e.g., 2 -> x is 2x y).
channelFigureSize = [900 700];
channelPlotAspect = 1;
channelTitleOverrides = {'SOX2', 'H2B'};  % {} auto, scalar text for all, or cell per channel
channelXLabel = 'time (hour) after 15min ActD priming';          % '' uses auto label
channelYLabelOverrides = {'GFP::SOX2', 'iRFP::H2B'}; % {} auto, scalar text for all, or cell per channel

% Optional slope reporting for channel dynamics.
% Fits per-position linear models after normalization/drift-correction,
% then reports mean slope +/- SE across positions in each condition.
showConditionSlopes = true;  % optional
slopeFitUseXLimits = true;    % fit only within xLimits when true
slopeMinPoints = 5;           % minimum timepoints per position trace for slope fit

% Cell-number dynamics plot settings
plotCellNumberDynamics = true;
cellNumberNormalization = 'none';  % 'none', 'baseline', or 'baseline_subtract'
cellNumberYLimits = [];            % [] for auto, or [min max]
cellNumberFigureSize = [900 700];
cellNumberPlotAspect = 1;
cellNumberTitle = '';   % '' uses auto title
cellNumberXLabel = '';  % '' uses auto label
cellNumberYLabel = '';  % '' uses auto label

%% Basic setup and validation
meta.timeInterval = timeInterval;

if numel(channelModes) ~= numel(meta.channelLabel)
    error('channelModes must have one entry per channel in meta.channelLabel.');
end

for mi = 1:numel(channelModes)
    validateMode(channelModes{mi});
end
validateMode(driftFitMode);

posPerConditionArray = meta.posPerCondition;
if isscalar(posPerConditionArray)
    posPerConditionArray = repmat(posPerConditionArray, [1, meta.nWells]);
end

tmax = meta.nTime;
timeUnit = strsplit(timeInterval, ' ');
if numel(timeUnit) < 2
    error('timeInterval must be in format ''<value> <unit>'', e.g. ''15 min''.');
end
dt = str2double(timeUnit{1});
if ~isfinite(dt) || dt <= 0
    error('timeInterval value must be a positive number.');
end
dtMinutes = dt * unitToMinutesFactor(timeUnit{2});
globalTimeUnit = canonicalTimeUnit(globalTimeUnit);
if strcmp(globalTimeUnit, 'min')
    timeScaleFromMinutes = 1;
elseif strcmp(globalTimeUnit, 'hour')
    timeScaleFromMinutes = 1/60;
else
    error('globalTimeUnit must be ''min'' or ''hour''.');
end
dt = dtMinutes * timeScaleFromMinutes;
t = ((1:tmax) - treatmentTime) * dt;
xData = t;
timeAxisLabel = globalTimeUnit;
slopeUnitLabel = ['/' globalTimeUnit];

if ~isempty(xLimits) && numel(xLimits) ~= 2
    error('xLimits must be [] or [xmin xmax].');
end

if ~isempty(legendNumColumns)
    if ~isscalar(legendNumColumns) || ~isfinite(legendNumColumns) || ...
            legendNumColumns < 1 || legendNumColumns ~= round(legendNumColumns)
        error('legendNumColumns must be [] or a positive integer.');
    end
end

validateFigureSizeOption(channelFigureSize, 'channelFigureSize');
validateFigureSizeOption(cellNumberFigureSize, 'cellNumberFigureSize');
validateAspectRatioOption(channelPlotAspect, 'channelPlotAspect');
validateAspectRatioOption(cellNumberPlotAspect, 'cellNumberPlotAspect');

validateScalarTextOption(channelXLabel, 'channelXLabel');
validatePerChannelTextOption(channelTitleOverrides, numel(meta.channelLabel), ...
    'channelTitleOverrides');
validatePerChannelTextOption(channelYLabelOverrides, numel(meta.channelLabel), ...
    'channelYLabelOverrides');
validateScalarTextOption(cellNumberTitle, 'cellNumberTitle');
validateScalarTextOption(cellNumberXLabel, 'cellNumberXLabel');
validateScalarTextOption(cellNumberYLabel, 'cellNumberYLabel');

if slopeFitUseXLimits && ~isempty(xLimits)
    slopeFitMask = t >= xLimits(1) & t <= xLimits(2);
else
    slopeFitMask = true(size(t));
end

baselineIdx = baselineFrames(1):baselineFrames(end);
baselineIdx = baselineIdx(baselineIdx >= 1 & baselineIdx <= tmax);
if isempty(baselineIdx)
    error('baselineFrames does not overlap with valid frame range 1:%d.', tmax);
end

%% Build drift correction curve (optional)
hasValidDriftFit = false;
driftCorrectionCurve = ones(1, tmax);
driftInfo = '';

if enableDriftCorrection
    if ~strcmp(normalization, 'baseline')
        warning(['Drift correction is supported only when normalization is ' ...
            '''baseline''. Current setting: ''%s''. Skipping drift correction.'], ...
            normalization);
    else
        fitConditionIdx = resolveLabelOrIndex(driftFitCondition, ...
            meta.conditions, 'condition');
        fitChannelIdx = resolveLabelOrIndex(driftFitChannel, ...
            meta.channelLabel, 'channel');

        fitPositions = getConditionPositionIndices(meta, ...
            posPerConditionArray, fitConditionIdx);
        fitPosSubset = positions(fitPositions);

        fitTraces = collectConditionTraces(fitPosSubset, numel(fitPositions), ...
            fitChannelIdx, driftFitMode, SNRcutoff, minNCells, tmax);
        fitTraces = normalizeTraces(fitTraces, normalization, baselineIdx);

        fitMean = mean(fitTraces, 1, 'omitnan');
        validFit = isfinite(fitMean);
        fitWindowMask = true(size(t));
        if ~isempty(xLimits)
            if numel(xLimits) ~= 2
                error('xLimits must be [] or [xmin xmax].');
            end
            fitWindowMask = t >= xLimits(1) & t <= xLimits(2);
        end
        validFit = validFit & fitWindowMask;
        if sum(validFit) < driftMinPoints
            warning(['Only %d valid points available for drift fit, but ' ...
                'driftMinPoints=%d. Skipping drift correction.'], ...
                sum(validFit), driftMinPoints);
        else
            switch lower(driftFitModel)
                case 'linear'
                    coeffs = polyfit(t(validFit), fitMean(validFit), 1);
                    fitSlope = coeffs(1);
                    fitIntercept = coeffs(2);
                    tCenter = mean(t(validFit), 'omitnan');

                    % Keep the same value at fit-window center while scaling slope.
                    fitIntercept = fitIntercept + ...
                        (1 - driftSlopeScale) * fitSlope * tCenter;
                    fitSlope = fitSlope * driftSlopeScale;

                    driftCorrectionCurve = fitSlope * t + fitIntercept;
                otherwise
                    error('Unsupported driftFitModel: %s', driftFitModel);
            end

            corrBaseline = mean(driftCorrectionCurve(baselineIdx), 'omitnan');
            if ~isfinite(corrBaseline) || abs(corrBaseline) < driftEpsilon
                warning(['Drift fit baseline is invalid (%.3g). ' ...
                    'Skipping drift correction.'], corrBaseline);
            else
                driftCorrectionCurve = driftCorrectionCurve / corrBaseline;
                badMask = ~isfinite(driftCorrectionCurve) | ...
                    driftCorrectionCurve <= driftEpsilon;
                if any(badMask)
                    warning(['Drift correction curve has %d non-finite or ' ...
                        'small values; clamping to driftEpsilon=%.3g.'], ...
                        sum(badMask), driftEpsilon);
                    driftCorrectionCurve(badMask) = driftEpsilon;
                end
                hasValidDriftFit = true;
                driftInfo = sprintf('%s / %s / %s (%s fit)', ...
                    meta.conditions{fitConditionIdx}, ...
                    meta.channelLabel{fitChannelIdx}, ...
                    driftFitMode, lower(driftFitModel));
                disp(['Using drift correction from: ' driftInfo]);
                disp(sprintf(['Drift fit slope=%.3g%s (slope scale %.3g), ' ...
                    'curve range in fit window [%.4f, %.4f]'], ...
                    fitSlope, slopeUnitLabel, driftSlopeScale, ...
                    min(driftCorrectionCurve(validFit)), ...
                    max(driftCorrectionCurve(validFit))));
            end
        end
    end
end

%% Plot each channel
for ci = 1:numel(meta.channelLabel)
    channelName = meta.channelLabel{ci};
    mode = channelModes{ci};
    colors = lines(meta.nWells);

    fig = figure('Visible', 'off');
    configureFigureSize(fig, channelFigureSize);
    hold on;
    legendEntries = {};

    % Loop through each condition separately
    for wi = 1:meta.nWells
        conditionPositions = getConditionPositionIndices(meta, ...
            posPerConditionArray, wi);
        nPos = numel(conditionPositions);
        posSubset = positions(conditionPositions);

        allTraces = collectConditionTraces(posSubset, nPos, ci, mode, ...
            SNRcutoff, minNCells, tmax);
        allTraces = normalizeTraces(allTraces, normalization, baselineIdx);

        if hasValidDriftFit
            allTraces = bsxfun(@rdivide, allTraces, driftCorrectionCurve);
        end

        condSlope = NaN;
        condSlopeSE = NaN;
        nSlopeFits = 0;
        if showConditionSlopes
            [condSlope, condSlopeSE, nSlopeFits] = estimateConditionSlope( ...
                allTraces, t, slopeFitMask, slopeMinPoints);
        end

        % Compute mean and error
        yMean = mean(allTraces, 1, 'omitnan');
        if strcmp(showError, 'std')
            yErr = std(allTraces, 0, 1, 'omitnan');
        elseif strcmp(showError, 'sem')
            nValid = sum(~isnan(allTraces), 1);
            yErr = std(allTraces, 0, 1, 'omitnan') ./ sqrt(max(nValid, 1));
            yErr(nValid == 0) = NaN;
        else
            yErr = [];
        end

        % Plot mean line with error bars
        if ~isempty(yErr)
            errorbar(xData, yMean, yErr, 'Color', colors(wi,:), ...
                'LineWidth', 2, 'CapSize', 4);
        else
            plot(xData, yMean, 'Color', colors(wi,:), 'LineWidth', 2);
        end

        if showConditionSlopes
            if isnan(condSlope)
                legendEntries{end+1} = sprintf('%s (slope n/a)', ...
                    meta.conditions{wi}); %#ok<AGROW>
            elseif isnan(condSlopeSE)
                legendEntries{end+1} = sprintf('%s (s=%.3g%s, n=%d)', ...
                    meta.conditions{wi}, condSlope, slopeUnitLabel, ...
                    nSlopeFits); %#ok<AGROW>
            else
                legendEntries{end+1} = sprintf('%s (s=%.3g +/- %.3g%s)', ...
                    meta.conditions{wi}, condSlope, condSlopeSE, ...
                    slopeUnitLabel); %#ok<AGROW>
            end
        else
            legendEntries{end+1} = meta.conditions{wi}; %#ok<AGROW>
        end
    end

    % Format plot
    defaultXLabel = ['time (' timeAxisLabel ')'];
    xLabelStr = resolveScalarText(channelXLabel, defaultXLabel, ...
        'channelXLabel');
    xlabel(xLabelStr, 'FontSize', 24, 'FontWeight', 'Bold');

    defaultYLabel = channelName;
    yLabelStr = resolvePerChannelText(channelYLabelOverrides, ci, ...
        defaultYLabel, 'channelYLabelOverrides', numel(meta.channelLabel));
    ylabel(yLabelStr, 'FontSize', 24, 'FontWeight', 'Bold');
    [legendLocation, legendOrientation] = resolveLegendPlacement( ...
        legendPlacement, legendInsideLocation);
    lgd = legend(legendEntries, 'Location', legendLocation);
    applyLegendStyle(lgd, legendOrientation, legendNumColumns);

    titleDefault = channelName;
    if hasValidDriftFit
        titleDefault = [titleDefault ' (drift corrected)']; %#ok<AGROW>
    end
    titleStr = resolvePerChannelText(channelTitleOverrides, ci, ...
        titleDefault, 'channelTitleOverrides', numel(meta.channelLabel));
    title(titleStr);
    set(gca, 'FontSize', 24, 'FontWeight', 'bold', 'LineWidth', 2);

    % Apply axis limits
    if ~isempty(xLimits)
        xlim(xLimits);
    end
    if ~isempty(yLimits)
        if iscell(yLimits)
            ylim(yLimits{ci});
        else
            ylim(yLimits);
        end
    end
    applyPlotAspect(gca, channelPlotAspect);
    hold off;

    saveas(fig, fullfile(outputDir, [channelName '_dynamics.png']));
    close(fig);
    disp(['Saved: ' channelName '_dynamics.png']);
end

if hasValidDriftFit
    disp(['Applied global drift correction: ' driftInfo]);
end

%% Plot cell-number dynamics
if plotCellNumberDynamics
    colors = lines(meta.nWells);
    fig = figure('Visible', 'off');
    configureFigureSize(fig, cellNumberFigureSize);
    hold on;
    legendEntries = {};

    for wi = 1:meta.nWells
        conditionPositions = getConditionPositionIndices(meta, ...
            posPerConditionArray, wi);
        nPos = numel(conditionPositions);
        allCellCounts = NaN(nPos, tmax);

        for pi = 1:nPos
            pos = positions(conditionPositions(pi));
            nvec = pos.ncells(:)';
            nt = min(numel(nvec), tmax);
            allCellCounts(pi, 1:nt) = nvec(1:nt);
        end

        allCellCounts = normalizeTraces(allCellCounts, ...
            cellNumberNormalization, baselineIdx);

        yMean = mean(allCellCounts, 1, 'omitnan');
        if strcmp(showError, 'std')
            yErr = std(allCellCounts, 0, 1, 'omitnan');
        elseif strcmp(showError, 'sem')
            nValid = sum(~isnan(allCellCounts), 1);
            yErr = std(allCellCounts, 0, 1, 'omitnan') ./ sqrt(max(nValid, 1));
            yErr(nValid == 0) = NaN;
        else
            yErr = [];
        end

        if ~isempty(yErr)
            errorbar(xData, yMean, yErr, 'Color', colors(wi,:), ...
                'LineWidth', 2, 'CapSize', 4);
        else
            plot(xData, yMean, 'Color', colors(wi,:), 'LineWidth', 2);
        end
        legendEntries{end+1} = meta.conditions{wi}; %#ok<AGROW>
    end

    defaultCellXLabel = ['time (' timeAxisLabel ')'];
    cellXLabelStr = resolveScalarText(cellNumberXLabel, defaultCellXLabel, ...
        'cellNumberXLabel');
    xlabel(cellXLabelStr, 'FontSize', 24, 'FontWeight', 'Bold');

    if strcmp(cellNumberNormalization, 'baseline')
        defaultCellYLabel = 'Cell number (normalized)';
    elseif strcmp(cellNumberNormalization, 'baseline_subtract')
        defaultCellYLabel = 'Cell number (delta)';
    else
        defaultCellYLabel = 'Cell number';
    end
    cellYLabelStr = resolveScalarText(cellNumberYLabel, defaultCellYLabel, ...
        'cellNumberYLabel');
    ylabel(cellYLabelStr, 'FontSize', 24, 'FontWeight', 'Bold');
    [legendLocation, legendOrientation] = resolveLegendPlacement( ...
        legendPlacement, legendInsideLocation);
    lgd = legend(legendEntries, 'Location', legendLocation);
    applyLegendStyle(lgd, legendOrientation, legendNumColumns);
    defaultCellTitle = 'Cell number';
    cellTitleStr = resolveScalarText(cellNumberTitle, defaultCellTitle, ...
        'cellNumberTitle');
    title(cellTitleStr);
    set(gca, 'FontSize', 24, 'FontWeight', 'bold', 'LineWidth', 2);

    if ~isempty(xLimits)
        xlim(xLimits);
    end
    if ~isempty(cellNumberYLimits)
        ylim(cellNumberYLimits);
    end
    applyPlotAspect(gca, cellNumberPlotAspect);
    hold off;

    saveas(fig, fullfile(outputDir, 'nCells_dynamics.png'));
    close(fig);
    disp('Saved: nCells_dynamics.png');
end

disp(['All plots saved to: ' outputDir]);

function conditionPositions = getConditionPositionIndices(meta, ...
    posPerConditionArray, conditionIdx)
startPos = meta.conditionStartPos(conditionIdx);
nPos = posPerConditionArray(conditionIdx);
conditionPositions = startPos:(startPos + nPos - 1);
end

function allTraces = collectConditionTraces(posSubset, nPos, channelIdx, ...
    mode, SNRcutoff, minNCells, tmax)
allTraces = NaN(nPos, tmax);
for pi = 1:nPos
    pos = posSubset(pi);

    % Ensure timeTraces is available with the current SNR threshold.
    SNRcutoffVec = ones([numel(pos.dataChannels) 1]) * SNRcutoff;
    pos.makeAvgTimeTraces(SNRcutoffVec);

    traceY = getTraceByMode(pos, channelIdx, mode);
    traceY(pos.ncells < minNCells) = NaN;
    allTraces(pi, :) = traceY;
end
end

function traceY = getTraceByMode(pos, channelIdx, mode)
switch mode
    case 'N'
        traceY = pos.timeTraces.nucLevelAvg(:, channelIdx)';
    case 'C'
        traceY = pos.timeTraces.cytLevelAvg(:, channelIdx)';
    case 'N-bg'
        traceY = pos.timeTraces.nucLevelAvg(:, channelIdx)' - ...
            pos.timeTraces.background(:, channelIdx)';
    case 'C-bg'
        traceY = pos.timeTraces.cytLevelAvg(:, channelIdx)' - ...
            pos.timeTraces.background(:, channelIdx)';
    case 'N:C'
        traceY = pos.timeTraces.ratioAvg(:, channelIdx)';
    case 'C:N'
        traceY = 1 ./ pos.timeTraces.ratioAvg(:, channelIdx)';
    otherwise
        error('Unknown mode: %s', mode);
end
end

function traces = normalizeTraces(traces, normalization, baselineIdx)
if strcmp(normalization, 'baseline') || strcmp(normalization, ...
        'baseline_subtract')
    for ti = 1:size(traces, 1)
        baselineVal = mean(traces(ti, baselineIdx), 'omitnan');
        if ~isfinite(baselineVal)
            traces(ti, :) = NaN;
            continue;
        end

        if strcmp(normalization, 'baseline')
            if baselineVal == 0
                traces(ti, :) = NaN;
            else
                traces(ti, :) = traces(ti, :) / baselineVal;
            end
        else
            traces(ti, :) = traces(ti, :) - baselineVal;
        end
    end
elseif ~strcmp(normalization, 'none')
    error('Unknown normalization mode: %s', normalization);
end
end

function idx = resolveLabelOrIndex(value, labels, what)
if isnumeric(value)
    idx = value;
elseif ischar(value) || isstring(value)
    idx = find(strcmpi(char(value), labels), 1, 'first');
    if isempty(idx)
        error('Unknown %s name: %s', what, char(value));
    end
else
    error('%s must be a string label or numeric index.', what);
end

if idx < 1 || idx > numel(labels)
    error('%s index out of bounds: %d', what, idx);
end
end

function validateMode(mode)
validModes = {'N', 'C', 'N-bg', 'C-bg', 'N:C', 'C:N'};
if ~ismember(mode, validModes)
    error('Unknown mode: %s. Valid modes: %s', mode, strjoin(validModes, ', '));
end
end

function [slopeMean, slopeSE, nFits] = estimateConditionSlope(traces, t, ...
    fitMask, minPoints)
slopes = NaN(size(traces, 1), 1);
for ri = 1:size(traces, 1)
    y = traces(ri, :);
    valid = isfinite(y) & fitMask;
    if sum(valid) >= minPoints
        p = polyfit(t(valid), y(valid), 1);
        slopes(ri) = p(1);
    end
end

validSlopes = isfinite(slopes);
nFits = sum(validSlopes);
if nFits == 0
    slopeMean = NaN;
    slopeSE = NaN;
    return;
end

slopeMean = mean(slopes(validSlopes));
if nFits > 1
    slopeSE = std(slopes(validSlopes), 0) / sqrt(nFits);
else
    slopeSE = NaN;
end
end

function unit = canonicalTimeUnit(unitRaw)
u = lower(strtrim(char(unitRaw)));
if ismember(u, {'min', 'mins', 'minute', 'minutes'})
    unit = 'min';
elseif ismember(u, {'h', 'hr', 'hrs', 'hour', 'hours'})
    unit = 'hour';
else
    error('Unknown time unit: %s. Use min or hour.', char(unitRaw));
end
end

function factor = unitToMinutesFactor(unitRaw)
unit = canonicalTimeUnit(unitRaw);
if strcmp(unit, 'min')
    factor = 1;
else
    factor = 60;
end
end

function [loc, orient] = resolveLegendPlacement(placementRaw, insideLocRaw)
p = lower(strtrim(char(placementRaw)));
insideLoc = lower(strtrim(char(insideLocRaw)));
switch p
    case {'inside'}
        loc = insideLoc;
        orient = 'vertical';
    case {'right', 'rightoutside', 'eastoutside'}
        loc = 'eastoutside';
        orient = 'vertical';
    case {'bottom', 'bottomoutside', 'southoutside'}
        loc = 'southoutside';
        orient = 'horizontal';
    otherwise
        error(['Unknown legendPlacement: %s. Use inside, rightoutside, ' ...
            'or bottomoutside.'], char(placementRaw));
end
end

function applyLegendStyle(lgd, orientation, numColumns)
if strcmp(orientation, 'horizontal')
    set(lgd, 'Orientation', 'horizontal');
end

if ~isempty(numColumns)
    set(lgd, 'NumColumns', numColumns);
end
end

function configureFigureSize(fig, figSizePixels)
if isempty(figSizePixels)
    return;
end

set(fig, 'Units', 'pixels');
pos = get(fig, 'Position');
pos(3:4) = figSizePixels;
set(fig, 'Position', pos);
set(fig, 'PaperPositionMode', 'auto');
end

function applyPlotAspect(ax, aspectRatio)
if isempty(aspectRatio)
    return;
end
pbaspect(ax, [aspectRatio 1 1]);
end

function validateFigureSizeOption(figSize, what)
if isempty(figSize)
    return;
end
if ~isnumeric(figSize) || numel(figSize) ~= 2 || any(~isfinite(figSize)) || ...
        any(figSize <= 0)
    error('%s must be [] or [width height] with positive values.', what);
end
end

function validateAspectRatioOption(aspectRatio, what)
if isempty(aspectRatio)
    return;
end
if ~isscalar(aspectRatio) || ~isfinite(aspectRatio) || aspectRatio <= 0
    error('%s must be [] or a positive scalar (x:y ratio).', what);
end
end

function validateScalarTextOption(value, what)
if isempty(value)
    return;
end
if ~(ischar(value) || (isstring(value) && isscalar(value)))
    error('%s must be empty or a text scalar.', what);
end
end

function validatePerChannelTextOption(value, nChannels, what)
if isempty(value)
    return;
end
if ischar(value) || (isstring(value) && isscalar(value))
    return;
end
if ~iscell(value) || numel(value) ~= nChannels
    error('%s must be empty, a text scalar, or a cell array with %d entries.', ...
        what, nChannels);
end
for i = 1:nChannels
    v = value{i};
    if ~(isempty(v) || ischar(v) || (isstring(v) && isscalar(v)))
        error('%s{%d} must be empty or a text scalar.', what, i);
    end
end
end

function outText = resolveScalarText(value, defaultText, what)
if isempty(value)
    outText = defaultText;
    return;
end

if isstring(value)
    if ~isscalar(value)
        error('%s must be a text scalar.', what);
    end
    outText = char(value);
elseif ischar(value)
    outText = value;
else
    error('%s must be a text scalar.', what);
end
end

function outText = resolvePerChannelText(value, idx, defaultText, what, nChannels)
if isempty(value)
    outText = defaultText;
    return;
end

if ischar(value) || (isstring(value) && isscalar(value))
    outText = resolveScalarText(value, defaultText, what);
    return;
end

if ~iscell(value) || numel(value) ~= nChannels
    error('%s must be empty, text scalar, or cell array with %d entries.', ...
        what, nChannels);
end

entry = value{idx};
outText = resolveScalarText(entry, defaultText, ...
    sprintf('%s{%d}', what, idx));
end
