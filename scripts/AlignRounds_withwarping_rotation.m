% ---------------------------------------- script guide ----------------------------------------
% OVERVIEW
% This script aligns multi-round imaging data to a reference round (round 1).
% Note that the script will backup original image stack before doing alignment, so no need to manually backup and restore.
% Each round goes through up to 10 steps:
%   Step 1:  Global transform — rotation/translation/scaling via SURF feature
%            matching, intensity-based registration, or manual parameters
%   Step 2:  Polynomial warp — smooth global correction of residuals from step 1
%            (only when feature points are available)
%   Step 3:  2D translation — cross-correlation-based shift correction
%   Step 4:  Secondary refinement — another round of similarity/rigid transform
%            to catch residual scale/rotation after steps 1-3
%   Step 5:  PIV warping — local deformation correction at multiple scales
%   Step 6:  Save MIP overlay
%            ** mipOnly = true stops here (steps 7-10 skipped) **
%   Step 7:  Apply all transforms to the full DAPI z-stack
%   Step 8:  Z-shift — align z-stacks between rounds
%   Step 9:  Save cross-section overlay
%   Step 10: Apply all transforms to every channel, write aligned stacks & MIPs
%
% HOW TO USE
% 1. Set mipOnly = true (the pipline will only run until generating mip overlay, use it when tuning parameters)
% 2. Set rounds = 2:nrounds (or a subset, e.g. [4 7 8])
% 3. Choose parameters for each round (see below), run, inspect MIP overlays
% 4. Iterate: adjust parameters, re-run, until overlays look good
% 5. Set mipOnly = false, run to apply transforms to full z-stacks
%
% PARAMETERS — STEP 1 (global transform)
%   modeList{ri}:
%     'rigid'       — rotation + translation only (safest, no scaling)
%     'similarity'  — + uniform scaling
%     'affine'      — + anisotropic scaling
%     'projective'  — + perspective/tilt correction (most flexible, default)
%     Start with 'projective'. If it produces artifacts, try less flexible modes.
%   nFeaturePointsList(ri):
%     Number of SURF features to detect. Default 20000. Increase if matching
%     is poor, decrease if running out of memory.
%   rotationModels{ri}:
%     'FeaturePoints' — SURF-based (default, fast, usually best)
%     'Intensity'     — intensity-based registration (slow, use as fallback)
%   Manual override (last resort):
%     If automatic methods fail for a round, set parameters manually:
%       applyManualAngles(ri) = true;   rotationAngleList(ri) = <degrees>;
%       applyManualTranslation(ri) = true; translationList(ri,:) = [tx ty];
%       applyManualScale(ri) = true;    scaleList(ri) = <scale>;
%     Tip1: get rough angle from the preview images, then fine-tune translation.
%     Tip2: try increment/decrement rotation angle with a small amount each time to figure out how to set translation accordingly to prevent cutoffs 
%     Tip3: it maybe useful to disable step-3 (automatic 2D translation) when tuning numbers and uncomment once done
%       1. set disableAuto2Dtranslation = false
%       2. manually set angle, translation, and scale and examine mip overlay
%       3. once done, set disableAuto2Dtranslation = true
%
% PARAMETERS — STEP 2 (polynomial warp)
%   polyDegree: 2 (quadratic) or 3 (cubic). Usually does not have a big impact.
%   Skipped automatically when using intensity-based or manual alignment
%   (no feature points available). Also skipped if the fit fails.
%
% PARAMETERS — STEP 4 (secondary refinement)
%   secondaryModeList{ri}: same options as modeList. Default 'similarity'.
%   Uses the same method as step 1 (feature-point or intensity-based).
%   Skipped for rounds with manual alignment.
%   Check logged "residual scale" — should be close to 1.0 if step 1 worked well.
%
% PARAMETERS — STEP 5 (PIV warping)
%   EdgeLengthsList{ri}: vector of grid sizes for multi-scale PIV.
%     Larger values = broader/smoother correction. Smaller = more local.
%     Examples:
%       [200 100 50 10]  — default, good for most rounds
%       [100 50 10]      — if over-warped (too much local distortion)
%       [300 200 50 10]  — if under-warped (residual broad misalignment)
%       [200 200 50 10]  — repeating a scale can help convergence
%       []               — skip warping entirely for this round
%     If large edge lengths cause cutoff artifacts, use smaller values.
%   sigmaList(ri): smoothing parameter for PIV. Default 0.75.
%     Increase (e.g. 1.5) if warping introduces too much local distortion.
%
% PARAMETERS — STEP 8 (z-alignment)
%   zAlignWindow: [center_x_pct center_y_pct nx ny]
%     Region used to compute z-shift. Default [0.5 0.5 500 500] = center.
%   zTransformSize: block size for applying z-shift. Default [500 500].

clear; close all;
warning('off', 'MATLAB:nearlySingularMatrix');

scriptPath = pwd;
%all data folders are in a subfolder called 'data' of the directory containing the script
baseDir = fullfile(scriptPath);

% list folders containing image data
dirs = {
    '260225_ChickEmbryoGuojunx11_Rd1_SMAD23_TBXT_pERK'
    '260225_ChickEmbryoGuojunx11_Rd2_TBX6g_SOX2r_BCatm'
    '260226_ChickEmbryoGuojunx11_Rd3_SOX2r_FOXA2m_SOX17g'
    '260227_ChickEmbryoGuojunx11_Rd4_ISL1m_OTX2g_pAKTr'
    '260302_ChickEmbryoGuojunx11_Rd5_ECadm_SNAI1g_YAPr'
    '260303_ChickEmbryoGuojunx11_Rd6_SOX17g_PRDM1rat_GATA3r'
    '260304_ChickEmbryoGuojunx11_Rd7_FOXC2sheep_TFAP2Cm_CDX2r'
    '260305_ChickEmbryoGuojunx11_Rd8_SMAD23m_TBXTg_MIXL1r'
    '260306_ChickEmbryoGuojunx11_Rd9_PODXLm_VIMg_LEF1rSC'
};
dataDirs = fullfile(baseDir, dirs);
nrounds = length(dataDirs);

bare = 'stitched_p%.4d_w%.4d_t%.4d';
mipbare = 'stitched_MIP_p%.4d_w%.4d_t%.4d';

% load metadata from each round
channelLabel = cell(1, nrounds);
nucChannels = NaN(1, nrounds);
metas = cell(nrounds, 1);
for ri = 1:nrounds
    meta = load(fullfile(dataDirs{ri}, 'meta.mat'), 'meta');
    metas{ri} = meta.meta;
    channelLabel{ri} = metas{ri}.channelLabel;
    nucChannels(ri) = metas{ri}.nucChannel;
end
% Note: if npos is incorrect, set it manually (e.g. npos = 1;)
% npos = metas{1}.nPositions;
npos = 1;

% get file extension for MIPs
listing = dir(fullfile(dataDirs{1}, 'MIP', 'stitched_MIP_*'));
[~, ~, mipext] = fileparts(fullfile(listing(1).folder, listing(1).name));

% rounds to align (each new folder can be aligned as acquired)
% rounds = 2:nrounds;
rounds = 8:9;

% whether to do mip alignment only (useful for choosing parameter values before full run)
mipOnly = false;

% rotation parameters
% can choose between intensity based or feature point based, but feature point based usually work better
rotationModels = repelem({"FeaturePoints"}, nrounds);
% reduce this if memory is not enough
% increase if rotation/stretching/shrinking is inaccurate
nFeaturePointsList = repelem(20000, nrounds);
nFeaturePointsList(2) = 60000;
nFeaturePointsList(3) = 40000;
% whether to filter feature points, usually true
filterFeaturePoints = repelem(true, nrounds);
% mode for estimating rotation/translation/scaling
% three modes: rigid (least flexible), similarity (more flexible), affine (even more flexible), projective (most flexible)
% more flexible modes can lead to issues sometimes, if so, use less flexible modes
modeList = repelem({"projective"}, nrounds)'; 
modeList(2) = {"affine"};
modeList(7) = {"similarity"};
secondaryModeList = repelem({"similarity"}, nrounds)';

% polynomial warp parameters
% degree 2 = quadratic, degree 3 = cubic
% generally dont have huge impact, either 2 or 3 is fine
polyDegree = 3;

% warping parameters
% global warping switch
doWarping = true;
EdgeLengthsList = cell(1, nrounds);
sigmaList = repelem(0.75, nrounds);
for i = 1:nrounds
    EdgeLengthsList{i} = [200 100 50 10];
end
% sometimes different rounds might need different warping edge lengths
% for instance, to change alignemnt between round 1 and round 4, do EdgeLengthsList{4} = [100 50 10]
% large edge length gives more local stretch, so if over-stretched, do [100 50 10] or [50 10]
% if under-stretched, do [300 200 50 10] for instance, note that large edge length might lead to cutoffs
% if see cutoffs, try use smaller edge lengths, can also repeat same edge lengths multiple times, such as [200 200 50 10]
% set to [] to skip warping for that round
EdgeLengthsList{2} = [100 100 50 50];
sigmaList(2) = 1;
sigmaList(3) = 1;
EdgeLengthsList{4} = [100 50 10];
EdgeLengthsList{5} = [200 200 200 100];
EdgeLengthsList{6} = [200 200 100 50];
EdgeLengthsList{7} = [300 250 250 200 100 50];
sigmaList(7) = 1.5;
EdgeLengthsList{8} = [];
EdgeLengthsList{9} = [];

% manually adjust alignment
% sometiems if automatic alignment does not work, need to figure out manually
disableAuto2Dtranslation = false;
applyManualAngles = repelem(false, nrounds);
rotationAngleList = repelem(0, nrounds);
applyManualAngles(8) = true;
rotationAngleList(8) = -35;
applyManualAngles(9) = true;
rotationAngleList(9) = 57;
applyManualTranslation = repelem(false, nrounds);
translationList = zeros(nrounds, 2);
applyManualTranslation(8) = true;
translationList(8,:) = [-1250 1500];
applyManualTranslation(9) = true;
translationList(9,:) = [3750 -1250];
applyManualScale = repelem(false, nrounds);
scaleList = ones(nrounds, 2);
applyManualScale(8) = true;
scaleList(8) = 1.1;
applyManualScale(9) = true;
scaleList(9) = 1.07;

% z-align parameters
zAlignWindow = [0.5 0.5 500 500];   % [center_x_pct center_y_pct nx ny]
zTransformSize = [500 500];

% --- backup / restore original images before alignment ---
for ri = rounds
    backupDir = fullfile(dataDirs{ri}, 'stitched_backup');
    if ~exist(backupDir, 'dir')
        fprintf('Backing up round %d to %s\n', ri, backupDir);
        mkdir(backupDir);
        copyfile(fullfile(dataDirs{ri}, 'stitched_*.tif'), backupDir);
    else
        fprintf('Restoring round %d from %s\n', ri, backupDir);
        copyfile(fullfile(backupDir, 'stitched_*.tif'), dataDirs{ri});
    end
end

% iterate over positions
for pidx = 1:npos
    fprintf('position %d of %d\n', pidx, npos);

    % load reference DAPI stack (round 1)
    fname = fullfile(dataDirs{1}, [sprintf(bare, pidx-1, nucChannels(1), 0), '.tif']);
    fprintf('step 0/10: read in %s\n', fname);
    DAPI1 = loadTiffStack(fname);
    nz1 = size(DAPI1, 3);
    fprintf('step 0/10: image size (%dx%dx%d)\n', size(DAPI1, 1), size(DAPI1, 2), size(DAPI1, 3));
    mip1 = max(DAPI1, [], 3);
    mip1_adj = imadjust(mip1);
    lim1 = stitchedlim(mip1);

    for ri = rounds
        fprintf('round %d of %d\n', ri, nrounds);
        QCdir = fullfile(dataDirs{ri}, 'overlays');
        disp(QCdir);
        if ~exist(QCdir, 'dir'), mkdir(QCdir); end

        prefix = sprintf(bare, pidx-1, nucChannels(ri), 0);
        fname = fullfile(dataDirs{ri}, [prefix, '.tif']);
        fprintf('step 0/10: read in %s\n', fname);
        DAPI2 = loadTiffStack(fname);
        mip2 = max(DAPI2, [], 3);
        lim2 = stitchedlim(mip2);

        % pixel grid for interpolation
        m = size(mip2, 1);
        n = size(mip2, 2);
        [Xw, Yw] = meshgrid(1:n, 1:m);
        
        nFeaturePoints = nFeaturePointsList(ri);
        useIntensity = strcmp(rotationModels{ri}, "Intensity");

        % --- step 1/10: transform via feature matching or intensity ---
        fprintf('step 1/10: compute %s transform (position %d/%d, round %d/%d)\n', modeList{ri}, pidx, npos, ri, nrounds);
        referenceView = imref2d(size(mip1_adj));
        if applyManualAngles(ri) || applyManualTranslation(ri)
            tformRigid = rigidtform2d(0, [0 0]);
            tform = simtform2d(1, 0, [0 0]);
            if applyManualAngles(ri)
                tformRigid.RotationAngle = rotationAngleList(ri);
            end
            if applyManualTranslation(ri)
                tformRigid.Translation = translationList(ri,:);
            end
            if applyManualScale(ri)
                tform.Scale = scaleList(ri);
            end
            mip2 = imwarp(mip2, tformRigid, "OutputView", referenceView);
            inlierIdx = [];
        elseif ~useIntensity && ~(applyManualAngles(ri) || applyManualTranslation(ri))
            mip1_tmp = mip1_adj;
            mip2_tmp = imadjust(mip2);
            ptsOriginal  = detectSURFFeatures(mip1_tmp);
            ptsDistorted = detectSURFFeatures(mip2_tmp);
            if filterFeaturePoints(ri)
                ptsOriginal  = ptsOriginal.selectStrongest(nFeaturePoints);
                ptsDistorted = ptsDistorted.selectStrongest(nFeaturePoints);
            end
            [featuresOriginal,  validPtsOriginal] = extractFeatures(mip1_tmp, ptsOriginal);
            [featuresDistorted, validPtsDistorted] = extractFeatures(mip2_tmp, ptsDistorted);
            clear mip1_tmp mip2_tmp;
            index_pairs = matchFeatures(featuresOriginal, featuresDistorted);
            matchedPtsOriginal  = validPtsOriginal(index_pairs(:,1));
            matchedPtsDistorted = validPtsDistorted(index_pairs(:,2));
            rng(1);
            [tform, inlierIdx] = estgeotform2d(matchedPtsDistorted, matchedPtsOriginal, modeList{ri});
            fprintf('  matched features: %d, inliers: %d\n', size(index_pairs,1), sum(inlierIdx));
        else
            [optimizer, metric] = imregconfig('monomodal');
            optimizer.MaximumIterations = 300;
            tform = imregtform(mip2, mip1, modeList{ri}, optimizer, metric, 'PyramidLevels', 5);
            inlierIdx = [];
            fprintf('  intensity-based %s registration\n', modeList{ri});
        end
        mipa = imwarp(mip2, tform, "OutputView", referenceView);
        clear mip2;
        fprintf('step 1/10: transform done (position %d/%d, round %d/%d)\n', pidx, npos, ri, nrounds);

        % --- step 2/10: polynomial warp from feature-point residuals ---
        doPolyWarp = false;
        fprintf('step 2/10: compute polynomial warp degree %d (position %d/%d, round %d/%d)\n', polyDegree, pidx, npos, ri, nrounds);
        if ~isempty(inlierIdx) && ~applyManualAngles(ri) && ~applyManualTranslation(ri)
            inlierPtsOrig = matchedPtsOriginal(inlierIdx);
            inlierPtsDist = matchedPtsDistorted(inlierIdx);
            warpedPts = transformPointsForward(tform, inlierPtsDist.Location);
            try
                tformPoly = fitgeotrans(warpedPts, inlierPtsOrig.Location, 'polynomial', polyDegree);
                mipa = imwarp(mipa, tformPoly, "OutputView", referenceView);
                doPolyWarp = true;
                fprintf('step 2/10: polynomial warp done (position %d/%d, round %d/%d)\n', pidx, npos, ri, nrounds);
            catch
                disp("polynomial warping failed, proceed to next step (not a big problem)")
            end
        else
            disp("no feature points, skip polynomial warping (not a big problem)")
        end

        % --- step 3/10: 2D translation ---
        if ~disableAuto2Dtranslation
            fprintf('step 3/10: compute 2D translation (position %d/%d, round %d/%d)\n', pidx, npos, ri, nrounds);
            shiftyx = findImageShift(mip1, mipa, 'automatic');
            mipa = alignImage(mip1, mipa, shiftyx);
            fprintf('  shift: [%.1f, %.1f] px\n', shiftyx(1), shiftyx(2));
            fprintf('step 3/10: 2D translation done (position %d/%d, round %d/%d)\n', pidx, npos, ri, nrounds);
        else
            fprintf('skip step 3/10: 2D translation (position %d/%d, round %d/%d)\n', pidx, npos, ri, nrounds);
        end

        % --- step 4/10: refinement (residual scale/rotation) ---
        doSimWarp = false;
        fprintf('step 4/10: compute %s refinement (position %d/%d, round %d/%d)\n', secondaryModeList{ri}, pidx, npos, ri, nrounds);
        if ~applyManualAngles(ri) && ~applyManualTranslation(ri)
            if ~useIntensity
                % feature-point based refinement
                mip2_tmp = imadjust(mipa);
                ptsOrig2  = detectSURFFeatures(mip1_adj);
                ptsDist2  = detectSURFFeatures(mip2_tmp);
                if filterFeaturePoints(ri)
                    ptsOrig2 = ptsOrig2.selectStrongest(nFeaturePoints);
                    ptsDist2 = ptsDist2.selectStrongest(nFeaturePoints);
                end
                [fOrig2, vpOrig2] = extractFeatures(mip1_adj, ptsOrig2);
                [fDist2, vpDist2] = extractFeatures(mip2_tmp, ptsDist2);
                clear mip2_tmp;
                try
                    idx2 = matchFeatures(fOrig2, fDist2);
                    rng(1);
                    [tformSim, ~] = estgeotform2d(vpDist2(idx2(:,2)), vpOrig2(idx2(:,1)), secondaryModeList{ri});
                    fprintf('  residual scale: %.6f, rotation: %.4f deg\n', tformSim.Scale, tformSim.RotationAngle);
                    mipa = imwarp(mipa, tformSim, "OutputView", referenceView);
                    doSimWarp = true;
                    fprintf('step 4/10: refinement done (position %d/%d, round %d/%d)\n', pidx, npos, ri, nrounds);
                catch
                    disp("feature-point refinement failed, proceed to next step (not a big problem)")
                end
            else
                % intensity-based refinement
                try
                    [optimizer, metric] = imregconfig('monomodal');
                    optimizer.MaximumIterations = 300;
                    tformSim = imregtform(mipa, mip1, secondaryModeList{ri}, optimizer, metric);
                    mipa = imwarp(mipa, tformSim, "OutputView", referenceView);
                    doSimWarp = true;
                    fprintf('  intensity-based %s refinement\n', secondaryModeList{ri});
                    fprintf('step 4/10: refinement done (position %d/%d, round %d/%d)\n', pidx, npos, ri, nrounds);
                catch
                    disp("intensity-based refinement failed, proceed to next step (not a big problem)")
                end
            end
        end

        % --- step 5/10: PIV warping ---
        if doWarping
            EdgeLengths = EdgeLengthsList{ri};
            fprintf('step 5/10: compute 2D warping (position %d/%d, round %d/%d)\n', pidx, npos, ri, nrounds);
            nEdgeLengths = length(EdgeLengths);
            VXinterps = zeros([size(Xw), nEdgeLengths], 'single');
            VYinterps = zeros([size(Xw), nEdgeLengths], 'single');
            mipw = mipa;
            for i = 1:nEdgeLengths
                fprintf('> warp with edge length: %d\n', EdgeLengths(i));
                [X, Y, VX, VY, ~] = GetPIV3(mip1, mipw, EdgeLengths(i), sigmaList(ri));
                VXinterps(:,:,i) = interp2(X, Y, VX, Xw, Yw);
                VYinterps(:,:,i) = interp2(X, Y, VY, Xw, Yw);
                mipw = uint16(interp2(Xw, Yw, single(mipw), ...
                    Xw + VXinterps(:,:,i), Yw + VYinterps(:,:,i), 'linear'));
            end
            fprintf('step 5/10: 2D warping done (position %d/%d, round %d/%d)\n', pidx, npos, ri, nrounds);
        else
            mipw = mipa;
        end
        clear mipa;

        % --- step 6/10: save MIP overlay ---
        fprintf('step 6/10: save mip overlay (position %d/%d, round %d/%d)\n', pidx, npos, ri, nrounds);
        overlay = makeMIPOverlay(mip1, mipw, lim1, lim2, 0.075, 1, ri);
        savename = fullfile(QCdir, [prefix, '_MIPoverlay.png']);
        imwrite(overlay, savename);
        clear overlay mipw;
        fprintf('step 6/10: mip overlay saved (position %d/%d, round %d/%d)\n', pidx, npos, ri, nrounds);

        if mipOnly
            continue;
        end

        % --- step 7/10: transform DAPI z-stack ---
        fprintf('step 7/10: transform dapi z-stack (position %d/%d, round %d/%d)\n', pidx, npos, ri, nrounds);
        m = size(DAPI2, 1);
        n = size(DAPI2, 2);
        nz = size(DAPI2, 3);
        imclass = class(DAPI2);

        dapia = zeros(m, n, nz, imclass);
        for zi = 1:nz
            slice = imwarp(DAPI2(:,:,zi), tform, "OutputView", referenceView);
            if doPolyWarp
                slice = imwarp(slice, tformPoly, "OutputView", referenceView);
            end
            slice = alignImage(slice, slice, shiftyx);
            if doSimWarp
                slice = imwarp(slice, tformSim, "OutputView", referenceView);
            end
            dapia(:,:,zi) = slice;
        end
        clear DAPI2 slice;

        if doWarping
            for i = 1:size(VXinterps, 3)
                dapia = xyalignImageStack(dapia, [0 0], Xw, Yw, VXinterps(:,:,i), VYinterps(:,:,i));
            end
        end
        fprintf('step 7/10: dapi z-stack transformed (position %d/%d, round %d/%d)\n', pidx, npos, ri, nrounds);

        % --- step 8/10: z-shift ---
        fprintf('step 8/10: compute z-shift (position %d/%d, round %d/%d)\n', pidx, npos, ri, nrounds);
        xMin = max(floor(n * zAlignWindow(1)) - zAlignWindow(3), 1);
        yMin = max(floor(m * zAlignWindow(2)) - zAlignWindow(4), 1);
        xMax = min(xMin + zAlignWindow(3), n);
        yMax = min(yMin + zAlignWindow(4), m);
        [shiftz, scalez] = zAlignImageStacks(DAPI1(yMin:yMax, xMin:xMax, :), dapia(yMin:yMax, xMin:xMax, :));
        if scalez < 0.5 || scalez > 2.0 || isnan(scalez)
            fprintf('  WARNING: z-alignment failed (shift=%.2f, scale=%.4f), using defaults\n', shiftz, scalez);
            shiftz = 0;
            scalez = 1;
        end
        fprintf('  z-shift: %.2f, z-scale: %.4f\n', shiftz, scalez);

        % apply z shift
        dapia_new = zeros(m, n, nz1, imclass);
        for j = 1:zTransformSize(2):m
            for i = 1:zTransformSize(1):n
                xMin_b = i;
                xMax_b = min(i + zTransformSize(1) - 1, n);
                yMin_b = j;
                yMax_b = min(j + zTransformSize(2) - 1, m);
                dapia_new(yMin_b:yMax_b, xMin_b:xMax_b, :) = zShiftScale(dapia(yMin_b:yMax_b, xMin_b:xMax_b, :), shiftz, scalez, nz1);
            end
        end
        dapia = dapia_new;
        clear dapia_new;
        fprintf('step 8/10: z-shift done (position %d/%d, round %d/%d)\n', pidx, npos, ri, nrounds);

        % --- step 9/10: save cross-section overlay ---
        fprintf('step 9/10: save cross-section overlay (position %d/%d, round %d/%d)\n', pidx, npos, ri, nrounds);
        index = round(size(dapia, 1) / 2);
        crossSectionOverlay(DAPI1, dapia, metas{1}, index, lim1, lim2, 0.3, 1, ri);
        savename = fullfile(QCdir, [prefix, sprintf('_crossSectionX%d_overlay.png', index)]);
        saveas(gcf, savename);
        exportgraphics(gcf, savename, 'Resolution', 300);
        fprintf('step 9/10: cross-section saved (position %d/%d, round %d/%d)\n', pidx, npos, ri, nrounds);

        % --- step 10/10: transform all channels ---
        fprintf('step 10/10: transform all channels (position %d/%d, round %d/%d)\n', pidx, npos, ri, nrounds);
        nchannels = metas{ri}.nChannels;

        for ci = 1:nchannels
            fprintf('> transform channel %d\n', ci);
            fname = fullfile(dataDirs{ri}, [sprintf(bare, pidx-1, ci-1, 0), '.tif']);

            if ci - 1 == nucChannels(ri)
                % nuclear channel is already aligned — write directly
                writeTiffZStack(dapia, fname);
                MIP = max(dapia, [], 3);
                clear dapia;   % <-- free immediately after writing
            else
                img = loadTiffStack(fname);
                m = size(img, 1);
                n = size(img, 2);
                nz = size(img, 3);
                imclass = class(img);

                % projective + polynomial + shift + similarity
                for zi = 1:nz
                    slice = imwarp(img(:,:,zi), tform, "OutputView", referenceView);
                    if doPolyWarp
                        slice = imwarp(slice, tformPoly, "OutputView", referenceView);
                    end
                    slice = alignImage(slice, slice, shiftyx);
                    if doSimWarp
                        slice = imwarp(slice, tformSim, "OutputView", referenceView);
                    end
                    img(:,:,zi) = slice;
                end
                clear slice;

                % PIV warping
                if doWarping
                    for i = 1:size(VXinterps, 3)
                        img = xyalignImageStack(img, [0 0], Xw, Yw, VXinterps(:,:,i), VYinterps(:,:,i));
                    end
                end

                % z shift
                img_new = zeros(m, n, nz1, imclass);
                for j = 1:zTransformSize(2):m
                    for i = 1:zTransformSize(1):n
                        xMin_b = i;
                        xMax_b = min(i + zTransformSize(1) - 1, n);
                        yMin_b = j;
                        yMax_b = min(j + zTransformSize(2) - 1, m);

                        img_new(yMin_b:yMax_b, xMin_b:xMax_b, :) = ...
                            zShiftScale(img(yMin_b:yMax_b, xMin_b:xMax_b, :), shiftz, scalez, nz1);
                    end
                end
                img = img_new;
                clear img_new;
                writeTiffZStack(img, fname);
                MIP = max(img, [], 3);
                clear img;
            end

            % write MIP
            fprintf('> channel %d transformed\n', ci);
            mipname = fullfile(dataDirs{ri}, 'MIP', [sprintf(mipbare, pidx-1, ci-1, 0), mipext]);
            if strcmp(mipext, '.jpg')
                imwrite(im2double(MIP), mipname, 'Quality', 99);
            else
                imwrite(MIP, mipname);
            end
            clear MIP;
            fprintf('> channel %d saved\n', ci);
        end

        clear VXinterps VYinterps tformPoly tformSim;
        fprintf('step 10/10: all channels transformed and saved (position %d/%d, round %d/%d)\n', pidx, npos, ri, nrounds);
    end
    close all;
end

function IMw = xyalignImageStack(img, shiftyx, Xw, Yw, VXinterp, VYinterp)
    m = size(img, 1);
    n = size(img, 2);
    nz = size(img, 3);
    IMw = zeros(m, n, nz, 'uint16');
    parfor zi = 1:nz
        im = img(:,:,zi);
        ima = alignImage(im, im, shiftyx);
        imw = uint16(interp2(Xw, Yw, single(ima), Xw + VXinterp, Yw + VYinterp, 'cubic'));
        IMw(:,:,zi) = imw;
    end
end

function RGB = makeMIPOverlay(mip1, mip2, lim1, lim2, cfs, rid1, rid2)
    m = size(mip1, 1);
    n = size(mip1, 2);
    A = imadjust(mip1, lim1);
    B = imadjust(mip2, lim2);
    overlay = cat(3, A, B, A);
    text = {sprintf('round %d', rid1), sprintf('round %d', rid2)};
    ypos = round([m*(1 - 1.5*cfs); m*(1 - 1.5*cfs)]);
    xpos = round([0.005*n; 0.5*n]);
    RGB = insertText(overlay, [xpos, ypos], text, ...
        'BoxColor', 'black', 'TextColor', {'m','g'}, 'BoxOpacity', 0.3, 'FontSize', 96);
end

function crossSectionOverlay(img1, img2, meta1, index, lim1, lim2, cfs, rid1, rid2)
    xres = meta1.xres;
    zres = meta1.zres;

    cross1 = transpose(squeeze(img1(index,:,end:-1:1)));
    zsize1 = round(size(cross1, 1) * zres / xres);
    cross1 = imresize(cross1, [zsize1, size(cross1, 2)]);

    cross2 = transpose(squeeze(img2(index,:,end:-1:1)));
    zsize2 = round(size(cross2, 1) * zres / xres);
    cross2 = imresize(cross2, [zsize2, size(cross2, 2)]);

    if zsize1 > zsize2
        cross2 = [zeros(zsize1 - zsize2, size(cross2, 2), 'uint16'); cross2];
    elseif zsize2 > zsize1
        cross1 = [zeros(zsize2 - zsize1, size(cross1, 2), 'uint16'); cross1];
    end
    zsize = size(cross1, 1);
    n = size(cross1, 2);

    C1 = imadjust(cross1, lim1);
    C2 = imadjust(cross2, lim2);
    o1 = cat(3, C1, C2, C1);

    titles = {sprintf('round %d', rid1), sprintf('round %d', rid2), ...
        "\color{magenta}round 1 \color{green} round 2"};
    crosses = {C1, C2, o1};

    figpos = figurePosition([n, 2*zsize*length(crosses)]);
    figure('Position', figpos);
    for ii = 1:length(crosses)
        subplot_tight(length(crosses), 1, ii);
        imshow(crosses{ii});
        cleanSubplot;
        ypos = zsize * 0.5 * cfs;
        xpos = 0.01 * n * cfs;
        text(xpos, ypos, titles{ii}, ...
            'Color', 'w', 'FontUnits', 'normalized', 'FontSize', cfs, ...
            'FontWeight', 'bold');
    end
end

function writeTiffZStack(img, fname)
    % estimate file size in bytes
    fileSize = numel(img) * 2;  % uint16 = 2 bytes
    useBigTiff = fileSize > 3.5e9;  % switch at 3.5GB
    
    if useBigTiff
        t = Tiff(fname, 'w8');
        for zi = 1:size(img, 3)
            t.setTag('ImageLength', size(img, 1));
            t.setTag('ImageWidth', size(img, 2));
            t.setTag('Photometric', Tiff.Photometric.MinIsBlack);
            t.setTag('BitsPerSample', 16);
            t.setTag('SamplesPerPixel', 1);
            t.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
            t.setTag('Compression', Tiff.Compression.None);
            t.setTag('SampleFormat', Tiff.SampleFormat.UInt);
            t.write(img(:,:,zi));
            if zi < size(img, 3)
                t.writeDirectory();
            end
        end
        t.close();
    else
        for zi = 1:size(img, 3)
            if zi == 1
                mode = 'overwrite';
            else
                mode = 'append';
            end
            imwrite(img(:,:,zi), fname, 'WriteMode', mode);
        end
    end
end
