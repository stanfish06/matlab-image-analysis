function run_pipeline(start_from, opts)
%RUN_PIPELINE Master pipeline entry point.
%   RUN_PIPELINE() runs the full pipeline from step 0.
%   RUN_PIPELINE(START_FROM) skips completed steps and starts at START_FROM.
%   RUN_PIPELINE(START_FROM, OPTS) accepts options:
%       OPTS.dryRun - if true, report planned steps and missing inputs
%                     without running scripts or submitting SLURM jobs.
%
%   Run from inside a plate directory after preprocessing is complete.
%
% Run from inside a plate directory (after preprocess is done).
% Set start_from to skip completed steps.
%
% Steps:
%   0 = run_cellpose      (submit SLURM job for cellpose segmentation)
%   1 = clean_masks       (filter cellpose masks, save FinalSegmentation + masks.mat)
%   2 = quantify          (single-cell quantification -> positions.mat)
%   3 = extract_colonies  (colony picking -> positions_colony.mat)
%   4 = make_scatter      (scatter plots + population stats)
%   5 = make_montage      (montage visualizations)
%
% Usage:
%   run_pipeline(0);  % run full pipeline (including cellpose)
%   run_pipeline(1);  % skip cellpose, start from clean masks
%   run_pipeline(4);  % re-run from scatter (e.g. after adjusting thresholds)
%   run_pipeline(0, struct('dryRun', true)); % preflight without running

if nargin < 1 || isempty(start_from)
    start_from = 0;
end
if nargin < 2 || isempty(opts)
    opts = struct();
end
if ~isfield(opts, 'dryRun')
    opts.dryRun = false;
end

steps = {
    'run_cellpose'
    'clean_masks'
    'quantify'
    'extract_colonies'
    'make_scatter'
    'make_montage'
};

step_labels = {
    'Cellpose segmentation (SLURM)'
    'Clean cellpose masks'
    'Single-cell quantification'
    'Extract colonies'
    'Scatter plots'
    'Montage visualization'
};

% expected outputs per step (used to validate before running next step)
expected_outputs = {
    {'*_cp_masks.tif'}
    {'*_FinalSegmentation.tif', '*_masks.mat'}
    {'positions.mat', 'meta_combined.mat'}
    {'positions_colony.mat', 'meta_colony.mat'}
    {fullfile('scatter', 'scatterPercent.csv')}
    {'mip_montage'}
};

% expected inputs per step (checked before running)
expected_inputs = {
    {'meta.mat', '*_w*_t0000.tif'}
    {'*_cp_masks.tif', 'meta.mat'}
    {'*_masks.mat', 'meta.mat'}
    {'positions.mat', 'meta.mat'}
    {'positions_colony.mat', 'meta_colony.mat'}
    {'meta.mat', '*_w*_t0000.tif'}
};

%% run pipeline
% step indices: 0-based in the UI, 1-based in arrays
for i = (start_from + 1):length(steps)
    step_num = i - 1; % display as 0-based
    fprintf('\n========================================\n');
    fprintf('Step %d/%d: %s\n', step_num, length(steps) - 1, step_labels{i});
    fprintf('========================================\n');

    % check required inputs exist
    missing = check_files(expected_inputs{i});
    if ~isempty(missing)
        msg = sprintf('Missing required inputs for step %d (%s):\n  %s\nDid a previous step complete?', ...
            step_num, steps{i}, strjoin(missing, '\n  '));
        if opts.dryRun
            fprintf('%s\n', msg);
        else
            error('run_pipeline:MissingInputs', '%s', msg);
        end
    end

    if opts.dryRun
        fprintf('[dry run] Would run step %d (%s).\n', step_num, steps{i});
        continue
    end

    % cellpose is a SLURM job, handle separately
    if step_num == 0
        run_cellpose_step();
        continue
    end

    % run MATLAB step in isolated scope
    tic;
    run_step(steps{i});
    elapsed = toc;
    fprintf('Step %d complete (%.1f sec).\n', step_num, elapsed);

    % verify outputs were created
    missing = check_files(expected_outputs{i});
    if ~isempty(missing)
        warning('Expected outputs not found after step %d:\n  %s', step_num, strjoin(missing, '\n  '));
    end
end

fprintf('\n========================================\n');
if opts.dryRun
    fprintf('Pipeline dry run complete.\n');
else
    fprintf('Pipeline complete.\n');
end
fprintf('========================================\n');

end

%% helper functions
function run_step(script_name)
    feval(script_name, pwd);
end

function run_cellpose_step()
    % submit cellpose as a SLURM job and wait for it to finish
    [status, out] = system('sbatch run_cellpose.sh');
    if status ~= 0
        error('Failed to submit cellpose job:\n%s', out);
    end
    fprintf('%s', out);
    % extract job ID from "Submitted batch job 12345"
    tokens = regexp(out, 'Submitted batch job (\d+)', 'tokens');
    if isempty(tokens)
        error('Could not parse SLURM job ID from: %s', out);
    end
    job_id = tokens{1}{1};
    fprintf('Waiting for SLURM job %s to complete...\n', job_id);

    % poll until job is no longer in the queue
    while true
        [~, squeue_out] = system(sprintf('squeue -j %s -h 2>/dev/null', job_id));
        if isempty(strtrim(squeue_out))
            break
        end
        fprintf('.');
        pause(30);
    end
    fprintf('\nCellpose job %s finished.\n', job_id);

    % check if it succeeded by looking for output masks
    masks = dir('*_cp_masks.tif');
    if isempty(masks)
        error('Cellpose job completed but no *_cp_masks.tif found. Check cellpose_*.err for errors.');
    end
    fprintf('Found %d mask files.\n', length(masks));
end

function missing = check_files(patterns)
    missing = cell(size(patterns));
    nMissing = 0;
    for j = 1:numel(patterns)
        pat = patterns{j};
        if contains(pat, '*')
            % glob pattern
            matches = dir(pat);
            if isempty(matches)
                nMissing = nMissing + 1;
                missing{nMissing} = pat;
            end
        else
            % exact path
            if ~exist(pat, 'file') && ~exist(pat, 'dir')
                nMissing = nMissing + 1;
                missing{nMissing} = pat;
            end
        end
    end
    missing = missing(1:nMissing);
end
