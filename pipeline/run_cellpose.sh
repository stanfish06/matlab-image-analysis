#!/usr/bin/env bash
#SBATCH --job-name=cellpose
#SBATCH --output=cellpose_%j.out
#SBATCH --error=cellpose_%j.err
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=zyyu@umich.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --account=iheemske0
# GPU partition: gpu or gpu_mig40
# CPU partition: standard or largemem
#SBATCH --partition=gpu
#SBATCH --gpus=1

#=== Cellpose parameters ===
# Nuclear channel index (matches nucChannel in preprocess)
NUC_CHANNEL=0
# Pretrained model name or path
CELLPOSE_MODEL=olympus_scope_dapi_20x_epoch_2800
# Estimated cell diameter in pixels (0 = auto-estimate)
DIAMETER=26
# Run 3D segmentation (set to false for 2D)
DO_3D=true
# Anisotropy ratio (z-spacing / xy-spacing), only used with 3D
ANISOTROPY=6
# Cell probability threshold (-6 to 6, lower = more permissive)
CELLPROB_THRESHOLD=-1
# Flow error threshold (0-3, higher = more permissive)
FLOW_THRESHOLD=0.5
# Use GPU acceleration
USE_GPU=true
# Save output as tif (in addition to npy)
SAVE_TIF=true
# Force re-copy TIFs into tmp even if tmp already exists
FORCE_RECOPY=true

#=== Environment ===
CELLPOSE_VENV=~/cellpose-gpu-stan/bin/activate

#=== Setup ===
module load python/3.12
source "$CELLPOSE_VENV"

WORK_DIR=$(pwd)
# Match the nuclear-channel images by suffix only, so this works for stitched
# (stitched_p####_w####_t0000.tif) AND non-stitched z-stacks
# (<rawname>_p####_w####_t0000.tif). Masks/segmentations end in _cp_masks.tif /
# _FinalSegmentation.tif so they are not matched. Run on stitched OR
# non-stitched data (preprocess_images does one or the other), not both at once.
NUC_PATTERN=$(printf "*_w%04d_t0000.tif" "$NUC_CHANNEL")
TMP_DIR="${WORK_DIR}/cellpose_tmp"

# create tmp dir and copy nuclear channel TIFs
if [ -d "$TMP_DIR" ] && [ "$FORCE_RECOPY" = false ]; then
    echo "cellpose_tmp already exists, reusing (set FORCE_RECOPY=true to re-copy)"
else
    mkdir -p "$TMP_DIR"
    echo "Copying nuclear channel TIFs to cellpose_tmp"
    cp ${WORK_DIR}/${NUC_PATTERN} "$TMP_DIR"/
fi

N_FILES=$(ls "$TMP_DIR"/${NUC_PATTERN} 2>/dev/null | wc -l)
echo "Found $N_FILES nuclear channel TIFs (w$(printf '%04d' "$NUC_CHANNEL"))"

if [ "$N_FILES" -eq 0 ]; then
    echo "ERROR: no files matching ${NUC_PATTERN} in cellpose_tmp"
    exit 1
fi

#=== Build command ===
CMD="cellpose --dir $TMP_DIR --pretrained_model $CELLPOSE_MODEL --diameter $DIAMETER --cellprob_threshold $CELLPROB_THRESHOLD --flow_threshold $FLOW_THRESHOLD --verbose"

if [ "$USE_GPU" = true ]; then
    CMD="$CMD --use_gpu"
fi

if [ "$DO_3D" = true ]; then
    CMD="$CMD --do_3D --anisotropy $ANISOTROPY"
fi

if [ "$SAVE_TIF" = true ]; then
    CMD="$CMD --save_tif"
fi

echo "Running: $CMD"
eval "$CMD"

#=== Copy masks back ===
echo "Copying masks to ${WORK_DIR}"
cp "$TMP_DIR"/*_cp_masks* "$WORK_DIR"/ 2>/dev/null
cp "$TMP_DIR"/*_seg.npy "$WORK_DIR"/ 2>/dev/null

echo "Done. Masks in ${WORK_DIR}"
