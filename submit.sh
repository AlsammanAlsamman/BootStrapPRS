#!/bin/bash
#SBATCH --job-name=gwas_finemap
#SBATCH --output=logs/gwas_finemap_%j.out
#SBATCH --error=logs/gwas_finemap_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

# Create logs directory
mkdir -p logs

# Parse command line arguments
DRY_RUN=""
TARGETS=()
EXTRA_FLAGS=""
SNAKEFILE=""
CORES=""
JOBS=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run|-n)
            DRY_RUN="--dry-run"
            shift
            ;;
        --snakefile)
            SNAKEFILE="$2"
            shift 2
            ;;
        --cores)
            CORES="$2"
            shift 2
            ;;
        --jobs|-j)
            JOBS="$2"
            shift 2
            ;;
        --*)
            # Pass through other flags
            EXTRA_FLAGS="$EXTRA_FLAGS $1"
            if [[ $2 && ! $2 =~ ^-- ]]; then
                EXTRA_FLAGS="$EXTRA_FLAGS $2"
                shift 2
            else
                shift
            fi
            ;;
        *)
            # This is a target file
            TARGETS+=("$1")
            shift
            ;;
    esac
done

# Detect if this is an individual rule submission
INDIVIDUAL_RULE=false
if [[ -n "$SNAKEFILE" && "$SNAKEFILE" =~ ^rules/ ]]; then
    INDIVIDUAL_RULE=true
fi

# Set defaults based on submission type
if [[ "$INDIVIDUAL_RULE" == true ]]; then
    # Individual rule defaults
    DEFAULT_SNAKEFILE="$SNAKEFILE"
    DEFAULT_CORES="${CORES:-2}"
    DEFAULT_JOBS="${JOBS:-1}"
    DEFAULT_MEM="32G"
    DEFAULT_TIME="2:00:00"
    DEFAULT_CPUS=2
else
    # Full pipeline defaults
    DEFAULT_SNAKEFILE="${SNAKEFILE:-Snakefile}"
    DEFAULT_CORES="${CORES:-8}"
    DEFAULT_JOBS="${JOBS:-15}"
    DEFAULT_MEM="64G"
    DEFAULT_TIME="2-00:00:00"
    DEFAULT_CPUS=8
fi

# Default to 'all' if no targets specified and not individual rule
if [ ${#TARGETS[@]} -eq 0 ] && [[ "$INDIVIDUAL_RULE" == false ]]; then
    TARGETS=("all")
fi

# Set up environment
module load slurm
module load python/3.7.0

# Export Python path if needed
export PYTHONPATH="/s/nath-lab/alsamman/____MyCodes____/FineMappingSuite:$PYTHONPATH"

# Print submission info
echo "Submission type: $([ "$INDIVIDUAL_RULE" == true ] && echo "Individual Rule" || echo "Full Pipeline")"
echo "Snakefile: $DEFAULT_SNAKEFILE"
echo "Cores: $DEFAULT_CORES"
echo "Target(s): ${TARGETS[*]:-"(from snakefile)"}"
echo "Jobs: $DEFAULT_JOBS"
echo ""

# Adjust SLURM header for individual rules
if [[ "$INDIVIDUAL_RULE" == true ]]; then
    echo "Adjusting SLURM parameters for individual rule execution..."
fi

# Submit the workflow to SLURM
if [[ "$INDIVIDUAL_RULE" == true ]]; then
    # Individual rule submission - simpler configuration
    snakemake \
        --cluster "sbatch \
            --mem=$DEFAULT_MEM \
            --time=$DEFAULT_TIME \
            --job-name=rule_{rule}_{wildcards} \
            --cpus-per-task=$DEFAULT_CPUS \
            --output=logs/{rule}_{wildcards}_%j.out \
            --error=logs/{rule}_{wildcards}_%j.err" \
        --jobs $DEFAULT_JOBS \
        --latency-wait 60 \
        --keep-going \
        --rerun-incomplete \
        --snakefile "$DEFAULT_SNAKEFILE" \
        --cores $DEFAULT_CORES \
        --verbose \
        --printshellcmds \
        --stats logs/snakemake_individual_stats.json \
        $DRY_RUN \
        $EXTRA_FLAGS \
        "${TARGETS[@]}"
else
    # Full pipeline submission - original configuration
    snakemake \
        --cluster "sbatch \
            --mem={resources.mem_mb} \
            --time={resources.time} \
            --job-name={rule}_{wildcards} \
            --cpus-per-task={threads} \
            --output=logs/{rule}_{wildcards}_%j.out \
            --error=logs/{rule}_{wildcards}_%j.err" \
        --jobs $DEFAULT_JOBS \
        --latency-wait 120 \
        --keep-going \
        --rerun-incomplete \
        --configfile config/pipeline_config.yaml \
        --snakefile "$DEFAULT_SNAKEFILE" \
        --cores $DEFAULT_CORES \
        --verbose \
        --printshellcmds \
        --stats logs/snakemake_stats.json \
        $DRY_RUN \
        $EXTRA_FLAGS \
        "${TARGETS[@]}"
fi

# Print completion message
echo "Snakemake workflow submitted for target(s): ${TARGETS[*]}"
echo "Monitor jobs with: squeue -u $USER"