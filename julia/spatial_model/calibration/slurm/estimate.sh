#!/bin/bash
# TikTak estimation of the spatial block on one node (spatial_estimate.jl), with
# SLURM_CPUS_PER_TASK − 1 worker processes. Submit from the repository root; extra
# arguments pass through to spatial_estimate.jl, e.g.
#   sbatch julia/spatial_model/calibration/slurm/estimate.sh --n-samples 512 --max-evals 6000
#   for r in 0 0.05 0.10 0.1813 0.30; do
#       sbatch -J sp-rm$r julia/spatial_model/calibration/slurm/estimate.sh --fix r_m=$r; done
# The soft deadline --max-seconds is set 30 minutes inside the wall time; resume an
# interrupted run with --run-dir <dir> --resume. Partitions preempt by priority tier
# (econ-grad 220 > econ 200 > sscc, short 100) and requeue the job; a requeued job
# resumes its run from the TikTak journal when --run-dir names a directory holding
# one. Do not edit the model files while a run is queued: the problem id hashes them,
# and a changed id refuses to resume.
#SBATCH -J sp-estimate
#SBATCH -p econ-grad
#SBATCH -t 24:00:00
#SBATCH -c 64
#SBATCH --mem-per-cpu 2G
#SBATCH --export=ALL
#SBATCH -o julia/spatial_model/calibration/slurm/logs/%x-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
wall=$(squeue -h -j "$SLURM_JOB_ID" -o %L | awk -F'[-:]' '{n=NF; s=$n+60*$(n-1); if(n>2)s+=3600*$(n-2); if(n>3)s+=86400*$(n-3); print s}')
# After a requeue, resume only a run that exists: the --run-dir from the arguments
# must hold a journal (a job preempted before TikTak started simply starts again).
resume=()
if [ "${SLURM_RESTART_COUNT:-0}" -gt 0 ]; then
    run_dir=""
    args=("$@")
    for ((i = 0; i < ${#args[@]} - 1; i++)); do
        [ "${args[i]}" = "--run-dir" ] && run_dir="${args[i+1]}"
    done
    if [ -n "$run_dir" ] && [ -f "$run_dir/history.jsonl" ] && [[ " $* " != *" --resume "* ]]; then
        resume=(--resume)
    fi
fi
exec julia --project=julia/spatial_model/calibration/estimation julia/spatial_model/calibration/spatial_estimate.jl \
    --max-seconds $(( wall > 3600 ? wall - 1800 : wall / 2 )) "$@" "${resume[@]}"
