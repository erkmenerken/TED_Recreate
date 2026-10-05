#!/bin/bash
#SBATCH -J ted_mmseqs_server
#SBATCH -A kutem
#SBATCH --qos=kutem
#SBATCH -p kutem_gpu
#SBATCH -c 16
#SBATCH --mem=200G
#SBATCH -t 1-00:00:00
#SBATCH -o /scratch/erkmenerken22/TED_Recreate/logs/mmseqs_server-%j.out
# Keeps the TED lookup DB hot and answers searches over HTTP (see ted_recreate/server.py).
#   CPU mode (default): the k-mer index lives in the page cache (~RAM of the index size), no GPU needed.
#   GPU mode: sbatch --gres=gpu:tesla_a100:1 --mem=64G scripts/mmseqs_server.sbatch gpu
# The address is written to db/mmseqs/SERVER; ted_recreate.classify uses it automatically.
# Stop with: scancel <jobid>
MODE=${1:-cpu}
R=/scratch/erkmenerken22/TED_Recreate
cd $R
export TMPDIR=$R/db/tmp_server
mkdir -p $TMPDIR
exec $R/.venv/bin/python -m ted_recreate.server --mode $MODE --threads ${SLURM_CPUS_PER_TASK:-16}
