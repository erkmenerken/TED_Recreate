#!/bin/bash
#SBATCH -J ted_run
#SBATCH -A kutem
#SBATCH --qos=kutem
#SBATCH -p kutem_gpu
#SBATCH --gres=gpu:tesla_a100:1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH -o ted_run-%j.out
# Usage: sbatch scripts/run_ted.sbatch {chop|label|cath} proteins.fasta out.json [workdir]
#   chop : sequence -> ESMFold model -> Merizo/Chainsaw/UniDoc -> TED consensus domains
#   label: chop, then CATH label for every domain (exact TED match, structure, MMseqs2 lookup)
#   cath : the FASTA already holds domain sequences; label them
set -euo pipefail
R=/scratch/erkmenerken22/TED_Recreate
export HF_HOME=$R/models/hf TMPDIR=${TMPDIR:-$R/runs/tmp}
mkdir -p $TMPDIR
CMD=$1; FASTA=$(readlink -f $2); OUT=$(readlink -f -m $3); WORK=${4:-}
cd $R
$R/.venv/bin/python -m ted_recreate $CMD --fasta $FASTA --out $OUT ${WORK:+--workdir $WORK} --threads ${SLURM_CPUS_PER_TASK:-8}
