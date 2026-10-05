#!/bin/bash
#SBATCH -J ted_benchmark
#SBATCH -A kutem
#SBATCH --qos=kutem
#SBATCH -p kutem_gpu
#SBATCH --gres=gpu:tesla_a100:1
#SBATCH -c 8
#SBATCH --mem=96G
#SBATCH -t 10:00:00
#SBATCH -o 03_run-%j.out
R=/scratch/erkmenerken22/TED_Recreate
export HF_HOME=$R/models/hf TMPDIR=$R/benchmark/tmp PYTHONUNBUFFERED=1; mkdir -p $TMPDIR
cd $R/benchmark
$R/.venv/bin/python 03_run.py; echo "exit code $?"
