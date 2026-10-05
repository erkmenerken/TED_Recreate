#!/bin/bash
#SBATCH -J ted_bench_sample
#SBATCH -A cosbi
#SBATCH --qos=cosbi
#SBATCH -p cosbi
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -o 01_sample-%j.out
cd /scratch/erkmenerken22/TED_Recreate/benchmark
/scratch/erkmenerken22/TED_Recreate/.venv/bin/python 01_sample.py 150000
