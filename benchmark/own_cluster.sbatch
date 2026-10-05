#!/bin/bash
#SBATCH -J ted_owncluster
#SBATCH -A cosbi
#SBATCH --qos=cosbi
#SBATCH -p cosbi
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -o own_cluster-%j.out
cd /scratch/erkmenerken22/TED_Recreate/benchmark
# rows of TED's cluster table (rep, member, label, method) whose member is one of the benchmark domains
LC_ALL=C awk -F'\t' 'NR==FNR {want[$1]=1; next} ($2 in want)' ted_domain_ids.txt ../db/mmseqs/clusters.sorted.tsv > own_cluster.tsv
wc -l own_cluster.tsv
