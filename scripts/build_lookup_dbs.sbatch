#!/bin/bash
#SBATCH -J ted_build_dbs
#SBATCH -A kutem
#SBATCH --qos=kutem
#SBATCH -p kutem_gpu
#SBATCH -c 24
#SBATCH --mem=220G
#SBATCH -t 2-00:00:00
#SBATCH -o /scratch/erkmenerken22/TED_Recreate/logs/build_lookup_dbs-%j.out
# Builds the CATH-label lookups (CPU only, no GPU requested):
#   db/ted_md5/            exact md5 -> TED label table (365M domains)
#   db/mmseqs/ted100       MMseqs2 DB of the 120.7M TED-100 cluster representatives with label headers
#   db/mmseqs/ted100.idx   precomputed k-mer index for the CPU server (mmseqs createindex)
#   db/mmseqs/ted100_pad   padded DB for the GPU server (mmseqs makepaddedseqdb)
# Needs data/ted/ted_365m.domain_summary... , data/ted/ted_324m_seq_clustering.cathlabels.tsv.gz and
# data/foldseek_teddb/teddb (plain sequence records in summary row order; no index/lookup needed).
set -euo pipefail
R=/scratch/erkmenerken22/TED_Recreate
cd $R
PY=$R/.venv/bin/python
export TMPDIR=$R/db/tmp SORT_MEM=60G
mkdir -p $TMPDIR
for f in data/foldseek_teddb/teddb data/ted/ted_324m_seq_clustering.cathlabels.tsv.gz; do
  test -s $f || { echo "missing $f"; exit 1; }
done
gzip -t data/ted/ted_365m.domain_summary.cath.globularity.taxid.tsv.gz || { echo "summary download incomplete"; exit 1; }

echo "=== md5 index  $(date)"
test -s db/ted_md5/labels.json || $PY scripts/build_ted_md5_index.py
echo "=== mmseqs DB  $(date)"
$PY scripts/build_ted_mmseqs_db.py
echo "=== stats  $(date)"
awk '{n++; s+=$3-2} END {printf "representatives %d, residues %.2f G, mean length %.1f\n", n, s/1e9, s/n}' db/mmseqs/ted100.index
echo "=== GPU padded DB  $(date)"
test -s db/mmseqs/ted100_pad.dbtype || tools/bin/mmseqs-gpu makepaddedseqdb db/mmseqs/ted100 db/mmseqs/ted100_pad --threads 24
echo "=== CPU index  $(date)"
test -s db/mmseqs/ted100.idx.dbtype || tools/bin/mmseqs createindex db/mmseqs/ted100 $TMPDIR --threads 24 --split-memory-limit 150G
ls -la db/mmseqs db/ted_md5
echo "=== done  $(date)"
