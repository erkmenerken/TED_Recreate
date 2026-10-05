"""Paths shared by the TED re-creation code. Everything lives under TED_Recreate/."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Official TED segmentation code (psipred/ted-tools, ted_consensus_1.0) with UniDoc v20250514 installed
TED_CONSENSUS = ROOT / "vendor" / "ted-tools" / "ted_consensus_1.0"
RUN_SEGMENTATION = TED_CONSENSUS / "run_segmentation.sh"
MERIZO_SEARCH = ROOT / "vendor" / "merizo_search" / "merizo_search"   # Foldclass network + TM-align binary

VENV = ROOT / ".venv"
PYTHON = VENV / "bin" / "python"
TOOLS_BIN = ROOT / "tools" / "bin"
MMSEQS = TOOLS_BIN / "mmseqs"
MMSEQS_GPU = TOOLS_BIN / "mmseqs-gpu"      # MMseqs2-GPU build (makepaddedseqdb, gpuserver, --gpu 1)
FOLDSEEK = TOOLS_BIN / "foldseek"
FOLDSEEK_TED = TOOLS_BIN / "foldseek8"     # release 8-ef4e960, the closest release to TED's commit a435618
TMALIGN = TOOLS_BIN / "TMalign"

HF_HOME = ROOT / "models" / "hf"
ESMFOLD_MODEL = "facebook/esmfold_v1"

DATA = ROOT / "data"
DB = ROOT / "db"

# CATH 4.3 reference (scripts/build_cath_reference.py); S40 stands in for TED's unpublished SSG5 list
CATH_LABELS = DB / "cath" / "cath43_labels.tsv"
CATH_NAMES = DB / "cath" / "cath43_names.tsv"
CATH_FOLDSEEK_DB = DATA / "cath" / "foldseek_s40" / "cath43_s40"
FOLDCLASS_DB = DB / "cath" / "foldclass_s40"

# TED lookups
TED_MD5_INDEX = DB / "ted_md5"                    # scripts/build_ted_md5_index.py (365M domains)
TED_MMSEQS_DB = DB / "mmseqs" / "ted100"          # scripts/build_ted_mmseqs_db.sh (cluster-search DB)
