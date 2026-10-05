"""Getting a 3D model for a protein so it can be chopped the TED way.

TED chopped AlphaFold DB v4 models. For a UniProt accession we therefore fetch the AFDB model
(AFDB now serves v6 files; for models created in the v4 era the coordinates are identical).
For an arbitrary sequence there is no AFDB model, so we predict one with ESMFold. That is the one
step that cannot be identical to TED: ESMFold is not AlphaFold2, and its models (and hence the
domain boundaries) can differ, especially for low-confidence regions.
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import urllib.request
from pathlib import Path

from . import config

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
VALID_AA = set(THREE_TO_ONE.values())


def md5(seq: str) -> str:
    return hashlib.md5(seq.encode("utf-8")).hexdigest()


def clean_sequence(seq: str) -> str:
    seq = "".join(seq.split()).upper()
    bad = set(seq) - VALID_AA
    if bad:
        raise ValueError(f"sequence contains non-standard residues {sorted(bad)}; TED models only cover the 20 standard amino acids")
    return seq


def read_pdb_ca(pdb_path: str | Path, chain: str | None = None):
    """Return (sequence, resnums, plddt/bfactor per residue) from the CA atoms of one chain."""
    seq, resi, bfac = [], [], []
    seen = set()
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                break
            if not line.startswith("ATOM") or line[12:16].strip() != "CA":
                continue
            ch = line[21]
            if chain is None:
                chain = ch
            if ch != chain:
                continue
            key = (line[22:27])
            if key in seen:
                continue
            seen.add(key)
            seq.append(THREE_TO_ONE.get(line[17:20], "X"))
            resi.append(int(line[22:26]))
            bfac.append(float(line[60:66]))
    return "".join(seq), resi, bfac


def _atom_records_only(path: Path):
    """Keep ATOM records, one TER after them, then END, every record padded to 80 columns (the layout of AFDB
    files). The 2022 UniDoc binary TED used crashes on any line shorter than ~22 characters that it reads, e.g.
    ESMFold's 'PARENT N/A' header or a bare 'END' left after pdb_selres drops the TER record."""
    atoms = [l for l in path.read_text().splitlines() if l.startswith("ATOM")]
    last = atoms[-1]
    ter = f"TER   {int(last[6:11]) + 1:5d}      {last[17:20]} {last[21]}{last[22:26]}"
    path.write_text("\n".join(l.ljust(80) for l in atoms + [ter, "END"]) + "\n")


def normalise_pdb(src: str | Path, dst: str | Path) -> Path:
    """Make a single-chain PDB that the TED tools expect: first model, first chain renamed to A, no HETATM,
    residues renumbered from 1, only ATOM/TER/END records. For AFDB models only the header lines are dropped
    (the parsers read ATOM records only)."""
    src, dst = Path(src), Path(dst)
    seq, resi, _ = read_pdb_ca(src)
    chains = set()
    with open(src) as fh:
        for line in fh:
            if line.startswith("ATOM"):
                chains.add(line[21])
    if chains == {"A"} and resi == list(range(1, len(resi) + 1)):
        dst.write_bytes(src.read_bytes())
    else:
        first_chain = None
        with open(src) as fh:
            for line in fh:
                if line.startswith("ATOM"):
                    first_chain = line[21]
                    break
        bin_dir = config.VENV / "bin"
        cmd = (f"{bin_dir}/pdb_selmodel -1 '{src}' | {bin_dir}/pdb_selchain -{first_chain} | {bin_dir}/pdb_delhetatm | "
               f"{bin_dir}/pdb_rplchain -{first_chain}:A | {bin_dir}/pdb_reres -1 | {bin_dir}/pdb_tidy > '{dst}'")
        subprocess.run(cmd, shell=True, check=True)
    _atom_records_only(dst)
    return dst


def fetch_afdb_model(uniprot: str, outdir: str | Path) -> tuple[Path, dict]:
    """Download the AlphaFold DB model for a UniProt accession. Returns (path, api_record)."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot}",
                                 headers={"User-Agent": "ted-recreate"})
    with urllib.request.urlopen(req, timeout=60) as r:
        records = json.load(r)
    rec = next((r for r in records if r["entryId"].endswith("-F1")), records[0])
    # TED used v4 file names; keep that naming so outputs line up with TED ids
    name = f"{rec['entryId']}-model_v4"
    path = outdir / f"{name}.pdb"
    req = urllib.request.Request(rec["pdbUrl"], headers={"User-Agent": "ted-recreate"})
    with urllib.request.urlopen(req, timeout=120) as r:
        path.write_bytes(r.read())
    return path, rec


class ESMFold:
    """Thin wrapper around the HuggingFace ESMFold port. Loads once, folds many sequences."""

    def __init__(self, device: str | None = None, chunk_size: int | None = 64):
        os.environ.setdefault("HF_HOME", str(config.HF_HOME))
        import torch
        from transformers import AutoTokenizer, EsmForProteinFolding

        self.torch = torch
        self.device = device or ("cuda" if torch.cuda.is_available() else "cpu")
        self.tokenizer = AutoTokenizer.from_pretrained(config.ESMFOLD_MODEL)
        self.model = EsmForProteinFolding.from_pretrained(config.ESMFOLD_MODEL, low_cpu_mem_usage=True)
        self.model.eval()
        if self.device.startswith("cuda"):
            self.model.esm = self.model.esm.half()   # standard trick: ESM-2 trunk in fp16, folding head in fp32
        self.model = self.model.to(self.device)
        if chunk_size:
            self.model.trunk.set_chunk_size(chunk_size)

    def fold(self, seq: str, out_pdb: str | Path) -> dict:
        torch = self.torch
        with torch.no_grad():
            tok = self.tokenizer([seq], return_tensors="pt", add_special_tokens=False)["input_ids"].to(self.device)
            out = self.model(tok)
        pdb_str = self.model.output_to_pdb(out)[0]
        plddt = out["plddt"][0, :, 1].float().cpu().numpy()  # per-residue, CA atom
        scale = 100.0 if plddt.max() <= 1.0 else 1.0
        # rewrite B-factors so every atom carries its residue pLDDT on the 0-100 scale (as in AFDB files)
        lines, res_plddt = [], {}
        for line in pdb_str.splitlines():
            if line.startswith("ATOM"):
                r = int(line[22:26])
                val = float(plddt[r - 1]) * scale
                line = f"{line[:60]}{val:6.2f}{line[66:]}"
            lines.append(line)
        Path(out_pdb).write_text("\n".join(lines) + "\n")
        return {"mean_plddt": float(plddt.mean() * scale), "ptm": float(out["ptm"].item())}
