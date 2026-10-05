"""TED-style domain chopping: structure -> Merizo + Chainsaw + UniDoc -> consensus domains.

The segmentation itself is TED's own code, run unmodified (vendor/ted-tools/ted_consensus_1.0/
run_segmentation.sh), so for a given structure the output is what TED's released pipeline gives:
  1. Merizo on the full model (dual pass, iterative re-segmentation, UniDoc-style merging)
  2. UniDoc on the residues Merizo put in domains (Merizo strips the non-domain residues first)
  3. Chainsaw on the full model (+1 residue offset, Chainsaw is 0-based)
  each followed by filter_domains.py (drop segments < 5 residues, domains < 25 residues)
  4. consensus: pairwise IoU between all predicted domains, edges where IoU >= 0.7, connected
     components of size 3 = high, 2 = medium, 1 = low; final range = residues shared by all
     members of the component; then the same 5/25 filter.
Only high + medium domains are TED domains.

Note on numbering: TED's consensus code adds 1 to every residue index before intersecting
(domain_consensus.py: domstr_to_assignment(..., offset=1)), so consensus ranges sit one residue
after the parsers' ranges. The published TED database has the same shift, so `chopping` below is
reported in TED's convention; `chopping_parser_frame` undoes the shift.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field, asdict
from pathlib import Path

from . import config
from .structure import (ESMFold, clean_sequence, fetch_afdb_model, md5, normalise_pdb,
                        read_pdb_ca)

LEVELS = ("high", "medium", "low")


def parse_chopping(chopping: str) -> list[list[tuple[int, int]]]:
    """'11-41_290-389,54-288' -> [[(11, 41), (290, 389)], [(54, 288)]]; 'na'/'NULL'/'0' -> []"""
    if chopping in ("na", "NULL", "NO_SS", "0", "", None):
        return []
    doms = []
    for dom in chopping.split(","):
        segs = []
        for seg in dom.split("_"):
            a, b = seg.split("-")
            segs.append((int(a), int(b)))
        doms.append(segs)
    return doms


def segments_to_str(segs) -> str:
    return "_".join(f"{a}-{b}" for a, b in segs)


def extract(seq: str, segs) -> str:
    return "".join(seq[a - 1:b] for a, b in segs)


@dataclass
class Domain:
    ted_id: str
    consensus_level: str
    chopping: str                  # TED convention (same frame as the TED database)
    chopping_parser_frame: str     # shifted back by one residue = what the parsers actually agreed on
    nres: int
    num_segments: int
    sequence: str
    md5_domain: str
    plddt: float | None


@dataclass
class ChopResult:
    name: str
    sequence: str
    md5: str
    nres: int
    structure: str
    structure_source: str
    methods: dict = field(default_factory=dict)      # merizo/chainsaw/unidoc -> chopping string
    consensus: dict = field(default_factory=dict)    # high/medium/low -> chopping string
    domains: list = field(default_factory=list)      # high + medium domains, named TED01.. from N-terminus
    low_consensus: list = field(default_factory=list)
    warnings: list = field(default_factory=list)

    def residue_labels(self) -> list[int]:
        """Per-residue domain index (1..n, TED numbering order) with 0 for residues not in a TED domain."""
        lab = [0] * self.nres
        for i, d in enumerate(self.domains, 1):
            for a, b in parse_chopping(d.chopping)[0]:
                for r in range(a, b + 1):
                    lab[r - 1] = i
        return lab

    def to_dict(self) -> dict:
        d = asdict(self)
        d["residue_labels"] = self.residue_labels()
        return d


def _read_table(path: Path) -> dict:
    out = {}
    if path.exists():
        for line in path.read_text().splitlines():
            f = line.split("\t")
            if f and f[0] and not f[0].startswith("Done"):
                out[f[0]] = f
    return out


def run_ted_segmentation(pdbs: dict[str, Path], workdir: Path) -> dict[str, dict]:
    """Run TED's run_segmentation.sh on {name: pdb} and return the parsed per-target outputs."""
    workdir = Path(workdir)
    inp, out = workdir / "input", workdir / "out"
    inp.mkdir(parents=True, exist_ok=True)
    for name, pdb in pdbs.items():
        normalise_pdb(pdb, inp / f"{name}.pdb")
    venv_link = workdir / "ted_consensus"          # run_segmentation.sh activates ./ted_consensus
    if not venv_link.exists():
        venv_link.symlink_to(config.VENV)
    env = dict(os.environ)
    env["PATH"] = f"{config.VENV / 'bin'}:{env.get('PATH', '')}"
    proc = subprocess.run(["bash", str(config.RUN_SEGMENTATION), "-i", "input", "-o", "out"],
                          cwd=workdir, env=env, capture_output=True, text=True)
    if proc.returncode != 0:
        logs = "\n".join(f"--- {p.name}\n{p.read_text()[-2000:]}" for p in sorted(out.glob("*.log")))
        raise RuntimeError(f"run_segmentation.sh failed ({proc.returncode})\n{proc.stdout}\n{proc.stderr}\n{logs}")
    tables = {m: _read_table(out / f"chopping_{m}.txt") for m in ("merizo", "chainsaw", "unidoc")}
    cons = _read_table(out / "consensus.tsv")
    results = {}
    for name in pdbs:
        c = cons[name]
        results[name] = {
            "md5": c[1], "nres": int(c[2]),
            "methods": {m: tables[m][name][4] if name in tables[m] else "NULL" for m in tables},
            "consensus": {"high": c[6], "medium": c[7], "low": c[8]},
        }
    return results


def _assemble(name, seq, pdb_path, source, seg) -> ChopResult:
    _, _, plddt = read_pdb_ca(pdb_path)
    res = ChopResult(name=name, sequence=seq, md5=md5(seq), nres=len(seq), structure=str(pdb_path),
                     structure_source=source, methods=seg["methods"], consensus=seg["consensus"])
    if seg["md5"] != res.md5:
        raise RuntimeError(f"{name}: sequence in structure does not match input sequence")
    doms = []
    for level in ("high", "medium"):
        for segs in parse_chopping(seg["consensus"][level]):
            doms.append((segs[0][0], level, segs))
    doms.sort()
    for i, (_, level, segs) in enumerate(doms, 1):
        dseq = extract(seq, segs)
        idx = [r - 1 for a, b in segs for r in range(a, b + 1) if r <= len(seq)]
        res.domains.append(Domain(
            ted_id=f"{name}_TED{i:02d}", consensus_level=level, chopping=segments_to_str(segs),
            chopping_parser_frame=segments_to_str([(a - 1, b - 1) for a, b in segs]),
            nres=sum(b - a + 1 for a, b in segs), num_segments=len(segs), sequence=dseq,
            md5_domain=md5(dseq),
            plddt=round(sum(plddt[j] for j in idx) / len(idx), 4) if idx else None))
    res.low_consensus = [segments_to_str(s) for s in parse_chopping(seg["consensus"]["low"])]
    for m, c in seg["methods"].items():
        if c == "NO_SS":
            res.warnings.append(f"{m} failed on this structure (NO_SS); consensus used the other parsers only")
    return res


def _check_name(name: str):
    if not name or any(c in name for c in "/ \t,"):
        raise ValueError(f"bad name {name!r}: no slashes, spaces, tabs or commas")


def ted_chop_batch(items, workdir: str | Path | None = None, keep_workdir: bool = False,
                   device: str | None = None) -> list[ChopResult]:
    """Chop many proteins in one go (models are loaded once).

    items: list of dicts, each with one of
        {"sequence": "MKT...", "name": "optional_id"}   -> folded with ESMFold
        {"uniprot": "P69905"}                            -> AFDB model (what TED used)
        {"pdb": "path/to/model.pdb", "name": "..."}      -> your own model (first chain used)
    """
    tmp = Path(workdir) if workdir else Path(tempfile.mkdtemp(prefix="ted_chop_", dir=os.environ.get("TMPDIR")))
    tmp.mkdir(parents=True, exist_ok=True)
    struct_dir = tmp / "structures"
    struct_dir.mkdir(exist_ok=True)
    prepared = {}   # name -> (seq, pdb, source)
    esm = None
    try:
        for k, it in enumerate(items):
            if isinstance(it, str):
                it = {"sequence": it}
            if "pdb" in it:
                pdb = Path(it["pdb"])
                name = it.get("name") or pdb.stem
                _check_name(name)
                seq = read_pdb_ca(pdb)[0]
                prepared[name] = (seq, pdb, "user_pdb")
            elif "uniprot" in it:
                pdb, rec = fetch_afdb_model(it["uniprot"], struct_dir)
                name = it.get("name") or pdb.stem
                prepared[name] = (read_pdb_ca(pdb)[0], pdb, f"AFDB v{rec.get('latestVersion')} ({rec['entryId']})")
            else:
                seq = clean_sequence(it["sequence"])
                name = it.get("name") or f"query_{md5(seq)[:10]}"
                _check_name(name)
                if esm is None:
                    esm = ESMFold(device=device)
                pdb = struct_dir / f"{name}.pdb"
                info = esm.fold(seq, pdb)
                prepared[name] = (seq, pdb, f"ESMFold (mean pLDDT {info['mean_plddt']:.1f}, pTM {info['ptm']:.2f})")
        if esm is not None:
            del esm
            import torch
            torch.cuda.empty_cache()
        seg = run_ted_segmentation({n: p for n, (_, p, _) in prepared.items()}, tmp / "segmentation")
        results = []
        for name, (seq, pdb, source) in prepared.items():
            r = _assemble(name, seq, pdb, source, seg[name])
            if workdir is None and not keep_workdir and str(pdb).startswith(str(tmp)):
                r.structure = ""          # temporary model is deleted; pass workdir= to keep it
            results.append(r)
        return results
    finally:
        if not keep_workdir and workdir is None:
            shutil.rmtree(tmp, ignore_errors=True)


def ted_chop(sequence: str | None = None, *, uniprot: str | None = None, pdb: str | None = None,
             name: str | None = None, **kw) -> ChopResult:
    """Chop one protein into TED consensus domains. Give exactly one of sequence / uniprot / pdb."""
    given = [x is not None for x in (sequence, uniprot, pdb)]
    if sum(given) != 1:
        raise ValueError("give exactly one of sequence=, uniprot=, pdb=")
    item = {"sequence": sequence} if sequence else {"uniprot": uniprot} if uniprot else {"pdb": pdb}
    if name:
        item["name"] = name
    return ted_chop_batch([item], **kw)[0]
