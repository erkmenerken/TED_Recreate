"""CATH labels for domains, following how TED assigned them.

TED labelled its 324M TED-100 domains in three ways (Lau et al. 2024, Supp. Methods):
  1. Foldseek structure search against CATH 4.3 SSG5 representatives
       search:  -e 0.108662 -c 0.366757 --cov-mode 5 -s 10
       H-level: E < 0.019, TM-score > 0.56, coverage >= 0.367
       T-level: E < 0.108662, TM-score >= 0.42, coverage >= 0.786
  2. For domains Foldseek could not label: Foldclass embedding nearest neighbour (Merizo-search) in CATH,
     confirmed by TM-align > 0.5 normalised by the TED domain length -> T-level label
  3. MMseqs2 sequence clusters (>= 50% identity, >= 90% coverage): unlabelled domains in a cluster with
     labelled members can take the members' label ("annotations transferrable", 11.6M domains)

Here the same tiers are exposed for a new domain:
  exact     md5 of the domain sequence found in TED -> TED's published label (instant, needs db/ted_md5)
  structure Foldseek (TED thresholds) then Foldclass+TM-align, against CATH 4.3 S40 (public stand-in for
            the unpublished SSG5 list; same size, so E-values are on the same scale)
  sequence  MMseqs2 search against all TED-100 domain sequences (cluster search over TED's own 50% clusters);
            a hit that would have fallen into the same TED cluster (>= 50% id, >= 90% coverage) transfers its
            label; weaker homologs are reported as hints only
"""
from __future__ import annotations

import hashlib
import json
import os
import pickle
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field, asdict
from pathlib import Path

import numpy as np

from . import config

# --- TED's thresholds (Supp. Methods, "Domain classification using Foldseek") -----------------------------
FS_SEARCH_EVALUE = 0.108662
FS_SEARCH_COV = 0.366757
FS_SEARCH_COVMODE = 5
FS_SENSITIVITY = 10
H_EVALUE, H_TM, H_COV = 0.019, 0.56, 0.367
T_EVALUE, T_TM, T_COV = 0.108662, 0.42, 0.786
FOLDCLASS_TM = 0.5                       # TM-align normalised by the query (TED) domain
SEQ_TRANSFER_ID, SEQ_TRANSFER_COV = 0.5, 0.9   # TED's MMseqs2 clustering criterion

# Which Foldseek columns implement TED's "TM-score" and "coverage"? The paper does not say. Calibrated with
# scripts/calibrate_foldseek_rules.py on 387 TED-100 domains (runs/calibrate): TM = max(qtmscore, ttmscore) and
# coverage = length ratio (cov-mode 5, as in TED's search) reproduce TED's Foldseek outcome for 375/387 domains
# (H 205/207, T 15/23, 2/157 spurious) with foldseek 8; alntmscore 368/387, qtmscore 366, ttmscore 351.
FS_TM_COLUMN = os.environ.get("TED_FS_TM_COLUMN", "maxtmscore")
FS_COV_COLUMN = os.environ.get("TED_FS_COV_COLUMN", "lenratio")   # min(qlen,tlen)/max(qlen,tlen), cov-mode 5

FS_FORMAT = "query,target,evalue,bits,prob,alntmscore,qtmscore,ttmscore,qcov,tcov,fident,alnlen,qlen,tlen,lddt"


def md5(seq: str) -> str:
    return hashlib.md5(seq.encode("utf-8")).hexdigest()


def _cath_labels() -> dict:
    lab = {}
    with open(config.CATH_LABELS) as fh:
        next(fh)
        for line in fh:
            f = line.split("\t")
            lab[f[0]] = f[1]
    return lab


def _cath_names() -> dict:
    names = {}
    with open(config.CATH_NAMES) as fh:
        next(fh)
        for line in fh:
            node, _, name = line.rstrip("\n").split("\t")
            names[node] = name
    return names


def truncate(label: str, level: str) -> str:
    """'3.40.50.300', 'T' -> '3.40.50'"""
    parts = label.split(".")
    return ".".join(parts[:4 if level == "H" else 3])


# ---------------------------------------------------------------------------------------------------------
# Tier 1: exact sequence match against TED (md5_domain)
# ---------------------------------------------------------------------------------------------------------
class TedMd5Index:
    """Sorted md5 table built by scripts/build_ted_md5_index.py from the 365M-domain TED summary.

    Files (numpy, memory-mapped, so lookups do not load the table): md5_hi.npy, md5_lo.npy (uint64, sorted),
    label.npy (int32 code into labels.json), level.npy (uint8: 0 none, 1 T, 2 H), method.npy (uint8),
    ted_row.npy (int64 row into ted_ids.txt offsets) + ted_ids.offsets.npy / ted_ids.txt.
    """

    LEVELS = {0: "-", 1: "T", 2: "H"}

    def __init__(self, path: Path = None):
        path = Path(path or config.TED_MD5_INDEX)
        if not (path / "md5_hi.npy").exists():
            raise FileNotFoundError(f"TED md5 index not built yet ({path}); run scripts/build_ted_md5_index.py")
        self.path = path
        self.hi = np.load(path / "md5_hi.npy", mmap_mode="r")
        self.lo = np.load(path / "md5_lo.npy", mmap_mode="r")
        self.label = np.load(path / "label.npy", mmap_mode="r")
        self.level = np.load(path / "level.npy", mmap_mode="r")
        self.method = np.load(path / "method.npy", mmap_mode="r")
        self.ted_off = np.load(path / "ted_ids.offsets.npy", mmap_mode="r")
        self.ted_row = np.load(path / "ted_row.npy", mmap_mode="r")
        meta = json.load(open(path / "labels.json"))
        self.labels, self.methods = meta["labels"], meta["methods"]
        self._ids = open(path / "ted_ids.txt", "rb")

    def _ted_id(self, row: int) -> str:
        a, b = int(self.ted_off[row]), int(self.ted_off[row + 1])
        self._ids.seek(a)
        return self._ids.read(b - a).decode().rstrip("\n")

    def lookup(self, seq: str, max_ids: int = 20) -> list[dict]:
        h = md5(seq)
        hi, lo = np.uint64(int(h[:16], 16)), np.uint64(int(h[16:], 16))
        a = int(np.searchsorted(self.hi, hi, side="left"))
        b = int(np.searchsorted(self.hi, hi, side="right"))
        out = []
        for i in range(a, b):
            if self.lo[i] != lo:
                continue
            out.append({"ted_id": self._ted_id(int(self.ted_row[i])) if len(out) < max_ids else None,
                        "cath_label": self.labels[int(self.label[i])],
                        "level": self.LEVELS[int(self.level[i])],
                        "method": self.methods[int(self.method[i])]})
        return out


# ---------------------------------------------------------------------------------------------------------
# Tier 2: structure (Foldseek with TED's thresholds, then Foldclass + TM-align)
# ---------------------------------------------------------------------------------------------------------
def foldseek_search(pdb_dir: Path, workdir: Path, threads: int = 8, foldseek: Path = None) -> list[dict]:
    foldseek = foldseek or config.FOLDSEEK_TED
    out = workdir / "foldseek.m8"
    cmd = [str(foldseek), "easy-search", str(pdb_dir), str(config.CATH_FOLDSEEK_DB), str(out), str(workdir / "fs_tmp"),
           "-e", str(FS_SEARCH_EVALUE), "-c", str(FS_SEARCH_COV), "--cov-mode", str(FS_SEARCH_COVMODE),
           "-s", str(FS_SENSITIVITY), "--max-seqs", "1000", "--threads", str(threads),
           "--format-output", FS_FORMAT]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"foldseek failed:\n{p.stderr[-3000:]}")
    cols = FS_FORMAT.split(",")
    hits = []
    for line in out.read_text().splitlines():
        h = dict(zip(cols, line.split("\t")))
        for k in cols[2:]:
            h[k] = float(h[k])
        h["query"] = h["query"].removesuffix(".pdb")
        h["target"] = h["target"].removesuffix(".pdb")
        h["lenratio"] = min(h["qlen"], h["tlen"]) / max(h["qlen"], h["tlen"])
        h["maxtmscore"] = max(h["qtmscore"], h["ttmscore"])
        h["mincov"] = min(h["qcov"], h["tcov"])
        hits.append(h)
    return hits


def ted_foldseek_rule(hits: list[dict], cath_labels: dict, tm_col: str = None, cov_col: str = None) -> dict | None:
    """Apply TED's H then T thresholds to a query's Foldseek hits. Returns the assignment or None."""
    tm_col, cov_col = tm_col or FS_TM_COLUMN, cov_col or FS_COV_COLUMN
    best = {}
    for h in sorted(hits, key=lambda x: (x["evalue"], -x["bits"])):
        if h["evalue"] < H_EVALUE and h[tm_col] > H_TM and h[cov_col] >= H_COV:
            best.setdefault("H", h)
        if h["evalue"] < T_EVALUE and h[tm_col] >= T_TM and h[cov_col] >= T_COV:
            best.setdefault("T", h)
    for level in ("H", "T"):
        if level in best:
            h = best[level]
            return {"cath_label": truncate(cath_labels[h["target"]], level), "level": level, "method": "foldseek",
                    "evidence": {"cath_domain": h["target"], "evalue": h["evalue"], "tm": h[tm_col],
                                 "coverage": h[cov_col], "prob": h["prob"], "fident": h["fident"]}}
    return None


class Foldclass:
    """Merizo-search's Foldclass (EGNN embedding) nearest neighbour in CATH S40, confirmed with TM-align."""

    def __init__(self, device: str = "cpu"):
        import sys
        import torch
        sys.path.insert(0, str(config.MERIZO_SEARCH))
        from programs.Foldclass.nndef_fold_egnn_embed import FoldClassNet
        from programs.Foldclass import utils as fcutils
        self.torch, self.fcutils = torch, fcutils
        self.device = torch.device(device)
        net = FoldClassNet(128).eval().to(self.device)
        weights = config.MERIZO_SEARCH / "programs" / "Foldclass" / "FINAL_foldclass_model.pt"
        net.load_state_dict(torch.load(weights, map_location="cpu"), strict=False)
        self.net = net
        self.db = torch.load(config.FOLDCLASS_DB.with_suffix(".pt"), map_location="cpu").to(self.device)
        self.index = pickle.load(open(config.FOLDCLASS_DB.with_suffix(".index"), "rb"))
        self.tmalign = config.MERIZO_SEARCH / "programs" / "Foldclass" / "tmalign"

    def search(self, pdb: Path, cath_labels: dict, workdir: Path) -> dict:
        torch, fcu = self.torch, self.fcutils
        q = fcu.read_pdb(str(pdb), pdb_chain="A")
        with torch.no_grad():
            emb = self.net(torch.from_numpy(q["coords"]).unsqueeze(0).to(self.device))
            cos = torch.nn.functional.cosine_similarity(self.db, emb, dim=-1)
            score, idx = torch.max(cos, dim=0)
        name, coords, seq = self.index[int(idx)]
        dom = Path(name).stem
        tfile = fcu.write_pdb(str(workdir), coords, seq)
        qfile = fcu.write_pdb(str(workdir), q["coords"], q["seq"])
        p = subprocess.run([str(self.tmalign), qfile, tfile], capture_output=True, text=True)
        tm = fcu.extract_tmalign_values(p.stdout)
        for f in (qfile, tfile):
            os.remove(f)
        res = {"cath_domain": dom, "cosine": float(score), "q_tm": tm["qtm"], "t_tm": tm["ttm"], "rmsd": tm["rmsd"],
               "len_ali": tm["len_ali"]}
        if tm["qtm"] > FOLDCLASS_TM:
            return {"cath_label": truncate(cath_labels[dom], "T"), "level": "T", "method": "foldclass", "evidence": res}
        return {"cath_label": None, "level": "-", "method": "foldclass", "evidence": res}


# ---------------------------------------------------------------------------------------------------------
# Tier 3: MMseqs2 against TED's sequence-cluster representatives (the "MMseqs server")
# ---------------------------------------------------------------------------------------------------------
MM_FORMAT = "query,target,theader,pident,fident,alnlen,qstart,qend,qlen,tstart,tend,tlen,qcov,tcov,evalue,bits"


def _parse_rep_header(header: str) -> dict:
    """'REP TOP LEVEL n=.. nlab=.. share=.. own=LABEL:METHOD' (scripts/build_ted_mmseqs_db.py)"""
    parts = header.split()
    d = {"ted_rep": parts[0], "cluster_label": parts[1], "cluster_level": parts[2]}
    for kv in parts[3:]:
        k, v = kv.split("=", 1)
        d[k] = v
    d["n_members"], d["n_labelled"], d["label_share"] = int(d.pop("n")), int(d.pop("nlab")), float(d.pop("share"))
    d["rep_label"], d["rep_method"] = d.pop("own").split(":")
    return d


def mmseqs_search(seqs: dict, workdir: Path, threads: int = 8, sensitivity: float = 7.5, max_seqs: int = 300,
                  db_load_mode: int | None = None, db: Path = None, gpu: bool = False) -> dict[str, list[dict]]:
    """Search query sequences against the 120.7M TED-100 cluster representatives.

    db_load_mode 2 (mmap) + a pre-built index kept in the page cache (mmseqs touchdb) is the CPU "server" mode;
    gpu=True uses MMseqs2-GPU against the padded DB (and a running `mmseqs gpuserver` if TED_MMSEQS_GPU_SERVER=1).
    """
    db = Path(db or config.TED_MMSEQS_DB)
    if gpu:
        db = Path(str(db) + "_pad")
    if not Path(str(db) + ".dbtype").exists():
        raise FileNotFoundError(f"TED MMseqs2 DB not built yet ({db}); run scripts/build_lookup_dbs.sbatch")
    mm = str(config.MMSEQS_GPU if gpu else config.MMSEQS)
    qfa = workdir / "query.fasta"
    qfa.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    qdb, res, tsv, tmp = workdir / "qdb", workdir / "res", workdir / "res.tsv", workdir / "mm_tmp"
    if db_load_mode is None:
        db_load_mode = int(os.environ.get("TED_MMSEQS_DB_LOAD_MODE", "0"))

    def run(*a):
        p = subprocess.run([mm, *map(str, a)], capture_output=True, text=True)
        if p.returncode != 0:
            raise RuntimeError(f"mmseqs {a[0]} failed:\n{p.stdout[-2000:]}\n{p.stderr[-2000:]}")
    run("createdb", qfa, qdb, "-v", 1)
    extra = ["--gpu", 1, "--gpu-server", int(os.environ.get("TED_MMSEQS_GPU_SERVER", "0"))] if gpu else ["-s", sensitivity]
    run("search", qdb, db, res, tmp, *extra, "--max-seqs", max_seqs, "-e", "1e-3", "-a",
        "--db-load-mode", db_load_mode, "--threads", threads, "-v", 1)
    run("convertalis", qdb, db, res, tsv, "--format-output", MM_FORMAT, "--db-load-mode", db_load_mode,
        "--threads", threads, "-v", 1)
    cols = MM_FORMAT.split(",")
    out = {k: [] for k in seqs}
    for line in tsv.read_text().splitlines():
        h = dict(zip(cols, line.split("\t")))
        for k in cols[3:]:
            h[k] = float(h[k])
        h.update(_parse_rep_header(h.pop("theader")))
        h["lenratio"] = min(h["qlen"], h["tlen"]) / max(h["qlen"], h["tlen"])
        h["joins_ted_cluster"] = h["fident"] >= SEQ_TRANSFER_ID and h["lenratio"] >= SEQ_TRANSFER_COV
        out[h["query"]].append(h)
    return out


def _gpu_available() -> bool:
    """A GPU we may use for a one-off search: only inside a Slurm job (the login node has a small shared GPU)."""
    if not os.environ.get("SLURM_JOB_ID"):
        return False
    try:
        import torch
        return torch.cuda.is_available()
    except Exception:
        return False


def sequence_transfer(hits: list[dict]) -> dict | None:
    """TED-style transfer. TED's clusters are linclust clusters (--min-seq-id 0.5 -c 0.9 --cov-mode 5): a member
    aligns to its representative with >= 50% identity and the shorter sequence is >= 90% of the longer one.
    A query meeting that against a representative would have joined its cluster, and TED gave cluster members
    the labels of labelled cluster-mates. Among qualifying representatives the best bit score wins (as in
    linclust, a sequence joins one cluster). Returns the assignment, also when the cluster is unlabelled."""
    joins = [h for h in hits if h["joins_ted_cluster"]]
    if not joins:
        return None
    h = max(joins, key=lambda h: h["bits"])
    label = None if h["cluster_label"] == "-" else h["cluster_label"]
    return {"cath_label": label, "level": h["cluster_level"] if label else "-", "method": "mmseqs-transfer",
            "evidence": {"ted_rep": h["ted_rep"], "fident": h["fident"], "lenratio": h["lenratio"],
                         "qcov": h["qcov"], "tcov": h["tcov"], "evalue": h["evalue"], "bits": h["bits"],
                         "cluster_size": h["n_members"], "cluster_labelled": h["n_labelled"],
                         "cluster_label_share": h["label_share"], "rep_own_label": h["rep_label"],
                         "rep_label_method": h["rep_method"], "n_qualifying_clusters": len(joins)}}


# ---------------------------------------------------------------------------------------------------------
# Putting it together
# ---------------------------------------------------------------------------------------------------------
@dataclass
class LabelResult:
    name: str
    sequence: str
    cath_label: str | None          # final label ('3.40.50.300' = H level, '3.40.50' = T level) or None
    level: str                      # 'H', 'T' or '-'
    source: str                     # which tier decided: ted-exact / foldseek / foldclass / mmseqs-transfer / none
    name_of_label: str | None = None
    ted_exact: list = field(default_factory=list)
    structure: dict | None = None
    sequence_transfer: dict | None = None
    sequence_hits: list = field(default_factory=list)
    remote_hint: dict | None = None
    notes: list = field(default_factory=list)

    def to_dict(self):
        return asdict(self)


def write_domain_pdb(chain_pdb: str | Path, segments, out: Path):
    """Cut a domain out of a chain model (residue ranges as in the chopping string, TED convention)."""
    keep = set()
    for a, b in segments:
        keep.update(range(a, b + 1))
    with open(chain_pdb) as fh, open(out, "w") as fo:
        for line in fh:
            if line.startswith("ATOM") and int(line[22:26]) in keep:
                fo.write(line)
        fo.write("END\n")
    return out


class Classifier:
    """Load once, label many domains. Each tier is optional and is skipped (with a note) if its data are missing."""

    def __init__(self, tiers=("exact", "structure", "sequence"), device: str = "cpu", threads: int = 8,
                 mmseqs_db: Path = None, db_load_mode: int | None = None, use_server: bool = True):
        self.tiers, self.threads, self.use_server = tuple(tiers), threads, use_server
        self.mmseqs_db, self.db_load_mode = mmseqs_db, db_load_mode
        self.cath_labels, self.cath_names = _cath_labels(), _cath_names()
        self.notes = []
        self.md5 = None
        if "exact" in self.tiers:
            try:
                self.md5 = TedMd5Index()
            except FileNotFoundError as e:
                self.notes.append(str(e))
        self.foldclass = Foldclass(device=device) if "structure" in self.tiers else None

    def label(self, items: list[dict], workdir: str | Path | None = None) -> list[LabelResult]:
        """items: dicts with 'name', 'sequence' and optionally 'pdb' (domain structure, chain A)."""
        tmp = Path(workdir) if workdir else Path(tempfile.mkdtemp(prefix="ted_cath_", dir=os.environ.get("TMPDIR")))
        tmp.mkdir(parents=True, exist_ok=True)
        try:
            return self._label(items, tmp)
        finally:
            if workdir is None:
                shutil.rmtree(tmp, ignore_errors=True)

    def _label(self, items, tmp: Path) -> list[LabelResult]:
        results = {it["name"]: LabelResult(name=it["name"], sequence=it["sequence"], cath_label=None, level="-",
                                           source="none", notes=list(self.notes)) for it in items}
        # exact
        if self.md5 is not None:
            for it in items:
                results[it["name"]].ted_exact = self.md5.lookup(it["sequence"])
        # structure
        with_pdb = [it for it in items if it.get("pdb")]
        fs = {}
        if self.foldclass is not None and with_pdb:
            qdir = tmp / "domains"
            qdir.mkdir(exist_ok=True)
            for it in with_pdb:
                shutil.copy(it["pdb"], qdir / f"{it['name']}.pdb")
            for h in foldseek_search(qdir, tmp, threads=self.threads):
                fs.setdefault(h["query"], []).append(h)
            for it in with_pdb:
                r = results[it["name"]]
                a = ted_foldseek_rule(fs.get(it["name"], []), self.cath_labels)
                if a is None:
                    a = self.foldclass.search(qdir / f"{it['name']}.pdb", self.cath_labels, tmp)
                r.structure = a
        elif "structure" in self.tiers:
            for r in results.values():
                r.notes.append("no structure given: structure tier skipped")
        # sequence
        if "sequence" in self.tiers:
            from .server import remote_search, server_url
            try:
                seqs = {it["name"]: it["sequence"] for it in items}
                url = server_url() if self.use_server else None
                if url:
                    hits = remote_search(url, seqs)
                elif _gpu_available() and Path(str(self.mmseqs_db or config.TED_MMSEQS_DB) + "_pad.dbtype").exists():
                    # no server, but we are on a GPU node: one-off MMseqs2-GPU search (loads the 15 GB padded DB)
                    hits = mmseqs_search(seqs, tmp, threads=self.threads, db=self.mmseqs_db, gpu=True)
                elif self.mmseqs_db is not None or os.environ.get("TED_MMSEQS_LOCAL_CPU") == "1":
                    hits = mmseqs_search(seqs, tmp, threads=self.threads, db=self.mmseqs_db,
                                         db_load_mode=self.db_load_mode)
                else:
                    raise FileNotFoundError(
                        "sequence tier skipped: no MMseqs server is running and no GPU is available. Start one with "
                        "`sbatch scripts/mmseqs_server.sbatch` (or set TED_MMSEQS_LOCAL_CPU=1 to search locally on "
                        "CPU, which loads a ~100 GB index)")
                for name, hs in hits.items():
                    results[name].sequence_hits = sorted(hs, key=lambda h: -h["bits"])[:25]
            except FileNotFoundError as e:
                for r in results.values():
                    r.notes.append(str(e))
        # decide: TED's own label if the exact domain is in TED, else structure (TED's method), else sequence
        for r in results.values():
            exact_lab = [e for e in r.ted_exact if e["level"] != "-"]
            r.sequence_transfer = sequence_transfer(r.sequence_hits)     # always reported, used last
            if exact_lab:
                e = max(exact_lab, key=lambda e: e["level"] == "H")
                r.cath_label, r.level, r.source = e["cath_label"], e["level"], "ted-exact"
            elif r.structure and r.structure.get("cath_label"):
                r.cath_label, r.level, r.source = r.structure["cath_label"], r.structure["level"], r.structure["method"]
            else:
                if r.sequence_transfer and r.sequence_transfer["cath_label"]:
                    t = r.sequence_transfer
                    r.cath_label, r.level, r.source = t["cath_label"], t["level"], "mmseqs-transfer"
                elif r.sequence_transfer:
                    r.notes.append("joins a TED sequence cluster that TED left unlabelled")
            if r.source != "mmseqs-transfer":
                labelled = [h for h in r.sequence_hits if h["cluster_label"] != "-"]
                if labelled and r.cath_label is None:
                    h = labelled[0]
                    r.remote_hint = {"cath_label": h["cluster_label"], "level": h["cluster_level"],
                                     "ted_rep": h["ted_rep"], "fident": h["fident"], "qcov": h["qcov"],
                                     "evalue": h["evalue"],
                                     "note": "homolog below TED's cluster criterion; not used as the label"}
            if r.ted_exact and not exact_lab:
                r.notes.append("sequence is identical to TED domain(s) that TED left unlabelled")
            if r.cath_label:
                r.name_of_label = self.cath_names.get(r.cath_label)
        return list(results.values())


def cath_label(sequence: str | None = None, *, pdb: str | Path | None = None, name: str = "query",
               tiers=("exact", "structure", "sequence"), **kw) -> LabelResult:
    """CATH label for one domain. Give its sequence, and/or a PDB of the domain (chain A) for the structure tier."""
    if sequence is None and pdb is None:
        raise ValueError("give sequence= and/or pdb=")
    if sequence is None:
        from .structure import read_pdb_ca
        sequence = read_pdb_ca(pdb)[0]
    item = {"name": name, "sequence": "".join(sequence.split()).upper()}
    if pdb:
        item["pdb"] = str(pdb)
    return Classifier(tiers=tiers, **kw).label([item])[0]


def cath_label_chopped(chop_results, tiers=("exact", "structure", "sequence"), workdir=None, **kw):
    """Label every TED domain of ted_chop() results; domain structures are cut from the chain models."""
    from .chop import parse_chopping
    tmp = Path(workdir) if workdir else Path(tempfile.mkdtemp(prefix="ted_cath_", dir=os.environ.get("TMPDIR")))
    tmp.mkdir(parents=True, exist_ok=True)
    items = []
    for r in chop_results:
        for d in r.domains:
            it = {"name": d.ted_id, "sequence": d.sequence}
            if r.structure and Path(r.structure).exists():
                it["pdb"] = str(write_domain_pdb(r.structure, parse_chopping(d.chopping)[0], tmp / f"{d.ted_id}.pdb"))
            items.append(it)
    try:
        return Classifier(tiers=tiers, **kw).label(items, workdir=tmp / "work")
    finally:
        if workdir is None:
            shutil.rmtree(tmp, ignore_errors=True)
