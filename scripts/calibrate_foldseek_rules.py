"""Which Foldseek columns implement TED's 'TM-score' and 'coverage' thresholds?

TED (Supp. Methods): H if E < 0.019, TM > 0.56, coverage >= 0.367; T if E < 0.108662, TM >= 0.42,
coverage >= 0.786 -- without naming the Foldseek output columns. We cut TED domains (published chopping)
out of the AFDB models fetched for validation, search them with TED's search settings against CATH 4.3 S40,
apply every column combination, and score agreement with TED's own Foldseek labels.

Usage: python scripts/calibrate_foldseek_rules.py SUMMARY_TSV OUTDIR AFDB_DIR [AFDB_DIR ...]
"""
import collections
import itertools
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ted_recreate import config
from ted_recreate.chop import parse_chopping
from ted_recreate.classify import _cath_labels, foldseek_search, ted_foldseek_rule, write_domain_pdb

summary, outdir, afdb_dirs = sys.argv[1], Path(sys.argv[2]), [Path(p) for p in sys.argv[3:]]
outdir.mkdir(parents=True, exist_ok=True)
models = {p.stem: p for d in afdb_dirs for p in d.glob("*.pdb")}

# TED-100 domains of chains we have models for (the TED-redundant set was labelled differently: skip ids in it)
redundant = set()
red_list = config.DATA / "ted" / "ted_redundant_40m_domain_id.list.gz"
ted = {}
for line in open(summary):
    f = line.rstrip("\n").split("\t")
    chain = f[0].rsplit("_TED", 1)[0]
    if chain in models:
        ted[f[0]] = {"chain": chain, "chopping": f[3], "label": f[13], "level": f[14], "method": f[15]}
import gzip
with gzip.open(red_list, "rt") as fh:
    for line in fh:
        if line.strip() in ted:
            redundant.add(line.strip())
ted = {k: v for k, v in ted.items() if k not in redundant}
print(f"{len(ted)} TED-100 domains from {len(models)} models", flush=True)

qdir = outdir / "domains"
qdir.mkdir(exist_ok=True)
for tid, d in ted.items():
    write_domain_pdb(models[d["chain"]], parse_chopping(d["chopping"])[0], qdir / f"{tid}.pdb")

cath = _cath_labels()
results = {}
for name, binary in (("foldseek8", config.FOLDSEEK_TED), ("foldseek_latest", config.FOLDSEEK)):
    (outdir / name).mkdir(exist_ok=True)
    hits = collections.defaultdict(list)
    for h in foldseek_search(qdir, outdir / name, threads=8, foldseek=binary):
        hits[h["query"]].append(h)
    for tm_col, cov_col in itertools.product(["alntmscore", "qtmscore", "ttmscore", "maxtmscore"],
                                             ["lenratio", "qcov", "tcov", "mincov"]):
        c = collections.Counter()
        for tid, d in ted.items():
            a = ted_foldseek_rule(hits.get(tid, []), cath, tm_col=tm_col, cov_col=cov_col)
            truth = (d["level"], d["label"]) if d["method"] == "foldseek" else ("-", "-")
            pred = (a["level"], a["cath_label"]) if a else ("-", "-")
            c["n"] += 1
            c["exact"] += truth == pred
            c["level_match"] += truth[0] == pred[0]
            if truth[0] == "H":
                c["ted_H"] += 1
                c["H_recovered"] += pred == truth
            if truth[0] == "T":
                c["ted_T"] += 1
                c["T_recovered"] += pred == truth
            if truth[0] == "-":
                c["ted_none"] += 1
                c["false_label"] += pred[0] != "-"
        results[f"{name} tm={tm_col} cov={cov_col}"] = dict(c)
        print(f"{name:16s} tm={tm_col:10s} cov={cov_col:8s} exact={c['exact']}/{c['n']}  "
              f"H {c['H_recovered']}/{c['ted_H']}  T {c['T_recovered']}/{c['ted_T']}  "
              f"labelled-but-TED-not-by-foldseek {c['false_label']}/{c['ted_none']}", flush=True)
json.dump(results, open(outdir / "calibration.json", "w"), indent=1)
