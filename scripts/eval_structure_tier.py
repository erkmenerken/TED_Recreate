"""Run the full structure tier (Foldseek with TED rules, then Foldclass) on the calibration domains and
compare the final label with TED's published label (Foldseek H/T or Foldclass T)."""
import collections, json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ted_recreate.classify import Classifier

summary, domdir = sys.argv[1], Path(sys.argv[2])
ted = {}
for line in open(summary):
    f = line.rstrip("\n").split("\t")
    if (domdir / f"{f[0]}.pdb").exists():
        ted[f[0]] = {"label": f[13], "level": f[14], "method": f[15]}
from ted_recreate.structure import read_pdb_ca
items = [{"name": t, "sequence": read_pdb_ca(domdir / f"{t}.pdb")[0], "pdb": str(domdir / f"{t}.pdb")} for t in ted]
res = Classifier(tiers=("structure",), threads=8).label(items, workdir=Path(sys.argv[3]))
c = collections.Counter(); rows = []
for r in res:
    t = ted[r.name]
    truth = (t["level"], t["label"] if t["level"] != "-" else None, t["method"])
    pred = (r.level, r.cath_label, r.source)
    key = f"TED {t['method']}-{t['level']}"
    c[(key, "exact")] += (truth[0], truth[1]) == (pred[0], pred[1])
    c[(key, "same_T")] += bool(truth[1] and pred[1] and truth[1].split('.')[:3] == pred[1].split('.')[:3])
    c[(key, "n")] += 1
    c[("all", "exact")] += (truth[0], truth[1]) == (pred[0], pred[1]); c[("all", "n")] += 1
    rows.append({"ted_id": r.name, "ted": truth, "ours": pred, "structure": r.structure})
for k in sorted(c): print(k, c[k])
json.dump(rows, open(Path(sys.argv[3]) / "structure_tier_eval.json", "w"), indent=1, default=str)
