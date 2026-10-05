"""Run our TED chopping on AFDB models of random TED proteins and compare with the published TED domains.

Usage: python scripts/validate_chop_vs_ted.py SUMMARY_TSV OUTDIR [N_RANDOM] [N_START1]
SUMMARY_TSV: rows of ted_365m.domain_summary (a prefix is fine; the last chain is dropped as possibly incomplete).
"""
import hashlib, json, random, sys, collections
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ted_recreate.chop import ted_chop_batch, parse_chopping, extract
from ted_recreate.structure import fetch_afdb_model, read_pdb_ca

summary, outdir = sys.argv[1], Path(sys.argv[2])
n_random = int(sys.argv[3]) if len(sys.argv) > 3 else 150
n_start1 = int(sys.argv[4]) if len(sys.argv) > 4 else 40
outdir.mkdir(parents=True, exist_ok=True)

chains = collections.OrderedDict()
for line in open(summary):
    f = line.rstrip("\n").split("\t")
    chain = f[0].rsplit("_TED", 1)[0]
    chains.setdefault(chain, []).append({"ted_id": f[0], "md5_domain": f[1], "level": f[2], "chopping": f[3],
                                         "cath": f[13], "cath_level": f[14]})
chains.popitem()  # possibly truncated
rng = random.Random(0)
names = list(chains)
start1 = [c for c in names if any(d["chopping"].startswith("1-") for d in chains[c])]
pick = rng.sample(names, n_random) + rng.sample(start1, n_start1)
pick = list(dict.fromkeys(pick))

def get(chain):
    acc = chain.split("-")[1]
    try:
        pdb, rec = fetch_afdb_model(acc, outdir / "afdb")
        seq = read_pdb_ca(pdb)[0]
        ok = all(hashlib.md5(extract(seq, parse_chopping(d["chopping"])[0]).encode()).hexdigest() == d["md5_domain"]
                 for d in chains[chain])
        return chain, pdb, ok, rec.get("modelCreatedDate")
    except Exception as e:
        return chain, None, False, repr(e)[:80]

with ThreadPoolExecutor(8) as ex:
    fetched = list(ex.map(get, pick))
usable = [(c, p) for c, p, ok, _ in fetched if p is not None and ok]
print(f"picked {len(pick)}; AFDB model with identical sequence for {len(usable)}", flush=True)
json.dump([{"chain": c, "pdb": str(p) if p else None, "seq_ok": ok, "info": i} for c, p, ok, i in fetched],
          open(outdir / "fetched.json", "w"), indent=1)

res = ted_chop_batch([{"pdb": str(p), "name": c} for c, p in usable], workdir=outdir / "work", keep_workdir=True)
rows = []
for r in res:
    ted = sorted((d["level"], d["chopping"]) for d in chains[r.name])
    ours = sorted((d.consensus_level, d.chopping) for d in r.domains)
    ted_set = {c for _, c in ted}
    rows.append({
        "chain": r.name, "start1": r.name in start1,
        "exact": ours == ted,
        "same_ranges_any_level": sorted(c for _, c in ours) == sorted(ted_set),
        "same_ranges_shift_minus1": sorted(d.chopping_parser_frame for d in r.domains) == sorted(ted_set),
        "n_ted": len(ted), "n_ours": len(ours), "ted": ted, "ours": ours, "methods": r.methods,
    })
json.dump(rows, open(outdir / "comparison.json", "w"), indent=1)
for grp, sel in (("random", [x for x in rows if not x["start1"]]), ("start-at-1", [x for x in rows if x["start1"]])):
    if not sel:
        continue
    n = len(sel)
    print(f"[{grp}] n={n}  exact={sum(x['exact'] for x in sel)}  same ranges={sum(x['same_ranges_any_level'] for x in sel)}"
          f"  match if our output shifted -1={sum(x['same_ranges_shift_minus1'] for x in sel)}"
          f"  same domain count={sum(x['n_ted']==x['n_ours'] for x in sel)}")
