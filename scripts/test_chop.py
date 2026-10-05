"""Smoke tests: official examples via pdb=, one UniProt accession via AFDB, one raw sequence via ESMFold."""
import sys, json, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ted_recreate import ted_chop, ted_chop_batch
from ted_recreate import config

ex = config.TED_CONSENSUS / "example"
expected = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(config.TED_CONSENSUS / "example_output" / "consensus.tsv")}
t = time.time()
res = ted_chop_batch([{"pdb": str(p)} for p in sorted(ex.glob("*.pdb"))])
for r in res:
    e = expected[r.name]
    same = (r.consensus["high"], r.consensus["medium"], r.consensus["low"]) == (e[6], e[7], e[8])
    print(f"{r.name:28s} identical_to_official_example={same}  domains={[(d.ted_id.split('_')[-1], d.consensus_level, d.chopping) for d in r.domains]}")
print(f"examples: {time.time()-t:.0f}s")

t = time.time()
r = ted_chop(uniprot="A0A1X0B8H0")
print("AFDB A0A1X0B8H0:", r.structure_source, [(d.ted_id, d.consensus_level, d.chopping) for d in r.domains], f"{time.time()-t:.0f}s")

# the same protein from sequence only -> ESMFold model instead of the AlphaFold2 model TED used
t = time.time()
r2 = ted_chop(sequence=r.sequence, name="A0A1X0B8H0_esmfold", workdir=Path(sys.argv[1]) / "esm_test")
print("ESMFold same seq:", r2.structure_source, [(d.ted_id, d.consensus_level, d.chopping) for d in r2.domains], f"{time.time()-t:.0f}s")
print("methods:", json.dumps(r2.methods))
print("TED DB for A0A1X0B8H0: medium 3-209, 216-499, 525-686")
