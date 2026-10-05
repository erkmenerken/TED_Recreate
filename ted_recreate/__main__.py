"""Command line for the TED re-creation.

  python -m ted_recreate chop  --fasta proteins.fa  [--out chop.json]      sequence -> ESMFold -> TED domains
  python -m ted_recreate chop  --uniprot P69905 Q9XYZ1                     AFDB model (what TED used) -> domains
  python -m ted_recreate chop  --pdb model.pdb                             your own model -> domains
  python -m ted_recreate cath  --fasta domains.fa  [--tiers exact,sequence] domain sequence -> CATH label
  python -m ted_recreate cath  --pdb domain.pdb                             domain structure -> CATH label
  python -m ted_recreate label --fasta proteins.fa                          chop, then label every domain

Run chop/label on a GPU node (ESMFold, Merizo, Chainsaw); see scripts/run_ted.sbatch.
"""
import argparse
import json
import os
import sys
from pathlib import Path


def read_fasta(path):
    recs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            recs[name] = ""
        elif name:
            recs[name] += line
    return recs


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m ted_recreate", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["chop", "cath", "label"])
    ap.add_argument("--fasta")
    ap.add_argument("--seq", help="a single sequence")
    ap.add_argument("--uniprot", nargs="+")
    ap.add_argument("--pdb", nargs="+")
    ap.add_argument("--tiers", default="exact,structure,sequence")
    ap.add_argument("--workdir", help="keep intermediate files (models, parser outputs) here")
    ap.add_argument("--out", help="JSON output (default: stdout summary only)")
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args(argv)

    items = []
    if a.fasta:
        items += [{"name": k, "sequence": v} for k, v in read_fasta(a.fasta).items()]
    if a.seq:
        items.append({"name": "query", "sequence": a.seq})
    for acc in a.uniprot or []:
        items.append({"uniprot": acc})
    for p in a.pdb or []:
        items.append({"pdb": p, "name": Path(p).stem})
    if not items:
        ap.error("give --fasta, --seq, --uniprot or --pdb")

    tiers = tuple(t for t in a.tiers.split(",") if t)
    out = {}
    tmp_work = None
    if a.command == "label" and not a.workdir:
        # the CATH structure tier cuts domains out of the chain models, so keep them until labelling is done
        import tempfile
        tmp_work = a.workdir = tempfile.mkdtemp(prefix="ted_label_", dir=os.environ.get("TMPDIR"))
    if a.command in ("chop", "label"):
        from ted_recreate.chop import ted_chop_batch
        res = ted_chop_batch(items, workdir=a.workdir, keep_workdir=bool(a.workdir))
        out["chop"] = [r.to_dict() for r in res]
        for r in res:
            doms = ", ".join(f"{d.ted_id.rsplit('_', 1)[-1]}:{d.chopping}({d.consensus_level[0]})" for d in r.domains)
            print(f"{r.name}\t{r.nres} aa\t{r.structure_source}\t{doms or 'no TED domains'}")
            for w in r.warnings:
                print(f"  WARNING {r.name}: {w}", file=sys.stderr)
        if a.command == "label":
            from ted_recreate.classify import cath_label_chopped
            labs = cath_label_chopped(res, tiers=tiers, threads=a.threads,
                                      workdir=Path(a.workdir) / "cath" if a.workdir else None)
            out["cath"] = [l.to_dict() for l in labs]
            for l in labs:
                print(f"{l.name}\t{l.cath_label or '-'}\t{l.level}\t{l.source}\t{l.name_of_label or ''}")
    else:
        from ted_recreate.classify import Classifier
        from ted_recreate.structure import read_pdb_ca
        cls_items = []
        for it in items:
            if "pdb" in it:
                cls_items.append({"name": it["name"], "sequence": read_pdb_ca(it["pdb"])[0], "pdb": it["pdb"]})
            elif "sequence" in it:
                cls_items.append({"name": it["name"], "sequence": "".join(it["sequence"].split()).upper()})
            else:
                ap.error("cath takes domain sequences (--fasta/--seq) or domain structures (--pdb)")
        labs = Classifier(tiers=tiers, threads=a.threads).label(cls_items, workdir=a.workdir)
        out["cath"] = [l.to_dict() for l in labs]
        for l in labs:
            print(f"{l.name}\t{l.cath_label or '-'}\t{l.level}\t{l.source}\t{l.name_of_label or ''}")
    if a.out:
        Path(a.out).write_text(json.dumps(out, indent=1, default=str))
        print(f"wrote {a.out}", file=sys.stderr)
    if tmp_work:
        import shutil
        shutil.rmtree(tmp_work, ignore_errors=True)


if __name__ == "__main__":
    main()
    # After ESMFold/CUDA work the interpreter can hang in teardown (seen in slurm jobs); results are written, so
    # flush and leave without waiting for library threads.
    import os
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(0)
