# TED_Recreate

A re-implementation of how TED (The Encyclopedia of Domains; Lau et al., Science 2024) chops proteins into
domains and assigns CATH labels, built on TED's own released code. It is meant for producing TED-style labels for
new sequences and for checking our domain-boundary and CATH models against TED.

- `docs/presentation_review.md`: how TED and CATH work, with corrections to the presentation script.
- `ted_recreate/`: the Python package. The two main functions are `ted_chop` and `cath_label`.
- `scripts/`: builders for the lookup databases, validation scripts and Slurm wrappers.

Python: `.venv/bin/python`, which is the conda env `ted` (torch 2.10, transformers 5.3). Run GPU work through
Slurm; the header of `scripts/run_ted.sbatch` has the settings.

## 1. Sequence → TED domains: `ted_chop`

```python
import sys; sys.path.insert(0, "/scratch/erkmenerken22/TED_Recreate")
from ted_recreate import ted_chop, ted_chop_batch

r = ted_chop(sequence="MKT...")          # folded with ESMFold, then chopped
r = ted_chop(uniprot="A0A1X0B8H0")      # AFDB model: the same input TED used
r = ted_chop(pdb="my_model.pdb")         # any single-chain model
for d in r.domains:
    print(d.ted_id, d.consensus_level, d.chopping, d.chopping_parser_frame, d.nres, d.plddt)
r.methods            # the three parsers' outputs: {'merizo': ..., 'chainsaw': ..., 'unidoc': ...}
r.residue_labels()   # per-residue domain index (0 = not in a TED domain), handy for training targets
```

What happens, following TED's released pipeline (`vendor/ted-tools/ted_consensus_1.0`, run unmodified apart
from the UniDoc build):

1. **Structure.** An AFDB model when you give a UniProt accession. AFDB now serves v6 files, whose coordinates are
   identical to v4 for 2022 models. Otherwise the sequence is folded with **ESMFold**.
2. **Merizo** (dual pass, iterative re-segmentation, UniDoc-style merge).
3. **UniDoc** on the residues Merizo put in domains, using the **2022 `UniDoc_structure` build**, which is the
   one that reproduces TED's published UniDoc outputs.
4. **Chainsaw**, with a +1 residue offset because Chainsaw is 0-based.
5. Each parser's output is filtered: segments < 5 residues and domains < 25 residues are dropped.
6. **Consensus.** Domains are linked when their IoU ≥ 0.7. Connected components of size 3/2/1 are
   high/medium/low. The range kept is the residues shared by the whole component. Only high and medium domains
   become TED domains, named TED01, TED02, … from the N-terminus.

Two boundary frames are reported:
- `chopping` is in **TED's convention**. It is +1 residue relative to what the parsers agreed, which is what the
  released code and the published TED-100 database contain.
- `chopping_parser_frame` undoes that shift.

`sequence` and `md5_domain` use `chopping`, just as TED's `md5_domain` does.

**How faithful is it?** On AFDB models of 209 random TED-100 chains, **202 chains (96.7%) and 381/387 domains
(98.4%)** come out exactly as published (`runs/validate/out2`). All remaining differences come from parser
builds: Chainsaw differs in 4 chains, UniDoc in 9, Merizo in 1. Whenever all three parser outputs match TED's
(196 chains), the consensus matches 196/196.

TED-redundant chains were made with a different consensus function (see `docs/presentation_review.md` §3), so
`ted_chop` reproduces TED-100, not TED-redundant.

**The one step that cannot match TED is folding a bare sequence.** ESMFold is not AlphaFold2, so the domains
can differ, especially in low-confidence regions. For UniProt proteins, prefer `uniprot=`.

## 2. Domain → CATH label: `cath_label`

```python
from ted_recreate.classify import cath_label, cath_label_chopped, Classifier
lab = cath_label("DOMAINSEQUENCE...")                  # sequence tiers only
lab = cath_label(pdb="domain.pdb")                     # + structure tier
labs = cath_label_chopped([r])                         # every domain of ted_chop() results
lab.cath_label, lab.level, lab.source, lab.name_of_label
```

The tiers follow how TED assigned labels (Supp. Methods). The final label comes from the first tier that
gives one:

| Tier | What it does | Data |
|---|---|---|
| `exact` | MD5 of the domain sequence looked up among all 365M TED domains; returns TED's published label | `db/ted_md5/` |
| `structure` | **Foldseek** (release 8, closest to TED's commit) against CATH 4.3 S40 with TED's thresholds: H if E < 0.019, TM > 0.56, cov ≥ 0.367; T if E < 0.108662, TM ≥ 0.42, cov ≥ 0.786. TM = max(qTM, tTM) and cov = length ratio (calibrated, see below). If nothing passes: **Foldclass** nearest neighbour + TM-align > 0.5 (normalised by the query) → T level | `data/cath/foldseek_s40`, `db/cath/` |
| `sequence` | **MMseqs2** against the 120.7M representatives of TED's own sequence clusters. If the query meets TED's clustering rule against a representative (≥ 50% identity, length ratio ≥ 0.9), it would have joined that cluster and gets the cluster's label, which is how TED spread labels. Weaker homologs are reported only as `remote_hint` | `db/mmseqs/ted100` |

Calibration and checks:
- TED searched 31,574 "SSG5" CATH representatives. That list is not public, so we use the CATH 4.3 S40 set
  (31,885 domains), which puts E-values on the same scale.
- On 387 TED-100 domains, Foldseek with these settings reproduces TED's Foldseek outcome for 375 (97%)
  (`runs/calibrate`).
- The whole structure tier reproduces TED's label for 353/387 (91%): H labels 205/207, TED-unlabelled domains
  correctly left unlabelled 115/123.

### The "MMseqs server"

The lookup DB holds the 120,748,700 representatives of TED's sequence clusters (14.6 G residues), each carrying
its cluster's label summary. Searching it cold means loading a ~100 GB k-mer index (CPU) or a 15 GB padded DB
(GPU) for every query. Keeping that in memory is what the "server" is:

```bash
sbatch scripts/build_lookup_dbs.sbatch      # once: db/ted_md5, db/mmseqs/ted100 (+ .idx CPU index, _pad GPU DB)
sbatch scripts/mmseqs_server.sbatch         # CPU server: preloads the index into RAM (touchdb) and answers over HTTP
sbatch --gres=gpu:tesla_a100:1 --mem=64G scripts/mmseqs_server.sbatch gpu   # or: MMseqs2-GPU gpuserver
```

The server writes its address to `db/mmseqs/SERVER`, and `cath_label` uses it automatically. Without a server,
the same search runs locally, just slower. You can also call it directly:
`python -c "from ted_recreate.server import server_url, remote_search; ..."`.

Measured on the full DB in GPU mode (one A100): a batch of 8 domains searched against all 120.7M
representatives in ~17 s, and 11 domains labelled with all three tiers in 28 s. The labels agreed with TED's
for every test domain, and a kinase with 30% of its residues randomised was still assigned 1.10.510.10.

Stop the server with `scancel` when you're done. It holds a lot of RAM (CPU mode) or a GPU (GPU mode).

How the DB was built (`scripts/build_ted_mmseqs_db.py`): TED publishes the cluster table and each domain's
`md5_domain`, but not the sequences. Those come from Foldseek's `teddb` (364,806,077 plain sequence records).
The builder walks `teddb` and the TED summary table together and synchronises them by MD5, so every sequence
is assigned to a TED id whose published `md5_domain` it matches. All 120,748,700 representatives were found,
with 0 unmatched records.

## 3. Command line / Slurm

```bash
python -m ted_recreate chop  --fasta proteins.fa --out chop.json
python -m ted_recreate label --uniprot P00533 --out labels.json          # chop + label every domain
python -m ted_recreate cath  --fasta domains.fa --tiers exact,sequence    # label domain sequences
sbatch scripts/run_ted.sbatch label proteins.fa out.json                  # the same on a GPU node
```

Example (`runs/e2e/`, one A100, ~1 min for AFDB input and ~5 min for two sequences through ESMFold):

```
$ python -m ted_recreate label --uniprot P00533        # EGFR, AFDB model
AF-P00533-F1-model_v4  1210 aa  AFDB v6   TED01:31-207(h), TED02:703-793(m), TED03:796-991(m)
..._TED01  3.80.20.20   H  foldseek  Receptor L-domain
..._TED02  3.30.200.20  H  foldseek  Phosphorylase Kinase; domain 1
..._TED03  1.10.510.10  H  foldseek  Transferase(Phosphotransferase) domain 1
TED database: 31-207 high 3.80.20.20 | 703-793 medium 3.30.200.20 | 796-989 medium 1.10.510.10
```

The same EGFR sequence through ESMFold gives 7 domains, because ESMFold's model differs from AlphaFold's. It
recovers the L domains, the cysteine-rich domains (2.10.220.10) and both kinase lobes.

## 4. Benchmark against TED

`benchmark/REPORT.md` has the full test on 600 random TED-100 proteins (1,037 published domains), with eight
figures in `benchmark/figures/`. In short:

| | AlphaFold model (TED's input) | Sequence only (ESMFold) |
|---|---|---|
| TED domains reproduced exactly | 98.0% | 21.7% |
| TED domains matched at overlap ≥ 0.8 | 98.7% | 79.8% |
| Boundaries within ±8 residues | 99.6% | 84.2% |
| TED-labelled domains recovered and given TED's label | 99.2% | 81.1% |
| ... the same, scored as if the protein were new to TED | 93.9% | 79.2% |

## 5. Layout

```
ted_recreate/        chop.py (ted_chop), structure.py (AFDB fetch, ESMFold, PDB normalisation),
                     classify.py (cath_label), server.py (MMseqs server), __main__.py (CLI)
vendor/ted-tools     TED's released code; programs/unidoc uses the 2022 UniDoc_structure build
                     (set TED_UNIDOC_BIN=UniDoc_struct for the released behaviour)
vendor/merizo_search Foldclass network + TM-align used by the structure tier
tools/bin            mmseqs, mmseqs-gpu, foldseek, foldseek8, TMalign
data/                TED Zenodo tables, Foldseek TED sequence DB, CATH 4.3 files
db/                  built lookups (cath/, ted_md5/, mmseqs/)
runs/                validation (validate/), calibration (calibrate/), end-to-end test (e2e/)
benchmark/           600-protein benchmark against TED: scripts, results/, figures/, REPORT.md
```

## 6. Setup on a new machine

The repository holds only our code, the benchmark and the docs. Not versioned (see `.gitignore`):

- `vendor/ted-tools`: `git clone https://github.com/psipred/ted-tools`, add the Merizo weights as its `setup.sh`
  does, unpack UniDoc into `ted_consensus_1.0/programs/unidoc` and copy `bin/UniDoc_structure` from the UniDoc
  2023 package (https://yanglab.qd.sdu.edu.cn/UniDoc/download/), then apply
  `patches/unidoc_wrapper_use_2022_binary.patch` to the wrapper script.
- `vendor/merizo_search`: https://github.com/psipred/merizo_search (Foldclass network and TM-align).
- `tools/bin`: `mmseqs`, `mmseqs-gpu`, `foldseek`, `foldseek8` (release 8-ef4e960), `TMalign`.
- `.venv`: Python 3.10 with torch, transformers, biopython, pdb-tools, networkx, natsort, einops,
  rotary-embedding-torch, pandas, matplotlib.
- `data/`: TED tables from Zenodo record 13908086, CATH 4.3 files, Foldseek's `teddb` sequence file.
- `db/`: built with `scripts/build_cath_reference.py` and `scripts/build_lookup_dbs.sbatch`.
- `models/hf`: `facebook/esmfold_v1` from Hugging Face.

Paths are set in `ted_recreate/config.py`.
