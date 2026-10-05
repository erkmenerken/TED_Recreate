## Numbers

Benchmark set: **600 TED-100 proteins** with **1,037 published TED domains**. Proteins that failed: 0 (AlphaFold-model mode), 0 (sequence-only mode).

### Domain chopping

| | AlphaFold model (TED's input) | Sequence only (ESMFold) |
|---|---|---|
| Proteins with every domain identical to TED | 96.8% | 10.3% |
| Proteins with the same number of domains | 98.8% | 82.0% |
| TED domains reproduced exactly | 98.0% | 21.7% |
| TED domains matched at overlap ≥ 0.9 | 98.6% | 68.5% |
| TED domains matched at overlap ≥ 0.8 | 98.7% | 79.8% |
| TED domains matched at overlap ≥ 0.5 | 99.3% | 91.5% |
| Our domains with a TED counterpart (overlap ≥ 0.8) | 99.2% | 78.9% |
| Our domains with no TED counterpart (overlap < 0.5) | 0.2% | 9.5% |
| Residues with the same in-domain / not-in-domain status | 99.1% | 89.3% |
| Boundaries identical to TED's | 99.2% | 47.1% |
| Boundaries within ±2 residues | 99.5% | 69.3% |
| Boundaries within ±8 residues | 99.6% | 84.2% |
| Boundaries within ±20 residues | 99.7% | 90.9% |
| Domains produced | 1,032 | 1,049 |
| Proteins where UniDoc failed | 0 | 0 |

ESMFold models vs AlphaFold models: median TM-score 0.85; 40.5% of proteins at TM ≥ 0.9, 16.8% below 0.5. Median ESMFold mean pLDDT 87.

### CATH labels, each route on TED's own published domains

| Route | Same label as TED | Same fold or better | Different label | No label | Correct when it gives a label | TED-unlabelled domains we label |
|---|---|---|---|---|---|---|
| exact lookup | 100.0% | 100.0% | 0.0% | 0.0% | 100.0% | 0.3% |
| structure | 93.5% | 94.6% | 3.1% | 2.3% | 95.7% | 6.4% |
| sequence | 97.3% | 98.1% | 0.8% | 1.1% | 98.4% | 9.8% |
| sequence (own cluster removed) | 83.9% | 85.3% | 1.5% | 13.2% | 96.7% | 8.4% |
| nearest labelled homolog (own cluster removed) | 91.2% | 93.1% | 4.3% | 2.6% | 93.6% | 51.9% |

Base: 740 TED-labelled domains (591 with a superfamily, 149 with a fold only) and 297 TED-unlabelled domains.

| Route | TED superfamily labels reproduced | TED fold-only labels reproduced |
|---|---|---|
| exact lookup | 591/591 (100.0%) | 149/149 (100.0%) |
| structure | 579/591 (98.0%) | 113/149 (75.8%) |
| sequence | 579/591 (98.0%) | 141/149 (94.6%) |
| sequence (own cluster removed) | 532/591 (90.0%) | 89/149 (59.7%) |
| nearest labelled homolog (own cluster removed) | 568/591 (96.1%) | 107/149 (71.8%) |

### End to end (chop, then label)

| | AlphaFold model (TED's input) | Sequence only (ESMFold) |
|---|---|---|
| TED-labelled domains recovered (overlap ≥ 0.8) | 99.2% | 85.8% |
| ... recovered **and** given TED's label | 99.2% | 81.1% |
| Label agreement among recovered TED-labelled domains | 100.0% | 94.5% |
| As if new to TED*: recovered **and** given TED's label | 93.9% | 79.2% |
| As if new to TED*: recovered and same fold or better | 95.0% | 81.6% |
| As if new to TED*: label agreement among recovered domains | 94.7% | 92.3% |
| Our domains decided by: exact lookup | 70.5% | 19.4% |
| Our domains decided by: Foldseek | 1.2% | 43.3% |
| Our domains decided by: Foldclass | 1.6% | 8.8% |
| Our domains decided by: sequence cluster | 2.3% | 3.9% |
| Our domains decided by: no label | 24.4% | 24.7% |

\* exact lookup switched off and the domain's own TED sequence cluster removed from the MMseqs2 database.

### Run time (one A100, 8 CPU cores)

| Stage | Total | Per protein or domain |
|---|---|---|
| Chop from AlphaFold models (600 proteins) | 20 min | 2.0 s per protein |
| ESMFold + chop from sequence (600 proteins) | 124 min | 12.4 s per protein |
| CATH labels, all three routes: TED's published domains (1037 domains) | 11.9 min | 0.69 s per domain |
| CATH labels, all three routes: our domains, AlphaFold models (1032 domains) | 11.3 min | 0.66 s per domain |
| CATH labels, all three routes: our domains, ESMFold models (1049 domains) | 11.5 min | 0.66 s per domain |
