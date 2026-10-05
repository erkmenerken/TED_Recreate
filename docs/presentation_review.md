# Review of the TED / CATH presentation script

Sources: Lau et al. 2024 (Science 386, eadq4946) and its Supplementary Materials; the Merizo, Chainsaw, UniDoc and
CATH (Orengo et al. 1997) papers in `Papers/`; TED's released code (`vendor/ted-tools`); the TED data on Zenodo
(record 13908086). Numbers were checked against Supp. Table S1. Where the paper and the code disagree, the code
(and the published database) wins, and that is flagged.

## 1. What to fix, claim by claim

| # | Script says | Verdict | Correct version |
|---|---|---|---|
| 1 | AFDB v4, ~214.7M predicted structures from UniProt | Correct | 214,683,829 models (UniProt release 2021_04). |
| 2 | Redundant sequences filtered → TED-100, 188.9M | Correct | 188,914,411 unique sequences. 38.9M targets share their sequence with another target; 13.2M of those are the kept representatives (the first target in a sorted list), and the other 25.8M form **TED-redundant**. |
| 3 | Redundant copies kept because identical sequences can have different AF models | Correct, add the evidence | Fig. S16: ~42% of identical-sequence groups (5.6M) have a max RMSD > 1 Å, up to ~65 Å. TED-redundant was chopped separately (40.4M domains). The CATH classification and the sequence clustering use **TED-100 only**. |
| 4 | "chopped these **sequences** into domains" | **Wrong word** | The three parsers chop the **3D models**. Every TED boundary is derived from the AlphaFold2 structure. This matters for our project: a sequence model is being trained to predict structure-derived boundaries. |
| 5 | 324.4M accepted TED-100 domains; +40.4M redundant = 364.8M | Correct | 324,389,697 = 195,229,271 high + 129,160,426 medium. Supp. Table S1 has a typo ("high 250,629,037"), and the Supp. text has an inconsistent 324,482,131. |
| 6 | Two operations in parallel | Correct | (a) MMseqs2 sequence clustering and (b) structure-based matching to CATH. The labels from (b) are then spread over the clusters from (a). |
| 7 | Cluster if ≥ 50% identity and alignment covers ≥ 90% of the shorter sequence | Nearly | The command is `mmseqs easy-linclust --min-seq-id 0.5 -c 0.9 --cov-mode 5`. cov-mode 5 means "the shorter sequence is ≥ 90% of the length of the longer one". This gives 120,748,700 clusters, ~81.7M of them singletons. |
| 8 | CATH = database of experimentally determined domains by structure and evolution | Correct | Semi-automatic classification of PDB domains. CATH 4.3 has 500,238 domains, 5 classes, 43 architectures, 1,472 topologies (folds) and 6,631 superfamilies. |
| 9 | C: mainly α, mainly β, α/β, few secondary structures | Correct | Classes 1–4, plus class 6 "special" in current releases (there is no class 5). |
| 10 | A: 3D arrangement of SSEs ignoring connectivity | Correct | Historically the one level assigned manually. |
| 11 | T: connectivity; "essentially the fold level" | Correct | CATH itself uses fold and topology interchangeably. |
| 12 | H: same label when evidence suggests common ancestry | Correct but vague | The evidence is high structural similarity **plus** sequence and/or functional similarity. |
| 13 | CAT vs CATH: similar structures ≠ common ancestry; proteins can evolve to fold similarly independently | Slightly overstated | Same T but different H means "**no convincing evidence** of homology". The similarity may be convergent, or the proteins may be homologs too divergent to prove it. CATH does not claim they are unrelated. |
| 14 | Worked example of a code | Add one | 3.40.50.300 = α/β (3) · 3-layer αβα sandwich (40) · Rossmann fold (50) · P-loop NTPase superfamily (300). |
| 15 | TED used CATH 4.3 | Correct | — |
| 16 | Foldseek turns 3D structure into a sequence-like representation | Correct | Foldseek's 3Di structural alphabet. |
| 17 | Foldseek searched "representative CATH domain structures" | Be specific | The targets were **31,574 CATH SSG5 representatives** (structural similarity groups: superfamily members within 5 Å). Search settings: `-e 0.108662 -c 0.366757 --cov-mode 5 -s 10`. |
| 18 | "very strong structural evidence → H, fold only → T" | Correct, add the thresholds | **H:** E < 0.019, TM > 0.56, coverage ≥ 0.367. **T:** E < 0.108662, TM ≥ 0.42, coverage ≥ 0.786. Thresholds were fitted with a genetic algorithm for **98% precision** on a SCOP-consistent CATH benchmark; recall is 0.59 (H) and 0.71 (T). |
| 19 | (missing) | **Add** | A second structural step, **Merizo-search / Foldclass**, labels domains Foldseek could not. An equivariant GNN embeds each domain, the nearest CATH domain is found by cosine similarity, and the match is confirmed if TM-align > 0.5 (normalised by the TED domain). This gives a **T-level label only**, and it is where most T labels come from: ~29.8M of 45.8M. |
| 20 | (missing, optional) | Add | Validation: CATH HMMs confirmed 88.5% of Foldseek's H assignments at superfamily level. |
| 21 | "194M H + 46M T + 73M not classified" | **Mixes two stages** (they sum to 313M, not 324M) | Direct structural labels: **193.9M H + 45.8M T** (16.0M Foldseek + 29.8M Foldclass) **+ 84.7M none = 324.4M**. The sequence clusters then add 11.6M: unlabelled domains in clusters with a labelled member inherit the label. After that, **251.3M domains (77%, in 78.9M clusters) are CATH-assignable** and **73.1M domains (in 41.9M clusters) are not**. |
| 22 | Unlabelled ≠ novel; filtered by quality, searched other DBs | Correct, add the specifics | The filters were applied to the representatives of the 41.9M unlabelled clusters: globularity (normalised Rg < 0.356 **and** packing density > 10.333, the 5th percentiles of H-labelled domains), ≥ 6 secondary-structure elements (STRIDE), and pLDDT80 ≥ 90. That leaves 8.6M clusters, which were clustered with Foldseek and searched against PDB, CATH, ECOD and SCOPe (TM > 0.56, 60% coverage). High-symmetry domains (SymD Z > 9) were set aside: **6,433 clusters**. A chopping-quality network, the DOM parser and final TM-align runs then leave **7,427 novel-fold clusters**. |
| 23 | Merizo: sequence + 3D coordinates, IPA, predicts domain membership and NDRs | Correct, refine | Inputs: one-hot sequence, Cα distance map, backbone frames and residue index. **No pLDDT.** An IPA encoder feeds a masked-transformer decoder that assigns residues to domains; a separate head predicts NDRs; a pIoU confidence is output. Trained on CATH 4.3 and fine-tuned on AFDB. Orientation invariance is a property of IPA in general, so cite AlphaFold2 for it. TED's version adds a second pass without NDRs, iterative re-segmentation of large domains, and UniDoc-style merging. |
| 24 | Chainsaw: pairwise co-membership probabilities, then "a global sequence partition" | Fix the second half | The ResNet input is the distance map plus STRIDE secondary structure; the output is an L×L "same domain" probability matrix. Domains are then assigned by a **greedy** search that best agrees with that matrix. It is not a guaranteed global optimum, and residues may stay unassigned (NDRs). Discontinuous domains are handled naturally. |
| 25 | UniDoc: distances + SS, top-down split then bottom-up merge | Correct, add that it is not ML | Cβ distances → contact probabilities → domain interaction scores. It splits recursively, never cutting inside a helix or strand, then merges fragments. UniDoc **cannot detect NDRs**, so TED runs it only on the residues Merizo put into domains. |
| 26 | Size filter: ≥ 25 amino acids | Incomplete | Segments < 5 residues are dropped **and** domains < 25 residues are dropped. The filter runs on each parser's output before consensus and again on the consensus. |
| 27 | Agreement = "overlap across at least 70% of their residue ranges" | **Imprecise** | The code uses **intersection-over-union ≥ 0.7**: shared residues divided by residues in either domain. A domain that sits inside a much larger one does not count as agreeing. |
| 28 | Graph, components of 3 / 2 / 1 → high / medium / low; only high and medium kept | Correct | High and medium domains cannot overlap. |
| 29 | Final boundary = intersection of the agreeing predictions | Correct as described in the paper | In the released code, and in the published TED-100 database, every consensus range is **shifted one residue towards the C-terminus**: `domstr_to_assignment(..., offset=1)`. No TED-100 domain starts at residue 1. TED-redundant used a different function without that shift (details in §3). |
| 30 | 195.2M high + 129.2M medium = 324.4M | Correct | — |

## 2. Corrected script

> For our project we want a model that, from the amino-acid sequence alone, predicts where the domain boundaries
> are and which structural family each domain belongs to. Our training data come from TED, The Encyclopedia of
> Domains.
>
> **Dataset.** TED starts from version 4 of the AlphaFold Protein Structure Database: 214.7 million predicted
> structures, one per UniProt sequence. Many sequences occur more than once, so TED keeps one representative per
> unique sequence. This gives TED-100, 188.9 million structures. The duplicates were not thrown away: identical
> sequences can have noticeably different AlphaFold models (about 42% of identical-sequence groups differ by more
> than 1 Å). So the remaining 25.8 million duplicate models form a separate set, TED-redundant, which was chopped
> the same way.
>
> **Chopping.** Every model is cut into domains by three independent parsers that work on the 3D structure:
> Merizo, Chainsaw and UniDoc. A domain is accepted only if the parsers agree. Domains supported by all three are
> "high consensus", domains supported by two are "medium consensus", and single-method predictions are discarded.
> This gives 324.4 million TED-100 domains (195.2M high + 129.2M medium), plus 40.4 million from TED-redundant:
> 364.8 million in total.
>
> **Classification.** Two things happen in parallel on the TED-100 domains:
> 1. MMseqs2 clusters the domain sequences: at least 50% identity, and the shorter sequence at least 90% of the
>    longer. That gives 120.7 million clusters.
> 2. Every domain structure is compared with CATH.
>
> **CATH** classifies experimentally determined protein domains into a four-level hierarchy:
> - **Class**: secondary-structure content (mainly α, mainly β, α/β, few secondary structures).
> - **Architecture**: how the helices and strands are arranged in 3D, ignoring connectivity; for example bundles,
>   barrels, sandwiches, rolls.
> - **Topology**: the fold, i.e. the same arrangement and the same connectivity of secondary-structure elements.
> - **Homologous superfamily**: domains with enough structural, sequence or functional evidence of common
>   ancestry.
>
> A code like 3.40.50.300 reads α/β · three-layer sandwich · Rossmann fold · P-loop NTPase superfamily. Sharing a
> topology but not a superfamily means CATH has no convincing evidence of homology. The shared fold may be
> convergent, or the proteins may be very distant relatives. TED used CATH 4.3.
>
> The structural comparison has two steps:
> 1. **Foldseek** encodes each structure as a sequence of 3Di structural letters, which makes searching fast. TED
>    searched every domain against 31,574 representative CATH domains. A strong match (E-value < 0.019,
>    TM-score > 0.56) earns the full four-level label (H). A weaker match that still covers most of the domain
>    (TM-score ≥ 0.42, coverage ≥ 0.786) only earns the fold (T). The thresholds were tuned for 98% precision.
> 2. For domains Foldseek could not place, a deep-learning embedding search (**Merizo-search / Foldclass**) finds
>    the nearest CATH domain, and TM-align confirms it. These matches only earn a fold (T) label.
>
> Together the structural steps label 193.9M domains at H and 45.8M at T, and leave 84.7M unlabelled. Labels are
> then spread through the sequence clusters: an unlabelled domain in a cluster that has labelled members inherits
> their label. This rescues another 11.6M domains. In total 251.3M domains (77%) can be placed in CATH, and 73.1M
> domains, in 41.9M clusters, cannot.
>
> Unlabelled does not mean novel. It can also mean a poor model, a non-globular or repetitive structure, or a
> relative too distant to detect. TED therefore kept only compact, well-folded, confident clusters (globularity,
> at least 6 secondary-structure elements, pLDDT80 ≥ 90). It searched those against PDB, CATH, ECOD and SCOPe,
> set aside highly symmetric repeats (6,433 clusters), removed badly chopped domains, and ended with 7,427 clusters
> of putative novel folds.
>
> **The three parsers.**
> - **Merizo** combines the sequence with the 3D backbone (distance map and residue frames) using invariant point
>   attention, and assigns every residue either to a domain or to "non-domain", which matters for the disordered
>   linkers common in AlphaFold models. It does not use pLDDT.
> - **Chainsaw** predicts, for every pair of residues, the probability that they belong to the same domain. It
>   uses the distance map plus secondary structure, then greedily finds the domain assignment that best matches
>   those probabilities. Discontinuous domains come out naturally, and residues can be left unassigned.
> - **UniDoc** is not machine learning. It turns inter-residue distances into contact scores and splits the chain
>   top-down, never inside a helix or strand, then merges fragments bottom-up. Because it cannot recognise
>   non-domain regions, TED runs it only on the residues that Merizo put into domains.
>
> **The consensus.** Each parser's output is first cleaned: segments shorter than 5 residues and domains shorter
> than 25 residues are removed. All predicted domains of a protein are then compared pairwise. Two predictions
> agree when their intersection-over-union is at least 0.7, so their start and end points do not need to be
> identical. Predictions are nodes in a graph, agreeing predictions are connected, and the size of each connected
> group gives the consensus level (3 = high, 2 = medium, 1 = low). Only high and medium domains are kept; they
> never overlap. The final boundaries are the residues shared by all agreeing predictions, and the same size filter
> is applied once more. Residues outside accepted domains stay unassigned. Overall TED-100 has 195.2M high-consensus
> and 129.2M medium-consensus domains, 324.4M in total.

## 3. Things we found by re-running TED's code (relevant for our model)

- **Off-by-one in the published boundaries.** TED-100 consensus ranges are the parsers' agreed ranges shifted by
  +1 residue. This comes from `domstr_to_assignment(..., offset=1)` in `domain_consensus.py`. We confirmed it in the
  released code, via the TED API, and across 2M database rows: 0 TED-100 domains start at residue 1, while
  129,876 start at residue 2. TED-redundant was made with another function
  (`calculate_domain_consensus_redundant`): no +1 shift, but Chainsaw's 0-based output was used without its +1
  correction. That reproduces 40/49 TED-redundant chains, versus 0/49 with the TED-100 path. **For training
  labels, keep the two subsets apart, or shift TED-100 boundaries back by one residue.**
- **UniDoc version.** TED's published UniDoc outputs are reproduced by UniDoc's 2022 static build
  (`UniDoc_structure`), not by the `UniDoc_struct` binary the released wrapper calls. On 106 chains, 97 raw UniDoc
  outputs are identical with the 2022 build versus 38 with `UniDoc_struct`.
- With those two points handled, our re-implementation reproduces the published TED-100 domains exactly for
  **202 of 209 random chains (381/387 domains, 98.4%)**, starting from the same AlphaFold models.
- `md5_domain` in TED's tables is the MD5 of the chain sequence cut with the published chopping string as written.
- Supp. Table S1 has a typo in the high-consensus count (use 195,229,271).
