# Juicebox curation worklist — v3 YaHS scaffolds (2026-10-16)

Files (server `01.Genome_Assembly/new-version/03_scaffold/primary_yahs/`):
```
scp <server>:/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/03_scaffold/primary_yahs/out_JBAT.hic .        # 379 MB
scp <server>:/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/03_scaffold/primary_yahs/out_JBAT.assembly .   # 17 KB
```
Open in Juicebox Assembly Tools 1.11.08 (File → Open .hic; Assembly → Import Map Assembly → .assembly). General instructions: `JUICEBOX_CURATION_GUIDE.md`.

Evidence used below: Hi-C contact maps (PretextMap, `figures/primary_yahs_FullMap.png`, per-scaffold snapshots) and whole-genome alignment to *Gopherus evgoodei* rGopEvg1 (`figures/dotplot_primary_yahs.png`, `figures/scaffold_to_GopEvg.tsv`). rGopEvg1 names 24 chromosomes; Testudinidae are 2n = 52, so 26 is the expected count and two microchromosomes are simply unnamed in the reference. Synteny is a *guide*: a real *M. impressa* rearrangement will disagree with *Gopherus* and still show continuous Hi-C contact. The Hi-C map decides; synteny only tells you where to look.

## A. Quick confirmation (scaffolds 1–11, 81–383 Mb)
Each maps 95–100% to one reference chromosome, collinear, no second hit. Contact maps continuous. Scroll each once at ~1 Mb resolution looking for a contact break at a contig boundary; none expected. Note orientation: s2, s3, s7–s10 align to the reference in reverse — leave as is; orientation is set in the naming step after curation.

## B. Scaffolds that need a decision

| scaffold | size | synteny | what to check in Juicebox | likely action |
|---|---|---|---|---|
| **s12** | 57.2 Mb | 58% chr14 + 32% chr21 — two reference chromosomes | Is there a contact break between the chr14-part and the chr21-part? Zoom to the contig boundary nearest the switch. | Break if contact drops to background → two chromosomes (brings count to 26 named + 1). Keep if contact is continuous → genuine fusion in *M. impressa*; record as such. |
| **s22** | 20.9 Mb | 42% chr3; s3 (217 Mb) already covers 99.9% of chr3 (225 Mb) | Does s22 show strong contact with either end of s3? Could be the missing ~8 Mb of chr3 plus unalignable sequence — or a retained haplotig of a chr3 region. | Join to s3 if end-contact is strong; move to debris if it mirrors a region already in s3 (parallel faint diagonal). |
| **s17** | 27.1 Mb | only 6% aligns (chr9/chr21 fragments) | Internal map: continuous diagonal (a real, repeat-rich chromosome, e.g. rDNA/satellite-bearing microchromosome) or a patchwork? | Keep if internally continuous; it is probably a genuine chromosome the reference lacks. |
| **s24** | 16.4 Mb | <1% aligns | Same test as s17. | Same. |
| **s14** | 38.4 Mb | 78% chr15 (rest unaligned) | Look for a contact break at the unaligned portion. | Likely fine. |
| **s16, s18, s23, s26** | 15–30 Mb | 75–89% to one chromosome | Quick scroll each. | Likely fine. |

## C. Small scaffolds that probably belong to a chromosome (join candidates)

| scaffold | size | aligns to | the chromosome is held by | action |
|---|---|---|---|---|
| s44 | 1.85 Mb | chr7 (99.6% of s44) | s7 | drag to the end of s7 whose contact it shows; if contact is weak on both ends it is a haplotig — debris |
| s52 | 1.51 Mb | chr23 (99.3%) | s25 | same test with s25 |
| s60 | 1.06 Mb | chr11 (80%) | s11 | same test with s11 |
| s35 | 3.66 Mb | chr1 (35%) | s1 | check contact with s1 ends; otherwise leave unplaced |
| s39 | 2.52 Mb | chr19 (21%) | s18 | check; otherwise leave unplaced |
| s58 | 1.17 Mb | chr9 (24%) | s9 | check; otherwise leave unplaced |

A small scaffold whose contact is spread along the *whole* length of a chromosome rather than concentrated at one end is a collapsed-repeat or haplotig fragment: move it to debris, do not join.

## D. Log every edit (this table goes into the manuscript as a supplementary table)

| scaffold | position (Mb) | issue | action taken | evidence |
|---|---|---|---|---|
| | | | | |

## E. When done
Export `Assembly → Export Assembly` as `curated.review.assembly`, copy it back to the server directory, and tell me — I apply it with `juicer post`, re-run QC and the contact map, then name and orient chromosomes (chr1–chrN by size; orientation matched to *G. evgoodei*) and proceed to annotation on the curated sequence.

---
## Resolved by Hi-C read-pair counts (2026-10-18; `scripts/genome_v3/curation_contacts.py`, 20,760,266 pairs; outputs `figures/curation_contacts_*.txt`)

| item | Hi-C evidence | verdict / Juicebox action |
|---|---|---|
| s12 seam (~33 Mb) | across-seam contact 6.6× above inter-chromosomal background; across/within ratio 0.32 at 2 Mb rising to 0.9 at ≥10 Mb (centromere-like boundary) | **keep; no cut** — one chromosome; chr14+chr21 fusion relative to *G. evgoodei* |
| s22 vs s3 | enrichment 0.84× (background), uniform along both | **keep separate**; *G. evgoodei* chr3 = s3 + s22 |
| s17 | only partner s14 at 3.8×, spread over all of both (compartment co-association) | keep; repeat-rich, no action |
| s24 | no partner ≥1.04×; continuous | keep; independent chromosome |
| s52 → s25 | 13×; 978 pairs in s25's first 5 Mb | **join to START of s25** |
| s60 → s11 | 4.6×; 385 at s11's start vs median 59 | **join to START of s11** |
| s44 → s7 | 3.0×; 318 at s7's start vs 83; from s44's first half | **join to START of s7, s44 inverted** |
| s39 → s5 (not s18) | 3.1× with s5, 230 at s5's end vs median 30; s18 flat | **join to END of s5** |
| s35 → s1 | 2.5×; 462/442/320/151 over s1's first 20 Mb (distance decay) | **join to START of s1** (confirm corner) |
| s27 → s14 | 4.8×; 500 at s14's start vs 157; from s27's second half; also 5.1× with s17 | probably **join to START of s14** (confirm corner) |
| s58 | 1,011 pairs only; weak hint of s5 end | leave unplaced unless corner unmistakable |
| s1 42.86–43.06 Mb | 250 kb with zero Hi-C reads inside ptg000003l; unique sequence, normal GC/MboI | no edit; HiFi coverage (job 94661, then 94667 on the release) confirms a heterozygous tandem duplication: 254 reads present, 239 at MAPQ 0; unfiltered depth 18.5× |

Chromosome-scale count after these actions: still 26.

---
**Status 2026-10-19:** all decisions above were applied by `apply_curation.py` (3 breaks, 3 joins), verified by remapping (job 94662: new junctions bg 57.8 / 39.4 / 16.1, no join < 5), and released as `Mimp_v3` (chr1–chr26 + 368 unplaced). This worklist is historical; the release-level evidence is in `PROTOCOL_v3.md` §3f–6 and Data file 3.
