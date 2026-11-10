# Manouria impressa genome v3 — complete protocol log

Every step run for the new assembly/annotation, in order, with the exact command, software version,
parameters, inputs, outputs and the numbers that came out. Written as the work is done, not afterwards.
This file is the source for the Methods, Table 1 (software), Table 2 (statistics) and the Code
Availability statement of the revised Data Descriptor. Mirrored in the analysis repository at
docs/genome_paper/PROTOCOL_v3.md.

Directory: /data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version
Compute: Slurm nodes the 64-CPU, 480 GB nodes (64 CPU, 480 GB each); env /data1/opt/conda/envs/assembly (pinned, see §0.2)
Started: 2026-10-15

## 0. Starting point and why a reassembly

- Previous release v2.4-20251210 (NGDC GWHHIVU00000000.1): 402 sequences (32 chr + 370 unplaced), 2,429,843,097 bp,
  scaffold N50 124.75 Mb, contig N50 86.9 Mb (from AGP; 502 components, 100 gaps of 100 bp, 10,120 N total),
  18,093 genes = 18,093 transcripts (no isoforms), 5,805 transcripts "repaired by ORF prediction" post hoc.
- The v2.4 package contains NO assembly/annotation pipeline scripts, no versions/parameters, no QC reports,
  and its METHODOLOGY.md (SALSA2) contradicts the manuscript (HiC-Pro + LACHESIS) and Table 1 (LACHESIS + YAGCloser).
  Manuscript statistics are internally inconsistent (three different N50 pairs; 240 vs 402 vs 81 sequences).
- prior peer review asked for HapHiC-based scaffolding on hifiasm output, two haplotypes, Juicebox curation.
  prior peer review asked for reconciled statistics, per-library read accounting, reproducible methods, sample provenance.
- Raw data available: PacBio HiFi (Revio, ccs 7.0.0, ≥4 SMRT cells; manuscript wrongly says Sequel II),
  Hi-C (MboI, 2 × ~3 GB gz), RNA-seq heart (IT3) + lung (IT4) PE150. No Illumina WGS reads exist on the server
  despite the manuscript section "Generation of short reads for genome correction".
- Decision (2026-10-15, authors): rebuild from raw reads with a fully scripted, versioned pipeline.

## 0.1 Raw data accounting
(job 94633, scripts/00_raw_stats.sbatch — results appended below when complete)

## 0.2 Software environment

Job 94633 (`scripts/00_raw_stats.sbatch`; samtools 1.21 view → awk over every read; zcat | awk over FASTQ):

| Library | Platform (from read headers) | Reads / pairs | Bases | Mean / N50 read length | Notes |
|---|---|---|---|---|---|
| PacBio HiFi | **Revio** (movies m84069_230831, m84129_231016, m84129_231026, m84129_231103; ccs 7.0.0; binding kit 102-739-100). The submitted manuscript said Sequel II | 3,618,287 reads | 68,548,442,464 (68.55 Gb) | 18,945 / 18,784 bp; max 65,258 | all reads rq ≥ 0.99 (Q20) already; 4 SMRT cells: 47.10 + 7.69 + 9.98 + 3.77 Gb |
| Hi-C | **MGI DNBSEQ** (read names E150015915…), MboI; 55% of reads carry the GATCGATC ligation junction (first 100k reads) | 32,770,909 pairs | 9.83 Gb (2 × 4.92) | 150 bp | manuscript's "724.3 million paired-end reads" is wrong; Table 2's 32.78 M spots / 9.83 Gb is right |
| RNA-seq heart (IT3) | **MGI DNBSEQ** (read names E200012209…; manuscript said NovaSeq — confirm) PE150 | 22,319,011 pairs | 6.70 Gb | 150 | |
| RNA-seq lung (IT4) | MGI DNBSEQ (as IT3) PE150 | 22,448,804 pairs | 6.73 Gb | 150 | |
| Illumina WGS | — | none on server | | | the submitted manuscript's §2.3 "short reads for genome correction" describes data that do not exist here; no short-read polishing was or will be done |

HiFi coverage = 68.55 Gb / 2.43 Gb assembly = **28.2×** (/ 2.13 Gb GenomeScope haploid estimate = 32.2×). The manuscript's "40.9×" is not reproducible from these numbers.
Vendor GenomeScope2 (k=21, from the submission folder): haploid length 2,127–2,133 Mb, heterozygosity 0.31–0.35%, repeats 468 Mb (22%), error 0.18%. Re-estimated in step 1.

## 0.2 Software environments (conda, /data1/opt/conda/envs/)

| env | tools (pinned) |
|---|---|
| assembly | hifiasm 0.25.0-r726, yahs 1.2.2, purge_dups 1.2.6, minimap2 2.28-r1209, samtools 1.21, bwa 0.7.18, chromap 0.2.6-r490, pairtools 1.1.0, seqkit 2.9.0, gfatools 0.5, bedtools 2.31.1, assembly-stats 1.0.1, PretextMap 0.1.9, PretextSnapshot 0.0.4, pigz |
| genomeqc | meryl 1.4.1, merqury 1.3, genomescope2 2.0 |
| busco | busco 5.8.0 (building) |
| compleasm | compleasm 0.2.6 (building) |
| annotation | BRAKER3 via Singularity container (to be pulled); RepeatModeler2/RepeatMasker (Dfam TETools container); eggNOG-mapper 2.1.12 (existing env) |

Rejected: one env containing QUAST (pins blast < 2.16) and BUSCO (dendropy pin) could not be solved — split as above. QUAST dropped; contiguity statistics come from assembly-stats and `scripts/genome_v3/asm_stats.py`.

## 1. HiFi read preparation and k-mer survey
Script `scripts/genome_v3/01_hifi_prep.sbatch` (job 94636). A first attempt (job 94634) wrote an empty FASTQ because `samtools fastq -0 /dev/null` routes *unpaired* reads to `-0`; recorded so nobody repeats it.
- `samtools view -e '[rq]>=0.99'` → `samtools fastq -T rq` (keeps the rq tag) → pigz
- Adapter screen: seqkit grep for the SMRTbell adapter (≤ 2 mismatches), the HiFiAdapterFilt criterion. Result: (appended)
- `meryl k=21 count` → GenomeScope2 `-p 2`: (appended)

## 2. Contig assembly
Script `scripts/genome_v3/02_hifiasm.sbatch` (job 94637, depends on 94636).
`hifiasm 0.25.0 -t 60 --h1 <Hi-C R1> --h2 <Hi-C R2> <HiFi>`, defaults otherwise. Outputs `IT1.hic.p_ctg` (primary), `IT1.hic.hap1.p_ctg`, `IT1.hic.hap2.p_ctg` (phased haplotypes). Statistics: (appended)

### 1 — results (job 94636, finished 2026-10-15 20:36 UTC)
- Reads with rq ≥ 0.99: 3,618,287 / 3,618,287 (all). Adapter screen removed **5 reads** (SMRTbell adapter, ≤ 2 mismatches) → 3,618,282 reads, 68,548,354,278 bp, mean 18,945, N50 18,784, max 65,258 bp.
- GenomeScope2 (k = 21, p = 2): haploid length **2,255.6–2,261.8 Mb**, heterozygosity 0.308–0.345%, repeat 596.8–598.4 Mb (26.4%), unique 1,658.8–1,663.4 Mb, error 0.166%, model fit 74.6–97.8%. HiFi coverage = 68.55 / 2.26 = **30.3×** (vs 28.2× of the final assembly length).
- hifiasm's own k-mer analysis: single homozygous peak at 27×, no heterozygous peak (`peak_hom: 27; peak_het: -1`) — a low-heterozygosity genome.

### 2 — results (job 94637, 4 h 52 min on 60 threads, peak RSS 118.6 GB, finished 2026-10-16 01:28 UTC)

| assembly | total bp | contigs | contig N50 | longest | mean |
|---|---|---|---|---|---|
| `IT1.hic.p_ctg` (primary) | 2,428,913,094 | 460 | **86,871,185** | 303,707,544 | 5.28 Mb |
| `IT1.hic.hap1.p_ctg` | 2,376,624,100 | 686 | 36,595,271 | 106,955,949 | 3.46 Mb |
| `IT1.hic.hap2.p_ctg` | 2,385,044,233 | 631 | 28,880,847 | 117,296,090 | 3.78 Mb |

The primary assembly reproduces the v2.4 contig set to the base pair in N50 (86,871,185 vs 86,873,400 from the v2.4 AGP) and longest contig (303.7 Mb): v2.4 was a hifiasm primary assembly of these reads, and the rebuild is reproducible with documented versions. The improvements of v3 therefore come from scaffolding, curation and annotation, not from contigging.

## 4 — baseline QC of v2.4 (jobs 94638 + 94640)
assembly-stats / asm_stats.py: 402 sequences, 2,429,843,097 bp, scaffold N50 124,747,853, contig N50 86.9 Mb, 112 gaps (100 × 100 bp + 12 × 10 bp), 10,120 N.
compleasm 0.2.10: vertebrata_odb12 S 98.82% D 0.94% F 0.12% M 0.12% (n 3,390); sauropsida_odb12 S 98.69% (n 6,038).
BUSCO 5.8.0 vertebrata_odb10 (miniprot): C 99.3% [S 97.8, D 1.6] F 0.3 M 0.4 (n 3,354) — 157 complete BUSCOs contain internal stop codons (E 4.7%). Resolves R2's "3,354 vs 3,295": 3,354 is the database size.
Merqury 1.3 (k = 21 HiFi db): **QV 67.43** (9,213 erroneous k-mers over 2,429,822,713 bp; error 1.8 × 10⁻⁷), k-mer completeness 95.43%. A first Merqury run (job 94638) read a half-built meryl db and a second (94639) had a wrong $MERQURY path — both discarded.

## 3. Scaffolding (submitted 2026-10-16 02:10 UTC)
- `03_scaffold_yahs.sbatch primary_yahs` (job 94641): chromap 0.2.6 `--preset hic -q 1 --remove-pcr-duplicates` → pairs → YaHS 1.2.2 `-e GATC`; Juicebox files via `juicer pre -a` + juicer_tools 1.22.01; Pretext map.
- `03b_scaffold_haphic.sbatch` (HapHiC commit 74f1cb2, `haphic pipeline <ctg> <bam> 26 --RE GATC`; 26 = Testudinidae 2n = 52): hap1 (94642), hap2 (94643), primary (94644) — prior peer review's request; YaHS vs HapHiC on the primary will be compared and the better-supported one curated.
- QC of primary contigs (94645) and Merqury of hap1 + hap2 together for phasing statistics (94646).

## Decisions (authors)
- **2026-10-16 — deposition:** v3 will be released as a new version of the existing NGDC Genome Warehouse record (GWHHIVU00000000 → .2) under BioProject PRJCA038476, not as a new accession. The chitin-metabolism manuscript's genome statistics will be updated to the v3 numbers once released.
- Expected improvement is in chromosome-scale structure, haplotypes, annotation and documentation — not in contig contiguity or base accuracy, which v2.4 already had (QV 67.4). The revised manuscript states this plainly.

### 3 — attempt 1 (jobs 94641–94644), discarded
- YaHS fed a BED derived from chromap `--pairs` output: links parsed (19,041,557 pairs; 50.3% inter-contig, 18.6% cis < 1 kb, 26.3% cis > 20 kb) but ten scaffolding rounds made **zero joins** (586 sequences = 460 contigs split by error-correction, no U rows in the AGP). Kept in `03_scaffold/primary_yahs_attempt1_bed/` for the record.
- `chromap --SAM -o /dev/stdout` aborts (core dump) on these data; it killed the HapHiC runs and the Pretext step.
- Fix (attempt 2, jobs 94647–94654): `03_hic_map.sbatch` maps Hi-C once per contig set to a SAM on disk → name-sorted BAM (YaHS input) + `-F 3340` filtered BAM (HapHiC input); YaHS run from BAM both with default contig error correction and with `--no-contig-ec`; HapHiC from the filtered BAM. Pretext maps regenerated from the final scaffold BAM in step 4.

### 3 — attempt 2 (jobs 94647–94654): Hi-C mapping works; HapHiC scaffolds the haplotypes; YaHS still reads no pairs
- chromap 0.2.6 to SAM on disk, per assembly: 46.27–46.38 M uniquely mapped reads, **41.56–41.66 M mappings passing filters** (MAPQ ≥ 1, PCR duplicates removed) — ~23 min each (jobs 94647/94648/94649).
- **HapHiC (hap1, hap2; 26 groups, 10 min each):**

| assembly | contigs → scaffolds | contig N50 → scaffold N50 | longest | N90 |
|---|---|---|---|---|
| hap1 | 686 → 462 | 36.6 Mb → **87.8 Mb** | 146.6 Mb | 15.9 Mb |
| hap2 | 631 → 403 | 28.9 Mb → **91.6 Mb** | 141.8 Mb | 18.6 Mb |

  So the Hi-C library carries ample long-range signal; the "50% inter-contig pairs" statistic reflects genuine inter-chromosomal contacts plus contig fragmentation, not noise.
- **YaHS (jobs 94650/94651) read 0 read pairs from the BAM**: chromap names mates `NAME/1` and `NAME/2`; YaHS pairs consecutive records only when QNAMEs are identical. The attempt-1 BED had the same suffix problem in a different form. Attempt 2 kept in `03_scaffold/primary_yahs_attempt2_bam_nopairs/`.
- HapHiC on the primary (94652) failed on a BAM-path lookup bug in my script (empty pattern); fixed by passing the mapping label explicitly.
- Attempt 3 (jobs 94655 YaHS, 94656 HapHiC, both on the primary): mate suffixes stripped from QNAME (`awk sub(/\/[12]$/,"",$1)`) before YaHS.

### 4 — v3 primary contigs QC (job 94645, in progress)
compleasm vertebrata_odb12 S 98.79% D 0.97% M 0.12% (n 3,390); sauropsida_odb12 S 98.66% D 1.05% M 0.20% (n 6,038) — identical to v2.4 within one gene, as expected for the same contigs.

### 3 — attempt 3 results (jobs 94655 YaHS, 94656 HapHiC on the primary; 2026-10-16 07:00 UTC)
YaHS read **20,760,223 read pairs** from 41,557,772 records (10,056,349 intra-contig + 10,703,874 inter-contig links) once mate suffixes were stripped.

| scaffolder | sequences | joins (gaps) | scaffold N50 | N90 | longest | ≥ 14 Mb scaffolds |
|---|---|---|---|---|---|---|
| **YaHS 1.2.2** (default, `-e GATC`) | 394 | 165 | **138,646,552** | 21.7 Mb | 382,984,538 | **26** (26th = 14.7 Mb; 27th = 8.8 Mb — a clean karyotype step, 2n = 52) |
| HapHiC 74f1cb2 (`pipeline … 26 --RE GATC`) | 404 | 56 | 106,232,630 | 13.9 Mb | 303,707,544 | no clean step |
| v2.4 deposited | 402 | 100 | 124,747,853 | — | 303,622,320 | 32 named "chromosomes" |

Note: the submitted manuscript's Table 2 values (scaffold N50 148.3 Mb, longest 384.3 Mb) match a YaHS-like scaffolding of these contigs far better than the deposited v2.4 FASTA does; the deposited version appears to have been re-split into 32 sequences afterwards.
Contact maps of both scaffold sets rendered with PretextMap for inspection (jobs 94657/94658, `03c_contactmap.sbatch`). Choice between YaHS and HapHiC made on the maps, then Juicebox curation.

### 4 — haplotype Merqury (job 94646)
hap1 + hap2 together: **QV 68.54** (13,988 erroneous k-mers over 4.76 Gb; error 1.4 × 10⁻⁷); k-mer completeness hap1 94.81%, hap2 95.00%, **both 99.52%**.

### 3 — decision (2026-10-16 09:00 UTC): curate the YaHS scaffolds
Contact maps (PretextMap, jobs 94657/94658; `figures/primary_*_FullMap.png`, per-scaffold snapshots in `contactmap_primary_yahs/snapshots_each/`): both scaffold sets clean at genome scale; YaHS scaffolds 1–3 (383 / 302 / 217 Mb) show continuous diagonals, two-arm structure in scaffold 2, no contact drop at any join. YaHS composition of scaffold_1: ptg000037l (41.1 Mb) + ptg000003l (303.7) + ptg000021l (26.1) + ptg000010l (6.7) + ptg000065l (5.4). Juicebox files for curation: `03_scaffold/primary_yahs/out_JBAT.hic` (379 MB), `out_JBAT.assembly`, `out_JBAT.liftover.agp`; apply edits with `juicer post` (see JUICEBOX_CURATION_GUIDE.md §5).
### 4 — v3 primary QC complete (job 94645): Merqury QV **67.54**, completeness 95.43%; BUSCO vertebrata_odb10 C 99.3% [S 97.8, D 1.6] — same as v2.4.
### 5 — repeats started (job 94659, `05_repeats.sbatch v3_yahs`): RepeatModeler2 −LTRStruct on the YaHS scaffolds; library will be re-applied to the curated assembly. *G. evgoodei* rGopEvg1 (GCF_007399415.2) downloading to /data1/reference/genomes for synteny and chromosome naming.

### 3e — synteny vs *G. evgoodei* rGopEvg1 (job 94660, `03d_synteny.sbatch`; minimap2 2.28 -x asm20, alignments ≥ 50 kb, MAPQ ≥ 30)
Reference: GCF_007399415.2, 24 named chromosomes. v3 scaffolds 1–11 (383–81 Mb) each map 95–100% to one reference chromosome (chr1–chr11, in order), collinear; s2, s3, s7–s10 in reverse orientation relative to the reference. Attention items handed to curation (`CURATION_WORKLIST_v3.md`): s12 = chr14 (58%) + chr21 (32%) (fusion vs misjoin, decided by Hi-C); s22 42% chr3 while s3 already covers chr3 (join or haplotig); s17 and s24 essentially unaligned (repeat-rich chromosomes or artefacts); join candidates s44→s7, s52→s25, s60→s11, plus s35/s39/s58. Dot plot `figures/dotplot_primary_yahs.png`; table `figures/scaffold_to_GopEvg.tsv`.
**Curation handed to the authors 2026-10-16 12:30 UTC.**

## 3f — quantitative curation (2026-10-18)
Rationale: manual Juicebox inspection was replaced by a scripted, reproducible test applied identically to every join, with the authors' eye-check of the s1 join as positive control. Two statistics from the scaffold-level Hi-C BAM (`scripts/genome_v3/join_qc.py`; 20,760,266 pairs, 100-kb bins, 2-Mb windows either side of each boundary):
- **ratio** = observed pairs across the boundary / expected for contiguous sequence from the genome-wide distance-decay curve measured inside hifiasm contigs;
- **bg** = observed / inter-scaffold background for the same windows.
Calibration on the 145 scorable YaHS joins: bg median 46.1, 5th percentile 11.3, 10th 25.8; the s1 join at 41,057,562 (eye-verified by the authors in Juicebox) scores ratio 0.62, bg 45.0. Decision rules: **break if bg < 5; accept a new join if bg ≥ 11** (5th percentile of YaHS's own joins).
Candidate adjacencies from `curation_contacts.py` (pairs between a small scaffold and every chromosome-scale scaffold, binned along the partner):

| edit | evidence | applied |
|---|---|---|
| s12 seam at ~33 Mb (chr14/chr21 synteny switch) | across-seam contact 6.6× above inter-chromosomal background; across/within 0.9 at ≥ 10 Mb | **no cut** — one chromosome (fusion relative to *G. evgoodei*) |
| s22 vs s3 | enrichment 0.84×; corner bg 2.4–3.9 (background) | **not joined** — separate chromosome |
| s17, s24 | no positional partner; internally continuous | kept |
| s52 → start of s25 | enrichment 13×; corner bg 58.6 | **joined**, orientation − (first half faces s25: 242 vs 123 contacts) |
| s60 → start of s11 | enrichment 4.6×; corner bg 34.9 | **joined**, orientation + |
| s44 → start of s7 | enrichment 3.0×; corner bg 13.6–14.9 (low edge of accepted) | **joined**, orientation +, flagged marginal |
| s39 → end of s5 (synteny had said s18) | enrichment 3.1× with s5, s18 0.0×; corner bg 7.7–8.6 | **not joined** — associated with s5's end but below threshold (sequence probably missing between) |
| s35 → start of s1 | enrichment 2.5×; corner bg 7.4–8.7 | not joined, same reasoning |
| s27 → s14, s58, s39 → s18 | corner bg ≤ 5 | not joined |
| s27 at 5,748,000 (ptg000019c \| ptg000109l) | bg 3.7; ptg000019c is a 6.08-Mb *circular* hifiasm contig (satellite array) that YaHS had cut into four scaffolds | **broken** |
| s59 at 927,537; s119 at 126,000 | bg 1.7 and 0.0 (1 and 0 read pairs) | **broken** |

`apply_curation.py` → `03_scaffold/curated/v3_curated.{fa,agp,names.tsv,edits.log}`: 394 sequences, 2,428,929,594 bp, scaffold N50 138,646,552, 165 gaps, 26 chromosome-scale scaffolds. Verification: Hi-C remapped to the curated FASTA and all joins re-scored (job 94662, `03e_verify_curation.sbatch`) — results appended below.

### HiFi coverage of the YaHS scaffolds (job 94661, `04c_hifi_coverage.sbatch`; minimap2 2.28 map-hifi; 10-kb windows, MAPQ ≥ 1)
Median depth 27.8× (mean 28.0×); 227 of 243,084 windows < 5% of median, 38 of them on chromosome-scale scaffolds; 153 windows > 2.5× median (collapsed repeats / organelle candidates). Only two chromosome-scale runs ≥ 50 kb: scaffold_1:42.87–43.13 Mb and scaffold_28:0.25–0.57 Mb.
**scaffold_1:42.87–43.13 Mb** (also zero in Hi-C at MAPQ ≥ 1; raised by the authors from Juicebox): 254 HiFi reads present, 239 at MAPQ 0. The 260-kb region aligns to hap2 contig h2tg000086l as a single ~130-kb unit twice over (19.59–19.74 Mb) and to hap1 as 260 kb: a **heterozygous tandem duplication** (2 copies on one haplotype, 1 on the other). Unfiltered depth ≈ 18.5× fits the 3-copies-over-2-haplotypes expectation (~21×), not a doubling artefact (~14×). No edit; reported in Technical Validation. scaffold_28:0.25–0.57 Mb shows the same signature (185 reads, 14 at MAPQ ≥ 1).

## 3g — Hi-C library statistics (job 94674, `03f_hic_stats.sbatch`; 2026-10-19)
32,770,909 pairs sequenced; 55,411,316 reads mapped (84.5 %), 46,271,164 uniquely (70.6 %); **20,760,724 valid pairs (63.4 %)** after PCR-duplicate removal, MAPQ ≥ 1 and proper-pair filtering (chromap `--preset hic -q 1 --remove-pcr-duplicates`); on the 460 primary contigs 48.4 % of valid pairs are cis (31.1 % ≥ 1 kb, 26.3 % ≥ 20 kb). Table: `06_release/Mimp_v3.hic_stats.tsv`.

## 7 — telomeres (job 94675, `09_telomere.py`, 50-kb terminal windows)
(TTAGGG)n/(CCCTAA)n arrays ≥ 300 bp at 5 of 52 chromosome ends (chr3, chr5, chr8, chr23 start; chr10 end; longest 1,134 bp), traces (≥ 60 bp) at 20 more, none at 27; no array ≥ 300 bp in any unplaced scaffold. Table: `06_release/Mimp_v3.telomeres.tsv`.

## 8 — software versions
`software_versions.sh` → `06_release/software_versions.tsv` (27 tools; conda package versions + container versions). compleasm: bioconda package 0.2.6, binary reports 0.2.10 — cite as 0.2.10.

## 9 — audit of 2026-10-19
Gaps found and closed/scheduled: haplotype FCS (94683/94684), mito screen of nuclear files (94679), functional annotation + annotation QC (94680), ncRNA (94685; new env `ncrna`: tRNAscan-SE 2.0.12, barrnap 0.9, Infernal 1.1.5), v2.4→v3 correspondence (94676). Full list: `docs/genome_paper_new_submission/AUDIT_2026-10-19.md`.

### 11 — contamination / organelle, corrections (2026-10-19)
- MitoHiFi: `findMitoReference.py` fails on compute nodes (no NCBI access). Reference fetched on the head node: *M. impressa* mitogenome **NC_011815.1** (16,642 bp, circular) → `07_contamination/mito_ref/`; MitoHiFi 3.2.1 run with `-f NC_011815.1.fasta -g NC_011815.1.gb -t 16 -o 2` (job 94686, `07c_mitohifi.sbatch`).
- FCS-GX: database (~470 GB) downloading on the head node; `fcs_gx_chain.sh` verifies every file size against the NCBI listing before writing `gxdb/.complete`; all FCS scripts require that flag. GX then runs on the primary, hap1 and hap2 in sequence (`07d_fcs_gx.sbatch`).
- Job 94659 (RepeatMasker of the pre-curation YaHS scaffolds) cancelled after 3 d 0.5 h: the library (2,492 families) was complete; the release is masked by job 94671 instead.
- Decision (authors, 2026-10-19): Illumina WGS reads exist but are not used, not deposited; no short-read polishing (not standard for HiFi; QV 67.5).

## 6 — release QC of Mimp_v3 (jobs 94663 QC, 94664 Merqury, 94665 contact map, 94666 synteny, 94667 HiFi coverage; 2026-10-19)
- Merqury (k = 21, HiFi reads): **QV 67.54** (8,979 erroneous k-mers / 2,428,901,914 bp; error 1.76 × 10⁻⁷), k-mer completeness **95.43 %**. Per-chromosome QV 60.8 (chr25) – 88.6 (chr15); most unplaced scaffolds have zero erroneous k-mers.
- compleasm 0.2.10 vertebrata_odb12 (n 3,390): S 98.82 % (3,350) D 0.94 % F 0.12 % M 0.12 %. BUSCO 5.8.0 vertebrata_odb10 (n 3,354): C 99.3 % [S 97.8 %, D 1.5 %] F 0.3 % M 0.4 %. Identical to the contig-level and v2.4 values, as expected for the same sequence.
- HiFi coverage (minimap2 map-hifi, 3,618,282 reads, 219,604 supplementary): 243,083 10-kb windows, median 27.8×, mean 28.0×; 222 windows < 5 % of median (34 on chromosomes); 152 windows > 2.5× median; 153 zero-coverage runs totalling 1,605,506 bp, largest chr1:42.87–42.99 and 43.01–43.13 Mb (the known heterozygous tandem duplication) and scaffold_1:6.20–6.52 Mb.
- Contact map `03_scaffold/contactmap_Mimp_v3/snapshots/Mimp_v3_FullMap.png`; synteny `03_scaffold/synteny_Mimp_v3/dotplot_Mimp_v3.png` (chr1–11 one-to-one with GopEvg chr1–11; chr12 = GopEvg chr14 + chr21; GopEvg chr3 = chr3 + chr22).
- Shortest release sequence is 1,000 bp (scaffold_368) — above the 200-bp INSDC/GWH minimum; nothing to drop.

## 15 — deposit files (2026-10-19, jobs 94673 haplotype release, 94676 v2.4 correspondence, 94686 MitoHiFi, 94687 mito screen)
- **Mimp_v3_hap1.fa**: 462 scaffolds (≥ 1 kb), 2,376,646,500 bp, scaffold N50 87,800,537, longest 146,557,687, 686 contigs (contig N50 36.6 Mb); 154 scaffolds align > 50 % to one primary chromosome (`Mimp_v3_hap1_to_primary.tsv`). compleasm vertebrata_odb12 S 98.67 % D 0.68 %; BUSCO odb10 C 98.8 % [S 97.5, D 1.3] F 0.4 M 0.9.
- **Mimp_v3_hap2.fa**: 403 scaffolds, 2,385,067,033 bp, N50 91,606,729, longest 141,836,076, 631 contigs; 156 assigned. compleasm S 99.03 % D 0.74 %; BUSCO running.
- md5: `06_release/release_md5.txt` (Mimp_v3.fa 2ea66ec4…, .agp e2ccbcb8…, hap1 cecdb401…, hap2 13077897…).
- **Mitogenome** (MitoHiFi 3.2.1, `-f/-g NC_011815.1 -t 16 -o 2`, HiFi reads): contig atg000011l, circular, rotated to tRNA-Phe, **18,133 bp**; 13 PCGs, 22 tRNAs, 2 rRNAs (`shared_genes.tsv`: every reference gene present). NC_011815.1 aligns in two blocks (1–4,870 and 5,123–15,530 of the reference); the extra 1.5 kb is a tandem array of four tRNA-Thr/tRNA-Pro units (15,357–16,326) between CYTB and the control region. Read support checked by job 94692 (`07e_mito_verify.sbatch`: depth along a doubled molecule; reads spanning 15,300 → past the origin). Release copies `06_release/Mimp_v3_mito.{fa,gb,annotation.png}`. The job exit status was "FAILED" only because of a `grep` with no match after the run; outputs are complete (script fixed).
- **Mito screen** (minimap2 asm10, MAPQ ≥ 30): no sequence of Mimp_v3, hap1 or hap2 has any mitogenome alignment ≥ MAPQ 30 → nothing to remove (`*.mito_hits.tsv`).
- **v2.4 → v3 correspondence** (minimap2 asm5, alignments ≥ 50 kb, MAPQ ≥ 30; `v2.4_to_v3_correspondence.tsv`): all 400 v2.4 sequences placed (≥ 90 % of length aligned for every v2.4 chromosome). The 32 v2.4 "chromosomes" are 26 chromosomes plus six fragments: v2.4 chr1 + chr16 + chr17 = v3 chr1; chr5 + chr12 + chr27 = chr3; chr7 + chr28 = chr6; chr9 + chr32 = chr9; v2.4 chr15 (86 % v3 chr14 + 13 % scaffold_3), chr26 (94 % chr24 + 7.5 % scaffold_48) and chr31 (91 % chr25 + 3 % scaffold_76) carry sequence that v3 leaves unplaced (haplotig/duplication decisions). 56 of the 368 v2.4 unplaced sequences lie inside v3 chromosomes; 43 v2.4 sequences are split (> 5 %) across two v3 sequences.
- ncRNA job 94685 failed in 12 s (`tRNAscan-SE --version` is not a valid flag; exit 1 under `set -e`); resubmitted as 94693 after the fix.

### 11e — mitogenome read support (job 94692, `07e_mito_verify.sbatch`)
HiFi reads aligned to a doubled mitogenome (minimap2 map-hifi, primary alignments, MAPQ ≥ 20): 2,074 reads; depth across the tRNA-Thr/Pro array + control region 1,766×. 1,526 reads span from 15,300 bp through the array and past the origin. Net indel length within the array per spanning read (`mito_array_copynumber_hist.tsv`): 0 bp in 677 reads (four units, as the consensus), +300 bp in 789 reads (five units), other values in < 4 %. → **length heteroplasmy of the repeat array (4 vs 5 units)**; deposit the four-unit consensus and state the heteroplasmy in Data validation. Chromosome fraction of the release: 2,284,523,576 bp in chr1–26 = 94.1 % of 2,428,929,594.
### 11a — FCS-adaptor re-run (jobs 94697–94699, `07a_fcs_adaptor.sbatch`)
Earlier adaptor runs (94668, 94683, 94684) produced no report: `av_screen_x -O` (should be `-o`); the error was in the log and masked by `&& echo`. Re-run on primary, hap1, hap2.
### 16c — haplotype divergence (job 94694, `11_hap_divergence.sbatch`), 12b — repeat landscape (job 94696 after masking), 15b — EBP notation 7.8.Q67

## 2026-10-20 00:45 UTC — results landed overnight
- **Haplotype divergence** (job 94694, minimap2 asm5 + paftools call, blocks ≥ 10 kb): hap1 vs hap2 — 2,214,797,928 bp aligned, **4,572,299 SNVs (2.06 / kb)**, 573,150 small indels (0.26 / kb), 50,175 SVs ≥ 50 bp. hap1 vs primary: 2,299,563,592 bp, 2,340,066 SNVs (1.02 / kb), 293,897 indels, 26,383 SVs; hap2 vs primary: 2,317,853,036 bp, 2,209,949 SNVs (0.95 / kb), 281,524 indels, 24,094 SVs. SNV + indel divergence between haplotypes = 0.23 %, consistent with the k-mer heterozygosity estimate (0.31–0.35 %, which also counts SV-borne k-mers); the primary is roughly equidistant from both haplotypes, as a mosaic should be.
- **hap2 BUSCO** (94670): C 99.0 % [S 97.9, D 1.2] F 0.4 M 0.5 (n 3,354).
- **ncRNA** (94693, tRNAscan-SE 2.0.12 `-E`; barrnap 0.9 euk): 24,584 raw tRNA predictions (24,182 non-pseudo with isotype) — the raw set includes tRNA-derived SINE copies; the EukHighConfidenceFilter set is produced by job 94700 and is the number for Table 3. rRNA (barrnap, incl. partial copies): 28S 127, 18S 106, 5.8S 79, 5S 193.
- **FCS-adaptor** (94697–94699): **0 flagged sequences** in Mimp_v3, hap1 and hap2.
- **FCS-GX database** complete (465 GB; two small files re-fetched after a failed transfer; every file size verified against NCBI); `.complete` written 00:46 UTC; FCS-GX running: primary 94701 → hap1 94702 → hap2 94703.
- RepeatMasker of the release (94671) still running at 14 h 46 min (second pass, Dfam Tetrapoda) on 16 cores; BRAKER3 (94672), repeat landscape (94696) and functional annotation (94680) wait on it.

## 12 — repeat masking of the release (job 94671, 15 h 06 min on 16 cores; landscape job 94696→rerun)
RepeatMasker 4.1.7-p1 (rmblastn 2.14.1), species library (2,492 families): **48.34 % masked** (1,174,255,184 bp); interspersed 47.05 %: retroelements 24.05 % (SINE 2.56, LINE 11.49 — CR1/L2 9.28, LTR 10.00 — Gypsy 6.15, ERV 2.31), DNA transposons 7.61 % (hAT 3.41, Harbinger 3.03), rolling-circle 0.12 %, unclassified 15.39 %; small RNA 1.08 %, satellites 0.05 %, simple repeats 0.52 %, low complexity 0.07 %. Second pass with Dfam 3.8 Tetrapoda on the masked output: +0.84 % (20,524,009 bp). Total ≈ 49.2 %. Landscape (`12_repeat_landscape.py`, divergence from the `.out` files, class totals scaled to the `.tbl`): LTR burst at 2–5 % divergence; broad old peak of unclassified/DNA at 15–30 %. Files: `06_release/repeat_landscape_Mimp_v3.{png,tsv}`, `repeat_summary.txt`; `05_annotation/Mimp_v3/01_repeats/{genome.softmasked.fa,repeats.gff,rm_custom/genome.fa.tbl}`.
BRAKER3 (94672) started 01:07 UTC: HISAT2 index built, IT3 alignment in progress.
tRNAscan-SE rerun with `--detail` (94704 → filter 94705): the EukHighConfidenceFilter needs the isotype-specific scores that only `--detail` writes.
## Manuscript figures/tables (2026-10-20)
`scripts/genome_v3/collect_manuscript_inputs.sh` → `09_manuscript_inputs/`; `scripts/genome_v3/make_manuscript_figs_tables.py` → `docs/genome_paper_new_submission/{figures,supplementary,manuscript/tables}`. Status list: `docs/genome_paper_new_submission/FIGURES_TABLES_STATUS.md`.

## 13 — BRAKER3 (2026-10-20)
- Attempt 1 (job 94672, 01:07–02:26 UTC) **failed** after HISAT2: `--busco_lineage=vertebrata_odb10` makes BRAKER run compleasm, which downloads its lineage from the internet; compute nodes have none (`compleasm_to_hints.py: return code 1`). HISAT2 alignments succeeded and are reused: **heart 97.11 %, lung 97.72 % overall alignment rate**. Attempt kept as `02_braker3_attempt1_failed_compleasm/`.
- Attempt 2 (job 94708, started 02:47 UTC): same command without `--busco_lineage` (standard BRAKER3: GeneMark-ETP + AUGUSTUS + TSEBRA with RNA-seq and OrthoDB v12 Vertebrata evidence). Gene-set completeness is measured afterwards by `06b_functional.sbatch` with the offline BUSCO lineage (job 94709, chained).
- Expected: GeneMark-ETP 10–20 h, AUGUSTUS training + prediction 6–12 h on 48 threads → gene models around 21 Oct midday UTC; functional annotation + ncRNA ~3 h later.

## 11f — FCS-GX results and decontamination (2026-10-20, jobs 94701–94703, 94710)
FCS-GX 0.5.5 (gxdb 2023-01-24, 709 Gbp) on the pre-screen releases: primary **33 unplaced scaffolds EXCLUDE, 11,526,814 bp** (12 *Parabacteroides gordonii*, 6 *Bacteroides neonati*, 4 *Lactococcus lactis*, 2 *Erysipelothrix aquatica*, 1 *Prevotella timonensis*, 1 phage, …; largest 4.36 Mb); hap1 48 EXCLUDE (11,545,238 bp) + 1 FIX (2,878 bp inside h1_scaffold_14, hard-masked); hap2 7 EXCLUDE (4,556,263 bp). No chromosome affected. Taxa are gut/post-mortem flora, consistent with liver from an animal that died naturally; they explain the 152 high-depth windows of §6. Applied with `gx clean-genome` (`13_fcs_clean.sbatch`): **Mimp_v3 361 sequences, 2,417,402,780 bp** (26 chr + 335 unplaced renumbered by length; contig N50 unchanged 86,871,185; 509 contigs); hap1 414 / 2,365,101,262; hap2 396 / 2,380,510,770. Pre-screen files kept as `*.pre-fcs.*`; contaminants in `07_contamination/<label>/gx_pre-fcs/contam.fa`; actions table `06_release/fcs_gx_actions_applied.tsv`. Re-screen of the cleaned files: jobs 94713–94715.
### QC of the cleaned primary (94716, 94717)
Merqury **QV 67.93** (8,173 erroneous k-mers / 2,417,377,800 bp; up from 67.54 — the bacterial scaffolds carried most of the "errors"), completeness **95.09 %** (down from 95.43 %: bacterial read k-mers are still in the read set). compleasm vertebrata_odb12 S 98.82 D 0.94 F 0.12 M 0.12; BUSCO odb10 C 99.3 [S 97.8 D 1.5] F 0.3 M 0.4 — unchanged. Haplotype divergence re-run (94720): hap1 vs hap2 4,572,296 SNVs (2.07 / kb) over 2,210,842,385 bp — unchanged.
## 12c — soft-masking corrected (job 94721)
`genome.softmasked.fa` from 05b carried only the Dfam-pass masking (0.85 %): RepeatMasker uppercases its input, so the second pass discarded the first. Correct union built with `bedtools maskfasta -soft` (species-library masked FASTA + Dfam-pass intervals), restricted to the cleaned release and renamed: **`Mimp_v3.softmasked.fa`, 361 sequences, 48.60 % soft-masked**; `Mimp_v3.repeats.gff` 5,637,211 features. BRAKER3 attempt 2 (94708, 6 h 45 min on the wrongly masked genome) cancelled; **attempt 3 = job 94722** on the cleaned, correctly masked genome (started 15:42 UTC after a 6-h queue wait caused by other users' jobs); functional annotation 94723 chained. 05b script fixed.

## 2026-10-21 01:30 UTC — disk-full incident and restart
`/data2` (6.3 TB NFS) reached 100 % at ~16:40 UTC on 20 Oct: the tRNAscan-SE re-run (94711) had found the previous run's output files, asked interactively "Overwrite/Append/Quit?", and — with no terminal — looped writing the prompt into `trnascan.log` until it reached **453 GB**. Casualties: BRAKER3 attempt 3 (94722, failed writing IT4.bam at 59 min), hap1/hap2 QC at the Merqury stage, FCS-GX re-screens (94713–94715), haplotype contact maps (chromap segfault on write). The compleasm/BUSCO parts of the hap1 QC completed before the failure (S 98.67 %, C 98.7 %). Recovery: log deleted; superseded intermediates removed (hifiasm `.tlb.bin` 65 GB, pre-curation HiFi BAM, pre-FCS and old Merqury/meryl databases, attempt-2 BRAKER working dir) → 711 GB free; `06c_ncrna.sbatch` now deletes old outputs first; `disk_watch.sh` logs free space every 10 min and alerts below 150 GB. Restarted 01:33 UTC: **BRAKER3 attempt 4 = 94757** → functional 94758; tRNAscan `-E -H --detail` 94759 → filter 94760; FCS-GX re-screen 94761–94763; haplotype Merqury 94764, hap2 QC 94765, haplotype contact maps 94766–94767.

## 2026-10-21 08:35 UTC
- **FCS-GX re-screen of the cleaned assemblies: 0 flagged sequences** in Mimp_v3 (94761), hap1 (94762), hap2 (94763). All three are deposit-ready on the contamination criterion.
- Cleaned primary AGP: **148 gaps** (17 of the 165 were inside removed bacterial scaffolds); 509 contigs.
- **Haplotype Merqury on the cleaned files** (94764): hap1 QV 69.12 (6,085 erroneous k-mers), hap2 68.51 (7,047), both 68.80; completeness 94.47 / 94.78 / **99.18 % combined**. hap2 compleasm S 99.00 D 0.74; BUSCO C 99.0 [S 97.8 D 1.2] (94765).
- **tRNA** (94759 `-E -H --detail`, 4 h; filter 94760): 24,407 predictions → **3,309 high-confidence** (10,997 pseudogenes; 10,101 below cutoffs); isotypes: Leu 569, Ser 519, Lys 348, Met 248, Val 204, Asn 192, Glu 190, Gly 174, Phe 165, Ala 144, Trp 107, Gln 94, Thr 83, Asp 80, Tyr 78, Ile 77, Pro 77, Cys 67, His 40, Arg 27, SeC 5 (per `trna_hc.log`). rRNA on the cleaned release (barrnap): 18S 92, 5.8S 79, 28S 99, 5S 181.
- Haplotype Hi-C maps (94766/94767) rendered: clean diagonals → Supplementary Fig S4.
- BRAKER3 attempt 4 (94757): GeneMark-ETP running since 03:08 UTC (7 h at this check).
- Table 2 regenerated from the cleaned files (`make_manuscript_figs_tables.py`).

## 2026-10-22 — repeat-fraction wording unified
"49.2 %" (48.34 % species-library pass + 0.84 % Dfam pass) double-counted the overlap between the passes; the correct, non-overlapping value is the union soft-masked over the cleaned release: **1,174,831,695 / 2,417,402,780 bp = 48.60 %**, used everywhere from now on. GenomeScope "Genome Repeat Length" 596.8–598.4 Mb (26.4 % of the 2,255.6–2,261.8 Mb haploid estimate; model fit 74.6–97.8 %, error 0.166 %) = genome fraction in 21-mers with copy number > 1 — a k-mer quantity, not comparable with homology-based masking (diverged copies share no identical 21-mers). The manuscript states both and explains the difference (Data validation, Repeats).

## 15d — structural differences vs *G. evgoodei* with Hi-C/contig evidence (job 94769, `16_synteny_rearrangements.py`)
Blocks ≥ 1 Mb: 22/26 chromosomes single-reference; chr17, chr25 no block (repeat-rich). Candidates and evidence on the Mimp_v3 side (all inside single HiFi contigs): chr12 GopEvg chr14|chr21 junction at 33.43 Mb in ptg000055l (Hi-C across 0.24 of contiguous expectation, 18× background — centromere-like, as in curation); chr2 165.36–166.82 Mb 1.5-Mb GopEvg-chr19 segment in ptg000008l (0.74 / 0.33; 53× / 24×); chr11 strand switch 7.44 Mb and 8-Mb order jump 10.33 Mb in ptg000002l (0.75, 0.50; 54×, 36×); chr12 order jump 40.42 Mb in ptg000063l (0.66; 47×); chr19 strand switch 21.92 Mb in ptg000121l (0.39; 28×). chr3 vs chr22 (GopEvg chr3 split): contact 0.84× background (curation_contacts, 2026-10-18) → separate molecules. No candidate at a contig junction. Supplementary Table S12. Direction of change not inferred.

## 14 — final gene set (2026-10-23/24; jobs 94773 decision, 94778 final, 94780 post)
BUSCO-protein (vertebrata_odb10) of the candidate sets: TSEBRA braker.gtf 17,392 genes C 89.6 [S 87.6 D 2.0]; GeneMark-ETP 27,692 C 90.2 [D 8.5]; AUGUSTUS-with-hints 26,279 **C 97.6 [S 95.8 D 1.9]**; braker + non-overlapping GeneMark 94.4; braker + non-overlapping AUGUSTUS 26,687 C 96.2; TSEBRA keep-all-AUGUSTUS 34,567 C 97.8 [D 19.1]. Pre-registered rule (max C, D < 4 %) → AUGUSTUS set. Evidence filter (≥ 1 RNA-seq read over the gene, eggNOG ortholog, or BUSCO match; eggNOG 22,136, RNA 22,490, BUSCO 3,399): **24,875 genes kept, 1,404 dropped**; 29,927 transcripts; IDs MIMP3_gNNNNN in chromosome order; GFF3 gene→mRNA→exon/CDS, 0 out-of-bounds, 0 duplicate IDs; proteins cleaned of gffread's "." and "*". Final: BUSCO-protein C 97.6 [S 95.8 D 1.9] F 1.4 M 1.0; eggNOG 22,136 (89.0 %; Pfam 20,957, GO 15,241, KEGG 15,259); RNA-seq ≥ 1 read 90.4 %, ≥ 10 reads 75.9 %; 23,000 genes on chromosomes; mean gene 33,694 bp, CDS 1,668, 9.82 exons, intron 4,418, median protein 380 aa, 4,355 single-exon transcripts. Release: `06_release/Mimp_v3.genes.gff3`, `.proteins.faa`, `.cds.fna`, `.eggnog.annotations.tsv`, `Mimp_v3.tRNA_hc.gff3` (3,309), `Mimp_v3.rRNA.gff3` (451), `gene_id_map.tsv`; md5 in `release_md5.txt`.
**GH18 locus on Mimp_v3**: chr4:139,053,434–139,249,474 (minus strand; = version-.1 chr3 12.42–12.63 Mb): MIMP3_g08468/g08469/g08470 (+ strand; = v2.4 088402, which was a fused model), g08472 (= L7), g08473 (L6), g08474 (L5), g08475 (L4), g08476 (L3, GT31-type), g08478 (= 088393 CHI3L-clade), g08479 (= 088392 CHIT1-clade); g08471, g08477 are additional non-GH18 models in the interval. All eight GH18 loci of the chitin-metabolism paper therefore have gene models in Mimp_v3 (`07_post_geneset/gh18_*.tsv`).
Incidents: 94774 failed (AUGUSTUS GTF gene lines carry a bare id); 94776 failed (gffread GFF3 lacked gene/mRNA); 94779 timed out (diamond rejected "." in proteins; error was hidden). All fixed in the scripts.

## 17. Mitogenome GFF3 (2026-11-08)
GWH requires a GFF3 annotation file per assembly; MitoHiFi writes only GenBank. `scripts/genome_v3/17_mito_gb2gff3.py` converts `06_release/Mimp_v3_mito.gb` → `06_release/Mimp_v3_mito.gff3` (seqid Mimp_v3_mito from the release FASTA; 43 genes = 13 CDS with transl_table=2, 2 rRNA, 28 tRNA; repeated tRNA-Thr2/Pro2 copies of the tandem array suffixed _2/_3 and annotated with a Note; control region 16,327–18,133 added as D_loop, inferred from position). Checks: unique IDs, all Parents resolve, every CDS length divisible by 3. md5 d78886f79baa5a320c7ffcc2802b6298 appended to `release_md5.txt`; file added to Data file 7 and to the GWH upload script (which now also uploads the final `Mimp_v3.genes.gff3` + proteins instead of the superseded 05_annotation path).

## 18. GWH validation of Batch0094178 and the annotation fix (2026-11-08/09)
GWH's automatic check (table2asn-style) rejected the first upload (`ngdc_forms/Batch0094178/Mimp_v3.fa.err`): (i) chromosome-assignment CSV used `no/no` — GWH wants `TRUE/FALSE` (fixed: chromosomes FALSE/FALSE, mitogenome TRUE/TRUE); (ii) 25,572 "duplicate feature IDs" — the multi-segment CDS features shared one ID (`tx.cds`), legal GFF3 but refused; (iii) 77 transcripts with CDS length not a multiple of 3; (iv) 33 StartCodon, 58 NoStop, 53 SuspiciousFrame (first CDS phase > 0), 33 BadProteinStart — all AUGUSTUS **partial models** at scaffold ends or gaps, which the validator treats as complete because the GFF3 carried no partial flags. Warnings only: 448 ShortExon, 11 FeatureCrossesGap, 5 non-consensus splice donors.
Fix (`scripts/genome_v3/06i_gff3_gwh_fix.py`, job 94783, 33 s with the genome in memory): CDS segments renamed `tx.cds1…n`; 139 transcripts (0.46 %) flagged `partial=true` with `start_range=.,X` / `end_range=X,.` (NCBI convention) on mRNA and CDS; the untranslated dangling bases of partial CDS trimmed — 52 at the 5\' end (phase → 0) and 27 at the 3\' end — so every CDS is a multiple of 3 with a consistent phase chain. **No protein sequence changed** (AUGUSTUS never translated those bases; `Mimp_v3.proteins.faa` md5 unchanged); `Mimp_v3.cds.fna` trimmed for the 28 affected transcripts to match. Gene/transcript counts unchanged (24,875 / 29,927). Validation: 0 duplicate IDs, 0 orphan parents, 0 CDS not %3, 0 phase-chain errors, CDS FASTA = spliced GFF3 CDS for all 29,927. Pre-fix files kept as `Mimp_v3.genes.pre-gwh.gff3`, `Mimp_v3.cds.pre-gwh.fna`; list of flagged models `Mimp_v3.genes.gwh_partial_models.tsv`. Re-upload: `Mimp_v3.genes.gff3` + both CSVs, then *FTP finish upload*.
**Round 2 (Batch0094178-batch2, 8 Nov 22:36):** duplicate-ID, CDS%3 and SuspiciousFrame errors gone; StartCodon 77 / NoStop 57 / BadProteinStart 77 remain — the validator ignores `partial=true`/`start_range`/`end_range`, i.e. GWH does not accept partial CDS at all (its GFF3 page defers to the NCBI GFF3 spec but the checker does not implement partiality). Decision: deposit a GFF3 **without the 139 partial transcripts** (`Mimp_v3.genes.gwh_deposit.gff3`: 24,750 genes / 29,788 transcripts; 125 genes lose their only transcript, 121 of them on unplaced scaffolds; none of the 8 GH18 loci affected) and archive the **complete** annotation (`Mimp_v3.genes.gff3`, flagged partials) in Data file 11 / Zenodo. The paper\'s statistics stay on the complete set; the Methods state the difference. If a third QC round still fails on the remaining warnings (ShortExon ×448, FeatureCrossesGap ×11, splice donors ×5) they are warnings and can be acknowledged in the submission form.

## 19. Mitogenome annotation curated for GenBase (2026-11-10)
GenBase rejected submission 0010787 (COX1 flagged 5' partial on a complete sequence; tRNA/rRNA without gene features; plus two form fields). Checking every CDS against the sequence and against the RefSeq record of the same species (NC_011815.1) showed that the MitoHiFi/MitoFinder calls needed four corrections: **COX1** starts at the in-frame GTG at 5404 (not 5422; GTG initiator, 515 aa as in the reference); **ND4** starts at the ATG at 10255 (not 10615; 459 aa, overlapping ND4L by 7 nt) and ends with an incomplete T stop at 11632 instead of running into tRNA-His; **ND3** carries the turtle/bird +1 frameshift at gene position 174 → join(9545..9718,9720..9894), 116 aa, incomplete T stop; **COX3** ends at 9476 with an incomplete T stop; ND1 is written with the incomplete TA stop at 3759 (reference convention). All 13 proteins now translate without internal stops at canonical lengths (ND1 323, ND2 346, COX1 515, COX2 228, ATP8 55, ATP6 227, COX3 261, ND3 116, ND4L 98, ND4 459, ND5 593, ND6 174, CYTB 378). tRNA/rRNA got gene features (trnF … trnP; 12S rRNA / 16S rRNA; the array units annotated individually); control region as D-loop. Files: `06_release/Mimp_v3_mito.gff3` (curated; the MitoHiFi call kept as `Mimp_v3_mito.mitohifi_raw.gff3`), `Mimp_v3_mito.tbl` (GenBase 5-column table), `Mimp_v3_mito.proteins.faa`; generator `scripts/genome_v3/17b_mito_annotation_curated.py`. The nuclear GFF3 needed no analogous curation beyond §18 (its partial models are genuine scaffold-end truncations).
