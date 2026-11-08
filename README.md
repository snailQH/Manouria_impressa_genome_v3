# *Manouria impressa* genome assembly Mimp_v3 — pipeline

The scripts, environments and protocol log behind the chromosome-level, haplotype-resolved genome assembly of the impressed tortoise (*Manouria impressa*, Testudinidae): PacBio HiFi + Hi-C of one adult female (IT-1), RNA-seq-supported gene annotation, two Hi-C-phased haplotype assemblies and the mitochondrial genome. Every statistic in the paper is computed by a script in this repository from the released files.

**How to cite the data:** the Zenodo archive of this repository (DOI pending; one DOI for Data files 1–13 and the pipeline).

**Paper:** Wang J, Li Q. A chromosome-level, haplotype-resolved genome assembly of the impressed tortoise, *Manouria impressa*. bioRxiv 2026, doi:10.1101/2025.03.10.642296 (version 2); Data Note submitted to *BMC Genomic Data*.

## The assembly in numbers
| | Mimp_v3 (primary) | hap1 | hap2 |
|---|---|---|---|
| Sequences / total length | 361 / 2,417,402,780 bp | 414 / 2,365,101,262 bp | 396 / 2,380,510,770 bp |
| Chromosomes | 26 (94.5 % of the assembly) | scaffold level | scaffold level |
| Contig N50 / scaffold N50 | 86.9 Mb / 138.6 Mb | — / 87.8 Mb | — / 91.6 Mb |
| Merqury QV / k-mer completeness | 67.9 / 95.1 % | 69.1 | 68.5 |
| compleasm (sauropsida) / BUSCO | S 98.8 % / C 99.3 % | combined 99.2 % | |
| Protein-coding genes | 24,875 (BUSCO-protein C 97.6 %) | | |
| Repeats | 48.6 % | | |
| Mitogenome | 18,133 bp, circular | | |

## Data
| data | where |
|---|---|
| PacBio HiFi reads (liver) | NGDC GSA CRA024571, run CRR1746698 (BioProject PRJCA038476, BioSample SAMC4973488) |
| Hi-C reads (liver, MboI) | GSA CRA024571, run CRR1746699 |
| RNA-seq heart (IT3) / lung (IT4) | GSA CRA024586, runs CRR1747397 / CRR1747398 (BioSamples SAMC4973532 / SAMC4976828) |
| Reference assembly + annotation | NGDC GWH GWHHIVU00000000.2 (version .1 = preliminary provider assembly) |
| hap1, hap2, mitogenome | GWH, same BioProject (accessions pending) |
| Data files 1–13 of the Data Note (extended description PDF, tables, curation evidence, QC, haplotype/mitogenome/synteny/repeat/annotation results, this pipeline) | `data_files/` here; archived on Zenodo with each release (DOI pending) |

## Layout
```
scripts/              the pipeline, numbered by step (Slurm job scripts + Python helpers):
                        00 raw-read statistics · 01 HiFi preparation · 02 hifiasm (Hi-C mode) · 03 Hi-C mapping, YaHS / HapHiC scaffolding,
                        contact maps, synteny, curation verification, Hi-C statistics · 04 QC (compleasm, BUSCO, Merqury, coverage) ·
                        05 repeats (RepeatModeler2 / RepeatMasker, soft-mask union) · 06 BRAKER3, functional + ncRNA annotation, gene-set
                        decision, final gene set · 07 FCS-adaptor / FCS-GX, MitoHiFi · 08 haplotype + mitogenome release · 09 telomeres ·
                        10 comparison with the deposited version .1 · 11 haplotype divergence · 12 repeat landscape · 13 FCS clean ·
                        14 Pretext renders · 16 synteny breakpoint evidence; helpers: apply_curation.py, curation_contacts.py, join_qc.py,
                        name_chromosomes.py, asm_stats.py; chains/watchers; software_versions.sh; collect_manuscript_inputs.sh
scripts/environments/ the eight conda environments (conda env export --no-builds) + container list
protocol/             PROTOCOL_v3.md — the chronological run log incl. incidents and decisions; genome_assembly_workflow.html (reusable
                      protocol with pitfalls); curation worklist and Juicebox guide
results/              the small output files every table and figure is built from (assembly stats, Merqury, compleasm/BUSCO, Hi-C
                      statistics, curation evidence, telomeres, FCS reports, gene-set statistics, software_versions.tsv)
supplementary/        Supplementary Tables S1–S12 of the Data Note (TSV)
data_files/           Data files 01–13 of the Data Note, one folder each with README.md (what, how produced, columns) and md5.txt; MANIFEST.md
```

## Running it
The protocol log is reproduced verbatim from the project; it occasionally refers to internal working documents (audit notes, author inputs) and to our cluster paths (`/data2/projects/…`), which are not part of this repository — every file needed to rerun or check the analysis is here. Set the project root at the top of each script for another system.
Slurm cluster; conda environments from `scripts/environments/`; containers braker3.sif (BRAKER 3.0.8), tetools.sif (RepeatModeler 2.0.5 / RepeatMasker 4.1.7-p1, Dfam 3.8), fcs-adaptor.sif (0.5.4), fcs-gx.sif (0.5.5), mitohifi.sif (3.2.1). Run the `.sbatch` scripts in numeric order; each states its inputs and outputs at the top and reads paths from one project root. The quantitative Hi-C curation (`curation_contacts.py` → `join_qc.py` → `apply_curation.py`) and the gene-set decision rule (`06f`) are documented in `protocol/PROTOCOL_v3.md`.

## Licence
MIT (`LICENSE`).
