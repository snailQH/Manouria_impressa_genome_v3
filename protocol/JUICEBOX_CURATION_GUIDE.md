# Manual Hi-C curation of the *M. impressa* v3 scaffolds in Juicebox

Why: automated scaffolders (YaHS, HapHiC) order and orient contigs from Hi-C contact frequency; they cannot see
every misjoin, missed join, inversion or retained haplotig. Manual inspection of the contact map is the standard
final step for reference-quality assemblies (VGP, Darwin Tree of Life, EBP) and is what prior peer review asked for.
For a hifiasm + Hi-C assembly with contig N50 ~ 90 Mb, expect a few edits, not dozens — but the pass must be
done and reported ("Scaffolds were manually inspected in Juicebox Assembly Tools; N misjoins and M missed joins
were corrected").

## 1. What you need (produced by `scripts/genome_v3/03_scaffold_yahs.sbatch`)

| file | from | purpose |
|---|---|---|
| `out_JBAT.hic` | `juicer pre -a` + `juicer_tools pre` | the contact map Juicebox displays |
| `out_JBAT.assembly` | `juicer pre -a` | contig layout Juicebox edits |
| `out_JBAT.liftover.agp` | `juicer pre -a` | needed later to apply edits |
| `yahs.bin`, contigs FASTA | YaHS | needed later to apply edits |

Location on the server: `01.Genome_Assembly/new-version/03_scaffold/<label>/`. Copy `out_JBAT.hic` and
`out_JBAT.assembly` to a workstation (the .hic is a few GB).

Software: **Juicebox Desktop 1.11.08** (the last version with Assembly Tools; later 2.x builds dropped them):
https://github.com/aidenlab/Juicebox/wiki/Download — "Juicebox with assembly tools". Needs Java 8–11 and ≥ 16 GB RAM
for a 2.4-Gb genome. On a Mac with 8 GB, use a workstation instead.

## 2. Load

1. Open Juicebox → *File → Open* → choose `out_JBAT.hic`.
2. *Assembly → Import Map Assembly* → choose `out_JBAT.assembly`. Scaffold boundaries appear as green lines,
   contig boundaries as blue.
3. Set normalisation to **Balanced** (bottom-left) and resolution to 1–2.5 Mb for the whole genome; zoom to
   100–250 kb for detail. Colour range: drag the slider until the diagonal is saturated and inter-scaffold
   background is faint.
4. *Assembly → Enable Edit Mode*. Save often (*Assembly → Export Assembly*), each time to a new file name.

## 3. What to look for — the four signatures

| signature on the map | meaning | fix |
|---|---|---|
| A scaffold's diagonal is interrupted: a block ends abruptly and the next block has **no contact** with it (contact "drops to background" at a contig boundary, often with a strong off-diagonal block elsewhere) | **misjoin** — two unrelated contigs joined | cut at the boundary, move the piece to where its off-diagonal block says it belongs |
| Two separate scaffolds whose *ends* show a strong square of contact with each other | **missed join** | drag one next to the other |
| Within a scaffold, a block whose contact with its neighbours forms a "bow-tie" / anti-diagonal pattern instead of continuing the diagonal | **inversion** of a contig | select the block, invert (right-click → invert / press *I*) |
| A faint second diagonal running parallel to the main one, or a small contig whose contacts mirror a region already present at full strength | **retained haplotig / duplicated region** | move it to "debris" (end of the assembly) so it becomes unplaced; it will be reported as a haplotypic duplicate |
| Strong contact between a small unplaced contig and one position on a chromosome | **unplaced contig that can be localised** | drag it in; orient by which end contacts more strongly |

Also check: the **number of chromosome-scale scaffolds** against the karyotype. Testudinidae are 2n = 52
(26 chromosome pairs; e.g. *Manouria emys*, *Gopherus*, *Chelonoidis*), so expect ~26 large scaffolds, with
the microchromosomes small. v2.4 reported 32 "chromosomes", which already exceeds 26 and suggests some
chromosome arms or microchromosomes were split — a specific thing to look for.

## 4. Editing operations (Juicebox Assembly Tools)

- **Select** a contig or run of contigs: click; shift-click to extend.
- **Cut** a contig at the cursor: hover inside the contig on the diagonal, press *Ctrl* (Windows/Linux) or
  *Cmd* (Mac) and click — a new blue boundary appears. Only cut where the contact map clearly breaks; hifiasm
  contigs are usually correct, so a cut inside a contig implies a hifiasm misassembly and should be noted.
- **Move**: drag the selection to a new position (a vertical green line shows the insertion point).
- **Invert**: with a selection, right-click → *Invert* (or *I*).
- **Translocate to debris**: drag to the end of the assembly, or right-click → *Move to debris*.
- **Chromosome boundaries**: *Assembly* edit mode treats each green block as a scaffold; to split a scaffold,
  place the cursor on the boundary and press *Shift*+click; to merge, delete the green line.
- **Undo** is limited; keep numbered exports (`curation_01.assembly`, `curation_02.assembly`, …).

Work chromosome by chromosome, largest first, and keep a log table (the Methods will cite it):

| scaffold | position (Mb) | issue | action | evidence |
|---|---|---|---|---|
| scaffold_3 | 112.4 | misjoin at contig boundary ptg000041l/ptg000207l | cut; ptg000207l moved to scaffold_11 end | off-diagonal block, 1-Mb resolution |

## 5. Apply the edits and regenerate the FASTA

On the server, in the same `03_scaffold/<label>/` directory, with the exported `curated.review.assembly`:

```bash
source /data1/opt/conda/etc/profile.d/conda.sh; conda activate /data1/opt/conda/envs/assembly
juicer post -o yahs_curated curated.review.assembly out_JBAT.liftover.agp <contigs.fa>
# -> yahs_curated.FINAL.fa, yahs_curated.FINAL.agp
```

Then re-run the QC and the contact map on the curated FASTA (`04_qc.sbatch` and the Pretext step of
`03_scaffold_yahs.sbatch`) — the post-curation map is the one that goes into the paper, and BUSCO/Merqury
must be re-reported for the curated sequence, not the automated one.

## 6. Name and orient chromosomes

After curation: sort chromosome-scale scaffolds by length, name `chr1…chrN`; orient each so that the
short arm is first where the Hi-C map shows a centromere (contact "bow-tie" at the centromere), and align
to a related chromosome-level assembly (e.g. *Gopherus evgoodei* rGopEvg1 or *Chelonoidis abingdonii*) with
minimap2 `-x asm20` to confirm synteny and to adopt consistent orientation; report the synteny dot plot as a
supplementary figure. Unplaced scaffolds keep their contig names with a `scaffold_` prefix.

## 7. What to write in the manuscript

"Scaffolds were manually inspected and curated in Juicebox Assembly Tools v1.11.08 using the Hi-C contact map
at 1 Mb–100 kb resolution. N misjoins were broken, M scaffolds were joined, K contigs were inverted and J
haplotypic duplicates were removed to unplaced sequences (Supplementary Table SX lists every edit with its
coordinates). The curated map (Fig. 2C) shows no off-diagonal interaction blocks between chromosomes. The
26 largest scaffolds correspond to the 2n = 52 karyotype of Testudinidae."

## 8. Alternative: PretextView

The Darwin Tree of Life uses PretextView (https://github.com/wtsi-hpag/PretextView) on `.pretext` maps, which
the pipeline also produces. It is lighter than Juicebox and edits are applied with `rapid_pretext2tpf` /
`pretext-to-tpf`. Juicebox is recommended here only because prior peer review named it; either is acceptable if the
Methods state which was used and list the edits.
