# Data file 07 mitogenome

Mitochondrial genome annotation (MitoHiFi) and read support for the tRNA-Thr/Pro tandem array.

Genome: *Manouria impressa* Mimp_v3 (BioProject PRJCA038476; BioSample SAMC4973488 = liver DNA libraries; RNA-seq BioSamples SAMC4973532 heart, SAMC4976828 lung). Produced by the pipeline in Data file 13 (scripts/genome_v3). mito_annotation.gff3/.tbl: curated annotation (17b_mito_annotation_curated.py; 13 CDS, 2 rRNA, 28 tRNA incl. the four-unit Thr/Pro array, control region) — COX1 GTG start at 5404, ND4 start 10255, ND3 +1 frameshift join(9545..9718,9720..9894), incomplete stops as transl_except, checked against NC_011815.1; .gb = the original MitoHiFi call. mito_array_copynumber_hist: net indel length (bp, 100-bp bins) within the array per spanning HiFi read, read count. 0 = four units (consensus), +300 = five units.

| file | bytes | md5 |
|---|---|---|
| mito_annotation.png | 80,980 | 44152fa36e508a576aacea812d0b9783 |
| mito_annotation.gff3 | 11,513 | 16f2fc66213d380b7ec95bda60932ee9 |
| mito_annotation.tbl | 5,241 | b06af075671149ea149453fcfb058666 |
| mito_annotation.gb | 35,052 | 01ec698bfec1cc400cc8e449dbc957e1 |
| mito_array_copynumber_hist.tsv | 111 | 6b9515f19561496b573cc134b25230e0 |
