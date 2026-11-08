# Data file 07 mitogenome

Mitochondrial genome annotation (MitoHiFi) and read support for the tRNA-Thr/Pro tandem array.

Genome: *Manouria impressa* Mimp_v3 (BioProject PRJCA038476; BioSample SAMC4973488 = liver DNA libraries; RNA-seq BioSamples SAMC4973532 heart, SAMC4976828 lung). Produced by the pipeline in Data file 13 (scripts/genome_v3). mito_annotation.gff3/.gb: the 43 genes (13 CDS, 2 rRNA, 28 tRNA incl. the four-unit Thr/Pro array) + control region, GFF3 (17_mito_gb2gff3.py from the MitoHiFi GenBank). mito_array_copynumber_hist: net indel length (bp, 100-bp bins) within the array per spanning HiFi read, read count. 0 = four units (consensus), +300 = five units.

| file | bytes | md5 |
|---|---|---|
| mito_annotation.png | 80,980 | 44152fa36e508a576aacea812d0b9783 |
| mito_annotation.gff3 | 11,303 | d78886f79baa5a320c7ffcc2802b6298 |
| mito_annotation.gb | 35,052 | 01ec698bfec1cc400cc8e449dbc957e1 |
| mito_array_copynumber_hist.tsv | 111 | 6b9515f19561496b573cc134b25230e0 |
