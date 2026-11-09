#!/bin/bash
# Gathers the small result files the manuscript figures/tables are built from into one folder on the server
# (new-version/09_manuscript_inputs/), computing the per-chromosome summaries that need the big files. Then scp the folder.
set -euo pipefail
NV=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version; O=$NV/09_manuscript_inputs; mkdir -p $O; cd $O
source /data1/opt/conda/etc/profile.d/conda.sh; conda activate /data1/opt/conda/envs/assembly
R=$NV/06_release
cp $R/Mimp_v3.fa.fai $R/Mimp_v3_hap1.fa.fai $R/Mimp_v3_hap2.fa.fai $R/Mimp_v3.agp $R/Mimp_v3.names.tsv $R/Mimp_v3.hic_stats.tsv $R/Mimp_v3.telomeres.tsv $R/software_versions.tsv $R/v2.4_to_v3_correspondence.tsv $R/Mimp_v3_hap1_to_primary.tsv $R/Mimp_v3_hap2_to_primary.tsv $R/release_md5.txt .
samtools faidx $NV/04_qc/v2.4_previous/Manouria_impressa_v2.genome.fasta 2>/dev/null || true; cp $NV/04_qc/v2.4_previous/Manouria_impressa_v2.genome.fasta.fai v2.4.fa.fai
cp $NV/00_data_qc/genomescope2_k21/{linear_plot.png,log_plot.png,transformed_linear_plot.png,summary.txt} . ; for f in linear_plot log_plot transformed_linear_plot; do mv $f.png genomescope_$f.png; done
cp $R/Mimp_v3.stats.tsv $R/Mimp_v3_hap1.stats.tsv $R/Mimp_v3_hap2.stats.tsv $R/fcs_gx_actions_applied.tsv . 2>/dev/null || true
cp $NV/04_qc/haplotypes_clean/merqury_run/merqury.qv merqury_haplotypes.qv 2>/dev/null; cp $NV/04_qc/haplotypes_clean/merqury_run/merqury.completeness.stats merqury_haplotypes.completeness.stats 2>/dev/null
cp $NV/05_annotation/Mimp_v3/04_ncrna/hc/trna_hc.log trna_hc.log 2>/dev/null; cp $NV/05_annotation/Mimp_v3/06_final_geneset/{gene_set_stats.tsv,eggnog_rate.txt,rnaseq_support_final.txt} . 2>/dev/null; cp $NV/05_annotation/Mimp_v3/06_final_geneset/busco_final/short_summary*.txt busco_protein_final.txt 2>/dev/null; cp $NV/05_annotation/Mimp_v3/07_post_geneset/{gene_density.tsv,gh18_best_hits.tsv,gh18_region_genes.tsv} . 2>/dev/null; cp $NV/05_annotation/Mimp_v3/05_annotation_check/busco_summary.txt busco_protein_candidates.txt 2>/dev/null; cp $NV/05_annotation/Mimp_v3/04_ncrna/rrna.gff rrna.gff 2>/dev/null
cp $NV/04_qc/Mimp_v3/merqury_run/merqury.asm1.spectra-cn.fl.png merqury_Mimp_v3_spectra-cn.png; cp $NV/04_qc/Mimp_v3/merqury_run/merqury.qv merqury_Mimp_v3.qv; cp $NV/04_qc/Mimp_v3/merqury_run/merqury.asm1.qv merqury_Mimp_v3_per_sequence.qv; cp $NV/04_qc/Mimp_v3/merqury_run/merqury.completeness.stats merqury_Mimp_v3.completeness.stats
cp $NV/03_scaffold/contactmap_Mimp_v3/snapshots/Mimp_v3_FullMap.png hic_map_Mimp_v3.png; for h in hap1 hap2; do cp $NV/03_scaffold/contactmap_Mimp_v3_$h/snapshots/Mimp_v3_${h}_FullMap.png hic_map_Mimp_v3_$h.png 2>/dev/null || true; done; cp $NV/03_scaffold/synteny_Mimp_v3/dotplot_Mimp_v3.png dotplot_Mimp_v3_vs_GopEvg.png
cp $NV/03_scaffold/curated/v3_curated.edits.log curation_edits.log; cp $NV/03_scaffold/contactmap_primary_yahs/join_qc_yahs.tsv join_qc_before_curation.tsv; cp $NV/03_scaffold/contactmap_v3_curated/join_qc_v3_curated.tsv join_qc_after_curation.tsv
cp $NV/04_qc/hificov_Mimp_v3/zero_coverage_runs.bed coverage_zero_runs.bed; cp $NV/04_qc/hap_divergence/hap_divergence_summary.tsv .
cp $NV/07_contamination/Mimp_v3/mito/final_mitogenome.annotation.png mito_annotation.png; cp $NV/07_contamination/Mimp_v3/mito/mito_array_copynumber_hist.tsv .; cp $NV/06_release/Mimp_v3_mito.gff3 mito_annotation.gff3; cp $NV/06_release/Mimp_v3.genes.gwh_partial_models.tsv gwh_partial_models.tsv; cp $NV/06_release/Mimp_v3_mito.gb mito_annotation.gb
for a in Mimp_v3 hap1_scaffolds hap2_scaffolds v2.4_previous; do for l in vertebrata_odb12 sauropsida_odb12; do cp $NV/04_qc/$a/compleasm_$l/summary.txt compleasm_${a}_$l.txt 2>/dev/null || true; done; cp $NV/04_qc/$a/busco_vertebrata/short_summary*.txt busco_${a}.txt 2>/dev/null || true; done
# per-chromosome SNV density from the haplotype-vs-primary calls
for h in hap1 hap2; do (echo -e "sequence\taligned_bp\tSNV\tvariants_all\tSNV_per_kb"; awk 'BEGIN{OFS="\t"} $1=="V"{v[$2]++; if(length($7)==1&&length($8)==1&&$7!="-"&&$8!="-")s[$2]++} $1=="R"{c[$2]+=$4-$3} END{for(k in c) if(c[k]>0) print k,c[k],s[k]+0,v[k]+0,1000*(s[k]+0)/c[k]}' $NV/04_qc/hap_divergence/${h}_vs_primary.var | sort -k1,1V) > snv_density_${h}_vs_primary.tsv; done
# chromosome table pieces: contigs and gaps per sequence from the AGP; v2.4 names per v3 chromosome
# contigs/gaps per release sequence: the release AGP maps chromosomes to curated scaffolds; the curated AGP maps scaffolds to contigs
cp $NV/03_scaffold/curated/v3_curated.agp curated.agp
(echo -e "sequence\tcontigs\tgaps"; awk -F'\t' 'BEGIN{OFS="\t"} FNR==NR && FNR>1 {rel[$2]=$1; next} /^#/{next} $5=="W"{c[rel[$1]]++} $5=="U"||$5=="N"{g[rel[$1]]++} END{for(k in c) if(k!="") print k,c[k],g[k]+0}' $R/Mimp_v3.names.tsv $NV/03_scaffold/curated/v3_curated.agp | sort -k1,1V) > agp_components.tsv
(echo -e "v3_chr\tv2.4_sequences_mapping_>=50%"; awk -F'\t' 'NR>1 && $3 ~ /^chr/ && $5>=50 && $1 !~ /^chrUn/ {a[$3]=a[$3] (a[$3]?", ":"") $1 " (" $5 "%)"} END{for(k in a) print k"\t"a[k]}' $R/v2.4_to_v3_correspondence.tsv | sort -k1,1V) > v24_names_per_v3_chr.tsv
cp $NV/05_annotation/Mimp_v3/04_ncrna/trna_hc_isotypes.txt $NV/05_annotation/Mimp_v3/04_ncrna/hc/trna_hc.summary . 2>/dev/null || true
ls -la | wc -l; echo "collected $(date)"
