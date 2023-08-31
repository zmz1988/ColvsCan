ln -s /SAN/mottlab/ArabidopsisAssembly/rsync_biomyriad/annotation_update/OMA_analysis/Can_new_annotation_clean1_uniprot_function.gff3 ./
ln -s /SAN/mottlab/ArabidopsisAssembly/rsync_biomyriad/annotation_update/OMA_analysis/Col_new_annotation_clean1_uniprot_function.gff3 ./

awk '$3=="gene"' Can_new_annotation_clean1_uniprot_function.gff3 | awk 'BEGIN{OFS="\t"} {print $1,$4,$5,$7,$9}' | sed 's/;uniprot.*//g' | sed 's/;low_identity.*//g' |  sed 's/ID=//g' > Can_new_annotation_clean1_uniprot_function.gene.interval.bed

awk '$3=="gene"' Col_new_annotation_clean1_uniprot_function.gff3 | awk 'BEGIN{OFS="\t"} {print $1,$4,$5,$7,$9}' | sed 's/;Name.*//g' | sed 's/ID=//g' > Col_new_annotation_clean1_uniprot_function.gene.interval.bed

sort -k1,1 -k2,2n ../Can_megalodon/modified_bases.5mC.bed > Can_megalodon_modified_bases.5mC.sorted.bed
sort -k1,1 -k2,2n ../Col_megalodon/modified_bases.5mC.bed > Col_megalodon_modified_bases.5mC.sorted.bed
sort -k1,1 -k2,2n Can_new_annotation_clean1_uniprot_function.gene.interval.bed > Can_new_annotation_clean1_uniprot_function.gene.interval.sorted.bed
sort -k1,1 -k2,2n Col_new_annotation_clean1_uniprot_function.gene.interval.bed > Col_new_annotation_clean1_uniprot_function.gene.interval.sorted.bed

bedtools map -a Can_new_annotation_clean1_uniprot_function.gene.interval.sorted.bed -b Can_megalodon_modified_bases.5mC.sorted.bed -c 11,11 -o mean,count > Can_gene_level_methylation
bedtools map -a Col_new_annotation_clean1_uniprot_function.gene.interval.sorted.bed -b Col_megalodon_modified_bases.5mC.sorted.bed -c 11,11 -o mean,count > Col_gene_level_methylation

bedtools map -a /SAN/mottlab/ArabidopsisAssembly/multi_omics_analysis_Col_Can/TE_distance/Col_gene_upstream_100bp.bed -b Col_megalodon_modified_bases.5mC.sorted.bed -c 11,11 -o mean,count > Col_gene_upstream_100bp_methylation
bedtools map -a /SAN/mottlab/ArabidopsisAssembly/multi_omics_analysis_Col_Can/TE_distance/Col_gene_upstream_500bp.bed -b Col_megalodon_modified_bases.5mC.sorted.bed -c 11,11 -o mean,count > Col_gene_upstream_500bp_methylation
bedtools map -a /SAN/mottlab/ArabidopsisAssembly/multi_omics_analysis_Col_Can/TE_distance/Col_gene_upstream_1000bp.bed -b Col_megalodon_modified_bases.5mC.sorted.bed -c 11,11 -o mean,count > Col_gene_upstream_1000bp_methylation
bedtools map -a /SAN/mottlab/ArabidopsisAssembly/multi_omics_analysis_Col_Can/TE_distance/Col_gene_upstream_2000bp.bed -b Col_megalodon_modified_bases.5mC.sorted.bed -c 11,11 -o mean,count > Col_gene_upstream_2000bp_methylation
bedtools map -a /SAN/mottlab/ArabidopsisAssembly/multi_omics_analysis_Col_Can/TE_distance/Col_gene_upstream_5000bp.bed -b Col_megalodon_modified_bases.5mC.sorted.bed -c 11,11 -o mean,count > Col_gene_upstream_5000bp_methylation

bedtools map -a /SAN/mottlab/ArabidopsisAssembly/multi_omics_analysis_Col_Can/TE_distance/Can_gene_upstream_100bp.bed -b Can_megalodon_modified_bases.5mC.sorted.bed -c 11,11 -o mean,count > Can_gene_upstream_100bp_methylation
bedtools map -a /SAN/mottlab/ArabidopsisAssembly/multi_omics_analysis_Col_Can/TE_distance/Can_gene_upstream_500bp.bed -b Can_megalodon_modified_bases.5mC.sorted.bed -c 11,11 -o mean,count > Can_gene_upstream_500bp_methylation
bedtools map -a /SAN/mottlab/ArabidopsisAssembly/multi_omics_analysis_Col_Can/TE_distance/Can_gene_upstream_1000bp.bed -b Can_megalodon_modified_bases.5mC.sorted.bed -c 11,11 -o mean,count > Can_gene_upstream_1000bp_methylation
bedtools map -a /SAN/mottlab/ArabidopsisAssembly/multi_omics_analysis_Col_Can/TE_distance/Can_gene_upstream_2000bp.bed -b Can_megalodon_modified_bases.5mC.sorted.bed -c 11,11 -o mean,count > Can_gene_upstream_2000bp_methylation
bedtools map -a /SAN/mottlab/ArabidopsisAssembly/multi_omics_analysis_Col_Can/TE_distance/Can_gene_upstream_5000bp.bed -b Can_megalodon_modified_bases.5mC.sorted.bed -c 11,11 -o mean,count > Can_gene_upstream_5000bp_methylation


bedtools map -a /SAN/mottlab/ArabidopsisAssembly/multi_omics_analysis_Col_Can/TE_distance/Col_gene_downstream_50bp.bed -b Col_megalodon_modified_bases.5mC.sorted.bed -c 11,11 -o mean,count > Col_gene_downstream_50bp_methylation
bedtools map -a /SAN/mottlab/ArabidopsisAssembly/multi_omics_analysis_Col_Can/TE_distance/Col_gene_downstream_100bp.bed -b Col_megalodon_modified_bases.5mC.sorted.bed -c 11,11 -o mean,count > Col_gene_downstream_100bp_methylation
bedtools map -a /SAN/mottlab/ArabidopsisAssembly/multi_omics_analysis_Col_Can/TE_distance/Col_gene_downstream_500bp.bed -b Col_megalodon_modified_bases.5mC.sorted.bed -c 11,11 -o mean,count > Col_gene_downstream_500bp_methylation
bedtools map -a /SAN/mottlab/ArabidopsisAssembly/multi_omics_analysis_Col_Can/TE_distance/Col_gene_downstream_1000bp.bed -b Col_megalodon_modified_bases.5mC.sorted.bed -c 11,11 -o mean,count > Col_gene_downstream_1000bp_methylation

bedtools map -a /SAN/mottlab/ArabidopsisAssembly/multi_omics_analysis_Col_Can/TE_distance/Can_gene_downstream_50bp.bed -b Can_megalodon_modified_bases.5mC.sorted.bed -c 11,11 -o mean,count > Can_gene_downstream_50bp_methylation
bedtools map -a /SAN/mottlab/ArabidopsisAssembly/multi_omics_analysis_Col_Can/TE_distance/Can_gene_downstream_100bp.bed -b Can_megalodon_modified_bases.5mC.sorted.bed -c 11,11 -o mean,count > Can_gene_downstream_100bp_methylation
bedtools map -a /SAN/mottlab/ArabidopsisAssembly/multi_omics_analysis_Col_Can/TE_distance/Can_gene_downstream_500bp.bed -b Can_megalodon_modified_bases.5mC.sorted.bed -c 11,11 -o mean,count > Can_gene_downstream_500bp_methylation
bedtools map -a /SAN/mottlab/ArabidopsisAssembly/multi_omics_analysis_Col_Can/TE_distance/Can_gene_downstream_1000bp.bed -b Can_megalodon_modified_bases.5mC.sorted.bed -c 11,11 -o mean,count > Can_gene_downstream_1000bp_methylation

awk '$3=="exon"' Can_new_annotation_clean1_uniprot_function.gff3 | awk 'BEGIN{OFS="\t"} {print $1,$4,$5,$7,$9}' | sed 's/;Parent.*//g' | sed 's/;low_identity.*//g' |  sed 's/ID=//g' > Can_new_annotation_clean1_uniprot_function.exon.interval.bed

awk '$3=="exon"' Col_new_annotation_clean1_uniprot_function.gff3 | awk 'BEGIN{OFS="\t"} {print $1,$4,$5,$7,$9}' | sed 's/;Parent.*//g' | sed 's/ID=//g' > Col_new_annotation_clean1_uniprot_function.exon.interval.bed

sort -k1,1 -k2,2n Col_new_annotation_clean1_uniprot_function.exon.interval.bed > Col_new_annotation_clean1_uniprot_function.exon.interval.sorted.bed
sort -k1,1 -k2,2n Can_new_annotation_clean1_uniprot_function.exon.interval.bed > Can_new_annotation_clean1_uniprot_function.exon.interval.sorted.bed

bedtools map -a Can_new_annotation_clean1_uniprot_function.exon.interval.sorted.bed -b Can_megalodon_modified_bases.5mC.sorted.bed -c 11,11 -o mean,count > Can_exon_level_methylation
bedtools map -a Col_new_annotation_clean1_uniprot_function.exon.interval.sorted.bed -b Col_megalodon_modified_bases.5mC.sorted.bed -c 11,11 -o mean,count > Col_exon_level_methylation

bedtools map -a /SAN/mottlab/ArabidopsisAssembly/multi_omics_analysis_Col_Can/Col_vs_Can_multi_factor/Can_new_annotation_clean1_intron.bed -b Can_megalodon_modified_bases.5mC.sorted.bed -c 11,11 -o mean,count > Can_intron_level_methylation

bedtools map -a /SAN/mottlab/ArabidopsisAssembly/multi_omics_analysis_Col_Can/Col_vs_Can_multi_factor/Col_new_annotation_clean1_intron.bed -b Col_megalodon_modified_bases.5mC.sorted.bed -c 11,11 -o mean,count > Col_intron_level_methylation

awk '$3=="mRNA"' Can_new_annotation_clean1_uniprot_function.gff3 | awk 'BEGIN{OFS="\t"} {print $1,$4,$5,$7,$9}' |sed 's/;Parent.*//g' |  sed 's/ID=//g' > Can_new_annotation_clean1_uniprot_function.transcripts.interval.bed

awk '$3=="mRNA"' Col_new_annotation_clean1_uniprot_function.gff3 | awk 'BEGIN{OFS="\t"} {print $1,$4,$5,$7,$9}' |sed 's/;Parent.*//g' |  sed 's/ID=//g' > Col_new_annotation_clean1_uniprot_function.transcripts.interval.bed

sort -k1,1 -k2,2n Col_new_annotation_clean1_uniprot_function.transcripts.interval.bed > Col_new_annotation_clean1_uniprot_function.transcripts.interval.sorted.bed
sort -k1,1 -k2,2n Can_new_annotation_clean1_uniprot_function.transcripts.interval.bed > Can_new_annotation_clean1_uniprot_function.transcripts.interval.sorted.bed

bedtools map -a Can_new_annotation_clean1_uniprot_function.transcripts.interval.sorted.bed -b Can_megalodon_modified_bases.5mC.sorted.bed -c 11,11 -o mean,count > Can_transcripts_level_methylation
bedtools map -a Col_new_annotation_clean1_uniprot_function.transcripts.interval.sorted.bed -b Col_megalodon_modified_bases.5mC.sorted.bed -c 11,11 -o mean,count > Col_transcripts_level_methylation

ln -s /SAN/mottlab/ArabidopsisAssembly/rsync_biomyriad/RNAseq_expression/Col_withCanRef/Col_withCanRef_kallisto_outdir/kallisto.gene.TMM.EXPR.matrix ./Col_withCanRef_kallisto.gene.TMM.EXPR.matrix
ln -s /SAN/mottlab/ArabidopsisAssembly/rsync_biomyriad/RNAseq_expression/Col_withCanRef/Col_withCanRef_kallisto_outdir/kallisto.gene.TPM.not_cross_norm ./Col_withCanRef_kallisto.gene.TPM.not_cross_norm
ln -s /SAN/mottlab/ArabidopsisAssembly/rsync_biomyriad/RNAseq_expression/Col_withCanRef/Col_withCanRef_kallisto_outdir/kallisto.isoform.TMM.EXPR.matrix ./Col_withCanRef_kallisto.isoform.TMM.EXPR.matrix
ln -s /SAN/mottlab/ArabidopsisAssembly/rsync_biomyriad/RNAseq_expression/Col_withCanRef/Col_withCanRef_kallisto_outdir/kallisto.isoform.TPM.not_cross_norm ./Col_withCanRef_kallisto.isoform.TPM.not_cross_norm

ln -s /SAN/mottlab/ArabidopsisAssembly/rsync_biomyriad/RNAseq_expression/Can_withColRef/Can_withColRef_kallisto_outdir/kallisto.gene.TMM.EXPR.matrix ./Can_withColRef_kallisto.gene.TMM.EXPR.matrix
ln -s /SAN/mottlab/ArabidopsisAssembly/rsync_biomyriad/RNAseq_expression/Can_withColRef/Can_withColRef_kallisto_outdir/kallisto.gene.TPM.not_cross_norm ./Can_withColRef_kallisto.gene.TPM.not_cross_norm
ln -s /SAN/mottlab/ArabidopsisAssembly/rsync_biomyriad/RNAseq_expression/Can_withColRef/Can_withColRef_kallisto_outdir/kallisto.isoform.TMM.EXPR.matrix ./Can_withColRef_kallisto.isoform.TMM.EXPR.matrix
ln -s /SAN/mottlab/ArabidopsisAssembly/rsync_biomyriad/RNAseq_expression/Can_withColRef/Can_withColRef_kallisto_outdir/kallisto.isoform.TPM.not_cross_norm ./Can_withColRef_kallisto.isoform.TPM.not_cross_norm

python3 combine_methylation_rna.py Can_gene_level_methylation Col_withCanRef_kallisto.gene.TMM.EXPR.matrix
python3 combine_methylation_rna.py Can_gene_level_methylation Col_withCanRef_kallisto.gene.TPM.not_cross_norm
python3 combine_methylation_rna.py Col_gene_level_methylation Can_withColRef_kallisto.gene.TMM.EXPR.matrix
python3 combine_methylation_rna.py Col_gene_level_methylation Can_withColRef_kallisto.gene.TPM.not_cross_norm



bedtools map -a Can_new_annotation_clean1_uniprot_function.gene.interval.sorted.bed -b Can_megalodon_modified_bases.5mC.sorted.bed -c 11,11,11,11,11,11 -o sum,min,max,median,mean,count > Can_gene_level_methylation_sum
bedtools map -a Col_new_annotation_clean1_uniprot_function.gene.interval.sorted.bed -b Col_megalodon_modified_bases.5mC.sorted.bed -c 11,11,11,11,11,11 -o sum,min,max,median,mean,count > Col_gene_level_methylation_sum
