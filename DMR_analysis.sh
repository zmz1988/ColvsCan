#job name
#$ -N motif
#name of logfile 
#$ -o motif.log
# Set the working directory 
#$ -wd /SAN/mottlab/ArabidopsisAssembly/Zhelyazkov/Ziming_analysis/Nanopore_analysis_all
#memory and runtime options 
#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=72:00:00 
#$ -S /bin/bash
#$ -j y
# -pe smp 2 

echo $SGE_TASK_ID $(date)

#mkdir DMR_calling
conda activate megalodon
#megalodon_extras per_read_text modified_bases ./Col_megalodon/
#megalodon_extras per_read_text modified_bases ./Can_megalodon_Col_reference/
#megalodon_extras per_read_text modified_bases ./Can_megalodon/
#megalodon_extras modified_bases create_motif_bed --motif CG 0 ./reference/Col_final_23022022_new_hapog_patch_hapog_reorder.fasta --out-filename Col_CG_sites.bed
#megalodon_extras modified_bases create_motif_bed --motif CH 0 ./reference/Col_final_23022022_new_hapog_patch_hapog_reorder.fasta --out-filename Col_CH_sites.bed

#megalodon_extras modified_bases create_motif_bed --motif CG 0 ./reference/Can_final_23022022_new_hapog.fasta --out-filename Can_CG_sites.bed
#megalodon_extras modified_bases create_motif_bed --motif CH 0 ./reference/Can_final_23022022_new_hapog.fasta --out-filename Can_CH_sites.bed

bedtools intersect -a ./Col_megalodon/modified_bases.5mC.bed -b Col_CG_sites.bed > CG_Col_modified.bed
bedtools intersect -a ./Can_megalodon_Col_reference/modified_bases.5mC.bed -b Col_CG_sites.bed > CG_Can_colref_modified.bed

bedtools intersect -a ./Col_megalodon/modified_bases.5mC.bed -b Col_CH_sites.bed > CH_Col_modified.bed
bedtools intersect -a ./Can_megalodon_Col_reference/modified_bases.5mC.bed -b Col_CH_sites.bed > CH_Can_colref_modified.bed

bedtools intersect -a ./Can_megalodon/modified_bases.5mC.bed -b Can_CG_sites.bed > CG_Can_modified.bed
bedtools intersect -a ./Can_megalodon/modified_bases.5mC.bed -b Can_CH_sites.bed > CH_Can_modified.bed

echo $SGE_TASK_ID $(date)




