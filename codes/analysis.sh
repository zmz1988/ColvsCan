# generate homology from the HOG group
# in /tmp/Ziming_analysis/annotation_update/OMA_synteny/Output/HOGFasta/ folder
python3 Get_gene_List_from_HOG.py > HOG_gene_list.txt
cp /tmp/Ziming_analysis/annotation_update/OMA_synteny/Output/HOGFasta/HOG_gene_list.txt ./

##################### starting here from the final analysis that are used #########
# in /tmp/Ziming_analysis/annotation_update/OMA_analysis


############# delete genes and transcripts that have the identical coordinates ########3
python3 delete_identical_gene_and_transcript_in_gff.py # to prepare the gene/trnascripts that need to be changed
# generated Col_identical_gene_to_delete_gff.txt, Can_identical_gene_to_delete_gff.txt, Col_identical_gene_transcript_name_change_gff.txt, Can_identical_gene_transcript_name_change_gff.txt
awk 'NF' Can_identical_gene_transcript_name_change_gff.txt > Can_identical_gene_transcript_name_change_gff1.txt  # delete the empty lines
awk 'NF' Can_identical_gene_to_delete_gff.txt >  Can_identical_gene_to_delete_gff1.txt
awk 'NF' Col_identical_gene_transcript_name_change_gff.txt > Col_identical_gene_transcript_name_change_gff1.txt
awk 'NF' Col_identical_gene_to_delete_gff.txt > Col_identical_gene_to_delete_gff1.txt
sed -i 's/ /\t/g' Col_identical_gene_transcript_name_change_gff1.txt
sed -i 's/ /\t/g' Can_identical_gene_transcript_name_change_gff1.txt
sed -i 's/ /\t/g' Can_identical_gene_to_delete_gff1.txt
sed -i 's/ /\t/g' Col_identical_gene_to_delete_gff1.txt
# change the gff file
python3 change_the_gff1.py  # generate two new .gff3 files Can_new_annotation.gff3 and Col_new_annotation.gff3
sort -k1,1 -k4,4 -k5,5 -k6,6 Can_new_annotation.gff3 | uniq > Can_new_annotation_clean.gff3
sort -k1,1 -k4,4 -k5,5 -k6,6 Col_new_annotation.gff3 | uniq > Col_new_annotation_clean.gff3
agat_sp_fix_features_locations_duplicated.pl --gff Can_new_annotation_clean.gff3 --out Can_new_annotation_clean1.gff3
agat_sp_fix_features_locations_duplicated.pl --gff Col_new_annotation_clean.gff3 -o Col_new_annotation_clean1.gff3
# note: if directly use agat to dedup the interpro.gff file, the duplicated genes will not be merged or deleted

grep 'case1 (removing mRNA isoform identic ): ' dedup_Can_new_annotation_gff.log | sed 's/case1 (removing mRNA isoform identic ): //g' | sed 's/c/C/g' | sed 's/,/\n/g' > dedup_Can_new_annotation_gff_transcripts_name.txt
grep 'case1 (removing mRNA isoform identic ): ' dedup_Col_new_annotation_gff.log | sed 's/case1 (removing mRNA isoform identic ): //g' | sed 's/c/C/g' | sed 's/,/\n/g' > dedup_Col_new_annotation_gff_transcripts_name.txt
python3 generate_changed_transcript_list_for_YongIn.py  # generate two files Col_changed_transcript_list.txt, Can_changed_transcript_list.txt

####################################### above was done in /tmp/Ziming_analysis/annotation_update/OMA_analysis ##########################

##################################### Gene Homology ############
###################### generate 1to1 genes list
python3 pyHam_related_analysis.py  # generate HOG_arabidopsis.txt, HOG_geneID_arabidopsis_1to1.txt, HOG_geneID_arabidopsis.txt, HOG_transcriptID_arabidopsis_1to1.txt, and two gene tree profiles 
python3 generate_1to1_geneID_file.py # to re-generate the HOG_geneID_arabidopsis_1to1.txt
sed -i '/^$/d' HOG_geneID_arabidopsis_1to1.txt
sed -i '/^$/d' HOG_transcriptID_arabidopsis_1to1.txt
python3 clean_oma_1to1.py > HOG_geneID_arabidopsis_1to1_clean.txt
sed -i 's/AT://' HOG_geneID_arabidopsis_1to1_clean.txt
sed -i 's/Cang://' HOG_geneID_arabidopsis_1to1_clean.txt
sed -i 's/Colg://' HOG_geneID_arabidopsis_1to1_clean.txt
awk '{print $4,$3,$2,$1}' HOG_geneID_arabidopsis_1to1_clean.txt > HOG_geneID_arabidopsis_1to1_clean_final.txt # final 27186 line or HOG

grep 'Cang' HOG_geneID_arabidopsis_1to1_clean_final.txt | wc -l # 25071 genes from Can 
grep 'Colg' HOG_geneID_arabidopsis_1to1_clean_final.txt | wc -l # 25121 genes from Col
grep 'AT' HOG_geneID_arabidopsis_1to1_clean_final.txt | wc -l # 25165 genes from Tair
awk '$3=="NA" && $4=="NA"' HOG_geneID_arabidopsis_1to1_clean_final.txt | wc -l # 268 Col specific genes
awk '$2=="NA" && $4=="NA"' HOG_geneID_arabidopsis_1to1_clean_final.txt | wc -l # 666 Can specific genes 
awk '$2=="NA" && $3=="NA"' HOG_geneID_arabidopsis_1to1_clean_final.txt | wc -l # 671 Tair specific genes

# to generate a statistics for the homology between Col, Can, and Tair
python3 statistics.py > statistics.txt  # number indicates that some are duplicated accession-specific genes

##################### Blast duplicated gene
python3 duplication_event.py  # generate input file HOG_genbeID_arabidopsis_duplication_event.txt for blast_duplicated_event.py; all the duplicated gene events were captured
python3 blast_duplicated_event.py  # generate fasta files for blast and HOG_arabidopsis_accession_specific_gene_with_duplication.txt, in which accession specific genes with duplication events were captured
python3 single_species_duplication_event.py  # generate HOG_genbeID_arabidopsis_single_species_duplication_event.txt, in which duplicated gene only in one accession

## NOTE: the blast_duplicated_event.py generates a file for accession specific genes that are duplicated, which genes will not be blasted. The single_species_duplication_event.py generates a list of genes that are not accession specific but only duplicated in one accession. These genes are blasted.

# to generate a specific HOG for each accession specific gene with duplication event
python3 hog_to_accession_specific_gene_with_duplication.py > HOG_arabidopsis_accession_specific_gene_with_duplication_final.txt  # produce file HOG_arabidopsis_accession_specific_gene_with_duplication_clean.txt

sed -i 's/AT://' HOG_arabidopsis_accession_specific_gene_with_duplication_final.txt
sed -i 's/Cang://' HOG_arabidopsis_accession_specific_gene_with_duplication_final.txt
sed -i 's/Colg://' HOG_arabidopsis_accession_specific_gene_with_duplication_final.txt
awk '{print $4,$3,$2,$1}' HOG_arabidopsis_accession_specific_gene_with_duplication_final.txt > HOG_arabidopsis_accession_specific_gene_with_duplication_anno.txt 

# run blast for duplicated genes (_db file feed to makedb, and _blast is for blast)
for file in $(find ./Blast_duplicated_gene -name "*_db*")
do
name=$(echo $file | sed 's/_db.*//')
blast=$(echo $file | sed 's/_db/_blast/')
makeblastdb -in $file -dbtype nucl -parse_seqids -out ${name}_database
blastn -db ${name}_database -query ${blast} -out ${name}_out -outfmt "6 qseqid sseqid evalue bitscore length pident qstart qend qlen sstart send slen" -max_hsps 1
done


for file in $(find ./Blast_duplicated_gene -name "*_self*")
do
name=$(echo $file | sed 's/_self.*//')
makeblastdb -in $file -dbtype nucl -parse_seqids -out ${name}_database
blastn -db ${name}_database -query ${file} -out ${name}_out -outfmt "6 qseqid sseqid evalue bitscore length pident qstart qend qlen sstart send slen" -max_hsps 1
done

################################## clean up blast results ##################
############### for the duplicated gene
####### for HOG s that has no inter-accession homologous or genes without riceprocal blast results
python3 screen_blast_duplicated_gene.py > duplicated_gene_no_homologous_blast.txt # in the output .txt file, lists the genes that have no homologous genes in other accessions
# genes with a HOG ID but can't find inter-accession homologs
grep 'not at all' duplicated_gene_no_homologous_blast.txt | sed 's/.*://' | sed 's/ has.*//g' > HOG_no_inter-accession_homology_duplicated_gene.txt  # 28 HOG outputs
python3 summarise_screen_result_duplicated_gene.py > non_reciprocal_result_after_screen_duplicated_gene.txt # in the .txt file listed genes that have no reciprocal blast results, but has blast results. wc -l results in 143 lines, and wc -w results in 573 words. So there are 573-143=430 genes that has no reciprocal blast results

# change the HOG number digit for Yong-in
python3 chang_HOG_number_digit.py


####### concatenate the reciprocal blast results
cat Blast_duplicated_gene/*_homology > reciprocal_duplicated_gene_homology.txt
python3 clean_reciprocal_duplicated_homology.py > reciprocal_duplicated_gene_homology_clean.txt
sed -i 's/AT://' reciprocal_duplicated_gene_homology_clean.txt
sed -i 's/Cang://' reciprocal_duplicated_gene_homology_clean.txt
sed -i 's/Colg://' reciprocal_duplicated_gene_homology_clean.txt
awk '{print $4,$3,$2,$1}' reciprocal_duplicated_gene_homology_clean.txt > reciprocal_duplicated_gene_homology_clean1.txt # final 1866 line, including some genes have been repeated multiple times, because another accession has multiple same copy of genes

########## generate the final homology file
# concatenate the 1to1 gene with the duplicated gene homology
cat HOG_geneID_arabidopsis_1to1_clean_final.txt reciprocal_duplicated_gene_homology_clean1.txt HOG_arabidopsis_accession_specific_gene_with_duplication_anno.txt  > HOG_gene_1to1_plus_reciprocal_duplicated.txt
sed -i 's/ /\t/g' HOG_gene_1to1_plus_reciprocal_duplicated.txt
grep 'Cang' HOG_gene_1to1_plus_reciprocal_duplicated.txt | sort -u -k3,3 | wc -l # 26141 genes from Can 
grep 'Colg' HOG_gene_1to1_plus_reciprocal_duplicated.txt | sort -u -k2,2 | wc -l # 26543 genes from Col
grep 'AT' HOG_gene_1to1_plus_reciprocal_duplicated.txt | sort -u -k4,4 | wc -l # 26177 genes from Tair
awk '$3=="NA" && $4=="NA"' HOG_gene_1to1_plus_reciprocal_duplicated.txt | wc -l # 289 Col specific genes
awk '$2=="NA" && $4=="NA"' HOG_gene_1to1_plus_reciprocal_duplicated.txt | wc -l # 693 Can specific genes 
awk '$2=="NA" && $3=="NA"' HOG_gene_1to1_plus_reciprocal_duplicated.txt | wc -l # 675 Tair specific genes

python3 take_out_identical_genes_HOG.py  # generate HOG_gene_1to1_plus_reciprocal_duplicated_annotation.txt, which has no identical gene records (identical in any accessions) (including those egens which are not identical or duplicated in other accessions)
wc -l HOG_gene_1to1_plus_reciprocal_duplicated_annotation.txt # 28142 lines
wc -l HOG_gene_1to1_plus_reciprocal_duplicated.txt # 29104, meaning >900 lines have identical genes involved
grep 'Cang' HOG_gene_1to1_plus_reciprocal_duplicated_annotation.txt | wc -l # 25804 genes from Can 
grep 'Colg' HOG_gene_1to1_plus_reciprocal_duplicated_annotation.txt | wc -l # 25987 genes from Col
grep 'AT' HOG_gene_1to1_plus_reciprocal_duplicated_annotation.txt | wc -l # 25915 genes from Tair

awk '{print $1}' reciprocal_duplicated_identical_gene_Col.txt | sed 's/_.*//g' | uniq -c | awk '$1>1' | wc -l # 156 genes in Col have identical duplication
wc -l reciprocal_duplicated_identical_gene_Can.txt # 339 can genes are among these cases
awk '{print $1}' reciprocal_duplicated_identical_gene_Can.txt | sed 's/_.*//g' | uniq -c | awk '$1>1' | wc -l # 66 genes in Can have identical duplication
wc -l reciprocal_duplicated_identical_gene_Tair.txt # 264 tair genes are among these cases
awk '{print $1}' reciprocal_duplicated_identical_gene_Tair.txt | sed 's/_.*//g' | uniq -c | awk '$1>1' | wc -l # 45 genes in Tiar have identical duplication


############ perform GO enrichment for duplicated genes
python3 grab_duplicated_gene_in_each_accession.py # produce the gene lists for each accession Can_duplicated_gene_list.txt (1536), Col_duplicated_gene_list.txt (1874), Tair_duplicated_gene_list.txt (1178 gene for tair) 
ln -s /tmp/Ziming_analysis/annotation_update/Can/blast/Can_braker_pasa_liftoff_coding_ID_tag_Interpro/GO.txt Can_Go.txt
ln -s /tmp/Ziming_analysis/annotation_update/Col/blast/Col_braker_pasa_liftoff_coding_ID_homology_interpro/GO.txt Col_Go.txt
python3 clean_GO_file.py # generate a new input of reference GO file at gene level
# the input interesting gene lists are reciprocal_duplicated_identical_gene_Tair.txt (or Can or Col)
# the input interesting gene lists are also Col_duplicated_gene_list.txt and Can_duplicated_gene_list.txt

#################################### for the transcripts
# prepare fast file for blast
ln -s /tmp/Ziming_analysis/annotation_update/OMA_synteny/DB/Can_new_annotation_clean1_amended1.splice ./
ln -s /tmp/Ziming_analysis/annotation_update/OMA_synteny/DB/Col_new_annotation_clean1_amended1.splice ./
python3 transcript_blast_prep.py
wc -l HOG_transcriptID_arabidopsis_1to1_noIsof.txt  # 10839 HOGs that has no isoform

###### blast

for file in $(find ./Blast_transcripts -name "*_db*")
do
name=$(echo $file | sed 's/_db.*//')
blast=$(echo $file | sed 's/_db/_blast/')
makeblastdb -in $file -dbtype nucl -parse_seqids -out ${name}_database
blastn -db ${name}_database -query ${blast} -out ${name}_out -outfmt "6 qseqid sseqid evalue bitscore length pident qstart qend qlen sstart send slen" -max_hsps 1
# out format "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
done

for file in $(find ./Blast_transcripts -name "*recip*")
do
name=$(echo $file | sed 's/_blastrecip.*//')
makeblastdb -in ${file} -dbtype nucl -parse_seqids -out ${name}_database_recip
blastn -db ${name}_database_recip -query ${file} -out ${name}_out_recip -outfmt "6 qseqid sseqid evalue bitscore length pident qstart qend qlen sstart send slen" -max_hsps 1
done

###### clean blast results
python3 screen_blast_transcripts.py > transcripts_no_homologous_blast.txt  # transcripts that has no homologous by blast

grep 'not at all' transcripts_no_homologous_blast.txt | sed 's/.*://' | sed 's/ has.*//g' > No_inter-accession_homology_transcripts.txt  # 2505 HOGs 
python3 summarise_screen_result_transcripts.py > non_reciprocal_result_after_screen_transcripts.txt 
grep -o Cang non_reciprocal_result_after_screen_transcripts.txt | wc -l # 5274 Can transcripts
grep -o Colg non_reciprocal_result_after_screen_transcripts.txt | wc -l # 5549 Col transcripts

cat Blast_transcripts/*_homology > reciprocal_transcripts_homology.txt
python3 clean_reciprocal_transcripts_homology.py > reciprocal_transcripts_homology_clean.txt
sed -i 's/AT://' reciprocal_transcripts_homology_clean.txt
sed -i 's/Cang://' reciprocal_transcripts_homology_clean.txt
sed -i 's/Colg://' reciprocal_transcripts_homology_clean.txt
awk '{print $4,$3,$2,$1}' reciprocal_transcripts_homology_clean.txt > reciprocal_transcripts_homology_clean1.txt # 
sed -i 's/^ //g' reciprocal_transcripts_homology_clean1.txt
sed -i 's/ /\t/g' reciprocal_transcripts_homology_clean1.txt 
wc -l reciprocal_transcripts_homology_clean1.txt  # 17341 transcripts
# so plus together the reciprocal results and no isoform results, there are 28180 transcripts that are sure for their homology

# for accession specific transcripts
awk '{print $1}' HOG_geneID_arabidopsis_1to1.txt | sort > HOG_gene_1to1
sed 's/_out.*//g' No_inter-accession_homology_transcripts.txt | sort > HOG_transcripts_no_interaccession_homolog # 2505 HOGs
comm -12 HOG_transcripts_no_interaccession_homolog HOG_gene_1to1 | wc -l # 2444 HOGs are from the 1to1 gene homoplogy

####### use protein sequence to do transcript homology
python3 transcript_blast_pep_prep.py  # generate the Transcript_homology_Gene_no_isoform.txt with genes that has no isoform in all accessions, as well all the fasta files in Blast_transcripts_pep
wc -l Transcript_homology_Gene_no_isoform.txt  # 11990 lines/genes has no isoform in all of the three accessions

for file in $(find ./Blast_transcripts_pep -name "*fasta")
do
name=$(echo $file | sed 's/.fasta//')
makeblastdb -in ${file} -dbtype prot -parse_seqids -out ${name}_database
blastp -db ${name}_database -query ${file} -out ${name}_out -outfmt "6 qseqid sseqid evalue bitscore length pident qstart qend qlen sstart send slen" -max_hsps 1
done


python3 screan_transcript_blastp_result.py > transcripts_no_homologous_blastp.txt # produce out_clean file in the blastp folder
python3 summarize_screen_transcript_blastp_result.py > non_reciprocal_result_after_screen_transcripts_pep.txt  # produce _homology file in blastp folder
grep -o Cang non_reciprocal_result_after_screen_transcripts_pep.txt | wc -l #8777 Can transcripts has no reciprocal results
grep -o Colg non_reciprocal_result_after_screen_transcripts_pep.txt | wc -l #9767 Col transcripts has no reciprocal results
cat Blast_transcripts_pep/*_homology > reciprocal_transcripts_pep_homology.txt
python3 clean_reciprocal_transcripts_pep_homology.py > reciprocal_transcripts_pep_homology_clean.txt
sed -i 's/Cang://' reciprocal_transcripts_pep_homology_clean.txt
sed -i 's/Colg://' reciprocal_transcripts_pep_homology_clean.txt
awk -v OFS='\t' '{print $1,$4,$3,$2}' reciprocal_transcripts_pep_homology_clean.txt > reciprocal_transcripts_pep_homology_clean1.txt
awk '$2=="NA"' reciprocal_transcripts_pep_homology_clean1.txt | wc -l  # 1953 transcripts only in Can and tair 
awk '$4=="NA"' reciprocal_transcripts_pep_homology_clean1.txt | wc -l  # 11014 transcripts only in Col and Can
awk '$3=="NA"' reciprocal_transcripts_pep_homology_clean1.txt | wc -l  # 3644 transcripts only in Col and tair
wc -l reciprocal_transcripts_pep_homology_clean1.txt  # 31832
awk '$2!="NA"' reciprocal_transcripts_pep_homology_clean1.txt | wc -l  # 29879 Col transcripts in the reciprocal results
awk '$4!="NA"' reciprocal_transcripts_pep_homology_clean1.txt | wc -l  # 20818 Can transcripts in the reciprocol results
awk '$3!="NA"' reciprocal_transcripts_pep_homology_clean1.txt | wc -l  # 28188 tair transcripts in the reciprocol results
awk '$2!="NA"' Transcript_homology_Gene_no_isoform.txt | wc -l  # 11560 Col
awk '$4!="NA"' Transcript_homology_Gene_no_isoform.txt | wc -l  # 11145 Can
awk '$3!="NA"' Transcript_homology_Gene_no_isoform.txt | wc -l  # 11266 tair
## to summarise: Col has 29879+11560=41439 transcripts have sure homology; Can has 20818+11145=31963 transcripts have sure homology; Tair has 28188+11266=39454 transcripts have sure homology

### NOTE: genes and HOGs that has no inter-accession blastp results at all
grep Attention transcripts_no_homologous_blastp.txt | wc -l  # 205 HOGs has no homology between the genes at all by preotein screen
grep Cang transcripts_no_homologous_blastp.txt > Can_transcript_no_homology_withinHOG.txt
wc -l Can_transcript_no_homology_withinHOG.txt # 2625 Can transcript has no homology in blastp results
grep Colg transcripts_no_homologous_blastp.txt > Col_transcript_no_homology_withinHOG.txt
wc -l Col_transcript_no_homology_withinHOG.txt  # 2679 Col transcript has no homology in blastp results
grep ARATH transcripts_no_homologous_blastp.txt > Tair_transcript_no_homology_withinHOG.txt
wc -l Tair_transcript_no_homology_withinHOG.txt  # 516 Tair transcript has no homology in blastp results

# accession specific transcripts



# genes without alternative splicing event
python3 generate_gene_list_splicing_event.py  # generate Gene_list_*.txt files
wc -l Gene_list_no_isoform_Col.txt # 16233
wc -l Gene_list_no_isoform_Can.txt # 16747
wc -l Gene_list_with_isoform_Col.txt # 12530
wc -l Gene_list_with_isoform_Can.txt # 11785

#################################################################################################################
####################### below codes has been run, but is not used in the final analysis ################
#### for the OMA group file
python3 main.py /tmp/Ziming_analysis/annotation_update/OMA_synteny/Output/OrthologousGroups.txt > OMAGroups_clean.txt
sed -i 's/transcript_id=//g' OMAGroups_clean.txt

python3 main_1_1.py OMAGroups_clean.txt > OMAGroups_Tair_Can_Col1.txt
python3 main_2_1.py OMAGroups_Tair_Can_Col1.txt > OMAGroups_Tair_Can_Col1_clean.txt
sed -i 's/_braker_pasa_liftoff_coding_ID_tag_function_amended1//g' OMAGroups_Tair_Can_Col1_clean.txt
sed -i 's/_braker_pasa_liftoff_coding_ID_homology_final_function_amended1//g' OMAGroups_Tair_Can_Col1_clean.txt

python3 main_2.py OMAGroups_Tair_Can_Col1_clean.txt > OMAGroups_Tair_Can_Col1_clean1.txt
# check whether there are two proteins asoigned to the the protein in different lines
awk '{print $2}' OrthologousGroups_Tair_Can_Col1_clean1.txt | sort |uniq -d | less #only NA shows up
awk '{print $3}' OrthologousGroups_Tair_Can_Col1_clean1.txt | sort | uniq -d | less # only NA
awk '{print $4}' OrthologousGroups_Tair_Can_Col1_clean1.txt | sort | uniq -d | less #only NA
# check how many transcripts are found homologous between Col and Can
awk '{print $3,$4}' OrthologousGroups_Tair_Can_Col1_clean1.txt | grep -wv NA | wc -l # 35850
# check how many transcripts are accession specific
awk '$2=="ARATa:NA" && $4=="Col:NA"' OrthologousGroups_Tair_Can_Col1_clean1.txt | wc -l # 3685 specific in Can
awk '$2=="ARATH:NA" && $3=="Can:NA"' OrthologousGroups_Tair_Can_Col1_clean1.txt | wc -l # 4129 specific in Col
awk '$$4=="Col:NA" && $3=="Can:NA"' OrthologousGroups_Tair_Can_Col1_clean1.txt | wc -l # 246 specific for Tair10
# check how many transcripts are in the pan-transcriptome file
awk '{print $2}' OrthologousGroups_Tair_Can_Col1_clean1.txt | grep -wv NA | wc -l # 26809 for Tair10
awk '{print $3}' OrthologousGroups_Tair_Can_Col1_clean1.txt | grep -wv NA | wc -l # 39787 for Can
awk '{print $4}' OrthologousGroups_Tair_Can_Col1_clean1.txt | grep -wv NA | wc -l # 41657 for Col

sed '1 i\PanNumber	Tair10_ID	Can_ID	Col_ID' OrthologousGroups_Tair_Can_Col1_clean1.txt > OrthologousGroups_Tair_Can_Col_annotation.txt
sed -i 's/ARATH://g' OrthologousGroups_Tair_Can_Col_annotation.txt
sed -i 's/Can://g' OrthologousGroups_Tair_Can_Col_annotation.txt
sed -i 's/Col://g' OrthologousGroups_Tair_Can_Col_annotation.txt
sed -i 's/t//g' OrthologousGroups_Tair_Can_Col_annotation.txt

awk 'BEGIN{FS=OFS="\t"} {print $3,$1,$4,$2}' OrthologousGroups_Tair_Can_Col_annotation.txt | awk '$1 != "NA"' > OrthologousGroups_Tair_Can_Col_annotation_Can.txt
awk 'BEGIN{FS=OFS="\t"} {print $4,$1,$3,$2}' OrthologousGroups_Tair_Can_Col_annotation.txt | awk '$1 != "NA"' > OrthologousGroups_Tair_Can_Col_annotation_Col.txt
#####


