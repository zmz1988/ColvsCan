#job name
#$ -N bismark
#name of logfile 
#$ -o bismark3.log
# Set the working directory 
#$ -wd /SAN/mottlab/ArabidopsisAssembly/Zhelyazkov/Ziming_analysis/bisulphite_analysis_new
# can also use #$-cwd for current working directory (where job was submitted
#memory and runtime options 
#$ -l tmem=7G
#$ -l h_vmem=7G
#$ -l h_rt=24:00:00 
#$ -S /bin/bash
#$ -j y
#Number of threads requested per job
#$ -pe smp 3 
#Number of jobs requested. Each will be given a variable $SGE_TASK_ID
#$ -t 1-4
#$ -tc 4

echo $SGE_TASK_ID $(date)

PATH=/home/zzhong/bin:$PATH
PATH=/share/apps/mottlab/miniconda3/bin/:$PATH
PATH=/share/apps/mottlab/miniconda3/envs/megalodon/bin/:$PATH

Sample=$(find /SAN/mottlab/ArabidopsisAssembly/Zhelyazkov/Ziming_analysis/bisulphite_analysis_new/data -mindepth 1 -type d | head -${SGE_TASK_ID}  | tail -1)
Maindir=/SAN/mottlab/ArabidopsisAssembly/Zhelyazkov/Ziming_analysis/bisulphite_analysis_new
SampleName=$(echo $Sample | cut -f 9 -d"/")
Reference=$(echo $SampleName | cut -f 1 -d"_")

#mkdir ./new_analysis/${SampleName}
cd ./new_analysis/$SampleName
conda activate megalodon

for file in $(find ${Sample}/ -name "*_1.fq.gz")
do
Repname=$(basename $file | sed 's/_1.fq.gz//')
R1=$(echo $file)
R2=$(echo $R1 | sed 's/_1.fq.gz/_2.fq.gz/')
######## run trim galore to trim off bad quality and adapters

/home/zzhong/TrimGalore-0.6.6/trim_galore --paired --clip_r1 10 --clip_r2 10 --fastqc -j 2 ${R1} ${R2} 2> ${Repname}_trim_galore.log

######## run bismark fro sequence mapping and methylation calling
bismark_genome_preparation --path_to_aligner /share/apps/mottlab/miniconda3/bin/ --bowtie2 --parallel 2 /SAN/mottlab/ArabidopsisAssembly/Zhelyazkov/Ziming_analysis/bisulphite_analysis_new/reference/${Reference}

# 1st bismark mapping
echo "Bismark mapping Start " $(date)

TrimR1=$(basename $R1 | sed 's/.fq.gz/_val_1.fq.gz/')
TrimR2=$(basename $R2 | sed 's/.fq.gz/_val_2.fq.gz/')
Lane=$(basename $R1 | cut -f 4,5 -d"_")

bismark -N 1 --genome /SAN/mottlab/ArabidopsisAssembly/Zhelyazkov/Ziming_analysis/bisulphite_analysis_new/reference/${Reference} -1 ${TrimR1} -2 ${TrimR2} --non_directional --score_min L,0,-0.4 --parallel 3 --rg_tag --rg_id ${Lane} --rg_sample ${SampleName} -un --temp_dir bismark_tmp --ambig_bam --nucleotide_coverage 2> bismark_mapping.log

/home/zzhong/TrimGalore-0.6.6/trim_galore --paired --hardtrim5 70 --fastqc -j 2 ${TrimR1}_unmapped_reads_1.fq.gz ${TrimR2}_unmapped_reads_2.fq.gz

bismark -N 1 --genome /SAN/mottlab/ArabidopsisAssembly/Zhelyazkov/Ziming_analysis/bisulphite_analysis_new/reference/${Reference}  -1 ${TrimR1}_unmapped_reads_1.70bp_5prime.fq.gz -2 ${TrimR2}_unmapped_reads_2.70bp_5prime.fq.gz --non_directional --score_min L,0,-0.4 --parallel 3 --rg_tag --rg_id ${Lane} --rg_sample ${SampleName} --unmapped --temp_dir bismark_tmp --ambig_bam --nucleotide_coverage 2> bismark_mapping_shoreten1.log

done

# bimark deduplicate
Bamfile=$(ls ./*_bismark_bt2_pe.bam | sed 's/$/ /g' | tr -d '\n')
deduplicate_bismark --paired -o ${SampleName} --bam --multiple ${Bamfile} 2> deduplicate_bismark.log

# methylation extraction  
bismark_methylation_extractor -p --no_overlap --parallel 3 --bedGraph --zero_based --CX_context --gzip --cytosine_report   --genome_folder /SAN/mottlab/ArabidopsisAssembly/Zhelyazkov/Ziming_analysis/bisulphite_analysis_new/reference/${Reference} ${SampleName}.multiple.deduplicated.bam 2> methylation_extractor.log

# generating a report for the mapping analysis

bismark2report ./

bismark2summary ./


echo $SGE_TASK_ID $(date)
 
 
