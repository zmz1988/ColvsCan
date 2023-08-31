for file in $(ls ../*.fasta)
do
name=$(basename $file | cut -d"_" -f1)
Subread=$(find /shared/ucl/depts/mottlab/Ziming_analysis/Founder_PacBio_DNA/ -name "${name}*")
ONTreads=$(find /shared/ucl/depts/mottlab/Ziming_analysis/Nanopore_supe_accuracy_reads/ -name "${name}*.fastq.gz")
minimap2 -t 20 -ax map-hifi $file $Subread --secondary=no > ${name}.hifi.sam 
samtools view -Sb ${name}.hifi.sam | samtools sort -@ 20 > ${name}.hifi.bam
samtools index ${name}.hifi.bam
minimap2 -t 20 -ax map-ont $file $ONTreads --secondary=no > ${name}.ont.sam
samtools view -Sb ${name}.ont.sam | samtools sort -@ 20 > ${name}.ont.bam
samtools index ${name}.ont.bam
bedtools genomecov -ibam ${name}.hifi.bam -bga > ${name}.hifi.coverage
bedtools genomecov -ibam ${name}.ont.bam -bga > ${name}.ont.coverage
done



