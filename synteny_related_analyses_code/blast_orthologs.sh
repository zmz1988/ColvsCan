#conda activate EDTA

#makeblastdb -in /tmp/Ziming_analysis/annotation/Col_final_23022022_new_hapog_patch_hapog_reorder.fasta -dbtype nucl -parse_seqids -out Col_genome_db

#makeblastdb -in /tmp/Ziming_analysis/annotation/Can_final_23022022_new_hapog.fasta -dbtype nucl -parse_seqids -out Can_genome_db

#blastn -db Can_genome_db -query Can_AT1G01370.gene.fasta -out Can_AT1G01370.gene_blast -outfmt 6
#blastn -db Can_genome_db -query Can_AT1G15660.gene.fasta -out Can_AT1G15660.gene_blast -outfmt 6
#blastn -db Can_genome_db -query Can_AT5G55820.gene.fasta -out Can_AT5G55820.gene_blast -outfmt 6
#blastn -db Can_genome_db -query 18s.fasta -out Can_18s_blast -outfmt 6
#blastn -db Can_genome_db -query 25s.fasta -out Can_25s_blast -outfmt 6
#blastn -db Can_genome_db -query 5.8s.fasta -out Can_5.8s_blast -outfmt 6

#blastn -db Col_genome_db -query Col_AT1G01370.gene.fasta -out Col_AT1G01370.gene_blast -outfmt 6
#blastn -db Col_genome_db -query Col_AT1G15660.gene.fasta -out Col_AT1G15660.gene_blast -outfmt 6
#blastn -db Col_genome_db -query Col_AT5G55820.gene.fasta -out Col_AT5G55820.gene_blast -outfmt 6
#blastn -db Col_genome_db -query 18s.fasta -out Col_18s_blast -outfmt 6
#blastn -db Col_genome_db -query 25s.fasta -out Col_25s_blast -outfmt 6
#blastn -db Col_genome_db -query 5.8s.fasta -out Col_5.8s_blast -outfmt 6

blastn -db Can_genome_db -query AT4G39200_S25_chr4.fasta -out AT4G39200_S25_chr4_Can -outfmt 6
blastn -db Can_genome_db -query AT2G21580_S25_chr2.fasta -out AT2G21580_S25_chr4_Can -outfmt 6
blastn -db Col_genome_db -query AT4G39200_S25_chr4.fasta -out AT4G39200_S25_chr4_Col -outfmt 6
blastn -db Col_genome_db -query AT2G21580_S25_chr2.fasta -out AT2G21580_S25_chr2_Col -outfmt 6

