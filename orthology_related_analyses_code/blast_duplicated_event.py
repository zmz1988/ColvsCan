import re
import os
from Bio import SeqIO


directory = "Blast_duplicated_gene"
parent_dir = "/tmp/Ziming_analysis/annotation_update/OMA_analysis_update/"
path = os.path.join(parent_dir, directory)
os.mkdir(path)

col_gene_fa = "/tmp/Ziming_analysis/annotation_update/OMA_analysis/Col_new_annotation_clean1.gene"
can_gene_fa = "/tmp/Ziming_analysis/annotation_update/OMA_analysis/Can_new_annotation_clean1.gene"
tair_gene_fa = "/tmp/Ziming_analysis/Emsemble_Tair10/Arabidopsis_thaliana.TAIR10.55.gene"
Duplication_gene_ID_file = "./HOG_genbeID_arabidopsis_duplication_event.txt"


def extract_fasta(transcript_list, fasta_file, name):
    file_name = name + ".fasta"
    path_name = os.path.join(path, file_name)
    with open(path_name, "w") as f:
        for item in transcript_list:
            for record in SeqIO.parse(fasta_file, 'fasta'):
                if item == record.id:
                    f.write(">" + str(record.id) + "\n")
                    f.write(str(record.seq) + "\n")
    return extract_fasta


def add_additional_fasta(transcript_list, fasta_file, name):
    file_name = name + ".fasta"
    path_name = os.path.join(path, file_name)
    with open(path_name, "a") as f:
        for item in transcript_list:
            for record in SeqIO.parse(fasta_file, 'fasta'):
                if item == record.id:
                    f.write(">" + str(record.id) + "\n")
                    f.write(str(record.seq) + "\n")
    return extract_fasta


with open(Duplication_gene_ID_file, 'r') as fin, \
        open ("HOG_arabidopsis_accession_specific_gene_with_duplication.txt", 'w') as fi:
    for row in fin:
        can_numb = len(re.findall("Cang", row))
        col_numb = len(re.findall("Colg", row))
        tair_numb = len(re.findall("AT", row))
        hog = row.split()[0]
        can_lst = re.findall("Cang[0-9]*", row)
        col_lst = re.findall("Colg[0-9]*", row)
        tair_lst = re.findall("AT[0-9A-Z_-]*", row)
        name = hog + "_db"
        blast_name = hog + "_blast"
        if col_numb == 0 and tair_numb == 0: # can specific
            fi.write(row)
        elif can_numb == 0 and tair_numb == 0: # col specific
            fi.write(row)
        elif can_numb == 0 and col_numb == 0: # tair specific
            fi.write(row)
        elif can_numb > 1 and col_numb <= 1 and tair_numb <= 1: # only Can has duplication
            extract_fasta(can_lst, can_gene_fa, name)
            extract_fasta(col_lst, col_gene_fa, blast_name)
            add_additional_fasta(tair_lst, tair_gene_fa, blast_name)
        elif can_numb <= 1 and col_numb > 1 and tair_numb <= 1: # only Col has duplication
            extract_fasta(col_lst, col_gene_fa, name)
            extract_fasta(can_lst, can_gene_fa, blast_name)
            add_additional_fasta(tair_lst, tair_gene_fa, blast_name)
        elif can_numb <= 1 and col_numb <= 1 and tair_numb > 1: # only Tair has duplication
            extract_fasta(tair_lst, tair_gene_fa, name)
            extract_fasta(can_lst, can_gene_fa, blast_name)
            add_additional_fasta(col_lst, col_gene_fa, blast_name)
        else: # multi-accession duplication
            self_name = str(hog + '_self')
            extract_fasta(can_lst, can_gene_fa, self_name)
            add_additional_fasta(col_lst, col_gene_fa, self_name)
            add_additional_fasta(tair_lst, tair_gene_fa, self_name)













