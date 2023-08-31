import sys
import os
import re
from Bio import SeqIO

directory = "Blast_transcripts"
parent_dir = "/tmp/Ziming_analysis/annotation_update/OMA_analysis_update/"
path = os.path.join(parent_dir, directory)
os.mkdir(path)

Col_trans = "Col_new_annotation_clean1_amended1.splice"
Can_trans = "Can_new_annotation_clean1_amended1.splice"
Col_trans_fa = "/tmp/Ziming_analysis/annotation_update/OMA_analysis/Col_new_annotation_clean1.mRNA"
Can_trans_fa = "/tmp/Ziming_analysis/annotation_update/OMA_analysis/Can_new_annotation_clean1.mRNA"
gene_ID_file = "./HOG_gene_1to1_plus_reciprocal_duplicated_annotation.txt"


def generate_transcript_list(string, b_f):
    for line in b_f.splitlines():
        if string in line:
            return ''.join(line)


def extract_fasta(transcript_list, fasta_file, name):
    file_name = str(name + ".fasta")
    path_name = os.path.join(path, file_name)
    with open(path_name, "w") as f:
        for item in transcript_list.split(";"):
            for record in SeqIO.parse(fasta_file, 'fasta'):
                if item.strip("t") == record.id:
                    f.write(">" + str(record.id) + "\n")
                    f.write(str(record.seq) + "\n")
    return extract_fasta


def add_additional_fasta(transcript_list, fasta_file, name):
    file_name = str(name + ".fasta")
    path_name = os.path.join(path, file_name)
    with open(path_name, "a") as f:
        for item in transcript_list.split(";"):
            for record in SeqIO.parse(fasta_file, 'fasta'):
                if item.strip("t") == record.id:
                    f.write(">" + str(record.id) + "\n")
                    f.write(str(record.seq) + "\n")
    return extract_fasta


with open(gene_ID_file, 'r') as fin, \
        open(Col_trans, 'r') as col_splice, \
        open(Can_trans, 'r') as can_splice, \
        open("HOG_transcriptID_arabidopsis_1to1_noIsof.txt", 'w') as fi:
    can_splice_db = can_splice.read()
    col_splice_db = col_splice.read()
    for row in fin:
        if "Cang" in row and "Colg" in row:
            hog, col, can, tair = row.split()
            name = hog + "_blast"
            db_name = hog + "_db"
            can_transcript = can + ".1t"
            col_transcript = col + ".1t"
            can_in_ca = can in can_splice_db
            col_in_co = col in col_splice_db
            if not can_in_ca and not col_in_co:
                fi.write('\t'.join(row.split()[0:3]) + '\n')  # for the genes that has no isoform at all
            elif can_in_ca and not col_in_co: # genes that Can has isoform but Col not
                can_trans_ls = generate_transcript_list(can, can_splice_db)
                extract_fasta(can_trans_ls.strip(), Can_trans_fa, db_name)
                extract_fasta(col_transcript, Col_trans_fa, name)
            elif not can_in_ca and col_in_co:  # genes that Can has no isoform but Col has
                col_trans_ls = generate_transcript_list(col, col_splice_db)
                extract_fasta(col_trans_ls.strip(), Col_trans_fa, db_name)
                extract_fasta(can_transcript, Can_trans_fa, name)
            else:  # genes that both Can and Col have isoforms
                col_trans_ls = generate_transcript_list(col, col_splice_db)
                can_trans_ls = generate_transcript_list(can, can_splice_db)
                ne_name = name + "recip"
                extract_fasta(can_trans_ls.strip(), Can_trans_fa, ne_name)
                add_additional_fasta(col_trans_ls.strip(), Col_trans_fa, ne_name)


