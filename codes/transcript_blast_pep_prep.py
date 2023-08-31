import os
from Bio import SeqIO
import re

directory = "Blast_transcripts_pep"
parent_dir = "/tmp/Ziming_analysis/annotation_update/OMA_analysis_update/"
path = os.path.join(parent_dir, directory)
os.mkdir(path)

col_trans = "/tmp/Ziming_analysis/annotation_update/OMA_analysis_update/Col_new_annotation_clean1_amended1.splice"
can_trans = "/tmp/Ziming_analysis/annotation_update/OMA_analysis_update/Can_new_annotation_clean1_amended1.splice"
tair_trans = "/tmp/Ziming_analysis/annotation_update/OMA_synteny/DB/ARATH.splice"
col_pep_fa = "/tmp/Ziming_analysis/annotation_update/OMA_analysis/Col_new_annotation_clean1.pep"
can_pep_fa = "/tmp/Ziming_analysis/annotation_update/OMA_analysis/Can_new_annotation_clean1.pep"
gene_id_file = "/tmp/Ziming_analysis/annotation_update/OMA_analysis_update/HOG_gene_1to1_plus_reciprocal_duplicated_annotation.txt"
tair_pep_fa = "/tmp/Ziming_analysis/annotation_update/OMA_synteny/DB/ARATH.fa"


def build_dictionary_4pep_tair(file):  # file is the protein fasta file
    output = []
    output1 = []
    for row in file:
        if re.search('>', row):
            transcript = row.split(' | ')[0].strip('>')
            output.append(transcript)
            gene = row.split(' | ')[1]
            output1.append(gene)
    pep_tair_dict = dict(zip(output1, output))
    return pep_tair_dict


def build_dictionary_4pep_can_col(file):
    output = []
    output1 = []
    for row in file:
        if re.search('>', row):
            transcript = row.split(' ')[0].strip('>')
            output.append(transcript)
            gene = row.split(' ')[1].split('=')[1]
            output1.append(gene)
    pep_tair_dict = dict(zip(output1, output))
    return pep_tair_dict


def select_transcript(gene, pep_dictionary):   # select the correct name of transcripts from fasta file
    if gene == 'NA':
        trans = 'NA'
    else:
        trans = pep_dictionary[gene]
    return trans


def generate_transcript_list(string, b_f):  # b_f is the _splice.read file
    for line in b_f.splitlines():
        if string in line:
            ge_lst = line.split(';')
            if 'Colg' in ge_lst[0] or 'Cang' in ge_lst[0]:
                new_ge_lst = [s.replace("t", "") for s in ge_lst]
                return new_ge_lst
            else:
                return ge_lst


def extract_fasta(gene, fasta_file, name):
    file_name = str(name + ".fasta")
    path_name = os.path.join(path, file_name)
    with open(path_name, "w") as f:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            if gene == record.id:
                f.write(">" + str(record.id) + "\n")
                f.write(str(record.seq) + "\n")
    return extract_fasta


def add_additional_fasta(gene, fasta_file, name):
    file_name = str(name + ".fasta")
    path_name = os.path.join(path, file_name)
    with open(path_name, "a") as f:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            if gene == record.id:
                f.write(">" + str(record.id) + "\n")
                f.write(str(record.seq) + "\n")
    return extract_fasta


with open(gene_id_file, 'r') as ge_homo, \
        open(col_trans, 'r') as col_splice, \
        open(can_trans, 'r') as can_splice, \
        open(tair_trans, 'r') as tair_splice, \
        open(col_pep_fa, 'r') as col_pep, \
        open(can_pep_fa, 'r') as can_pep, \
        open(tair_pep_fa, 'r') as tair_pep, \
        open("Transcript_homology_Gene_no_isoform.txt", 'w') as no_iso:
    can_splice_db = can_splice.read()
    col_splice_db = col_splice.read()
    tair_splice_db = tair_splice.read()
    tair_pep_dict = build_dictionary_4pep_tair(tair_pep)
    col_pep_dict = build_dictionary_4pep_can_col(col_pep)
    can_pep_dict = build_dictionary_4pep_can_col(can_pep)
    next(ge_homo)
    for row in ge_homo:
        hog, col, can, tair = row.split()
        if len(re.findall('NA', row)) < 2:
            name = hog + "_blast"
            can_in_ca = can in can_splice_db
            col_in_co = col in col_splice_db
            tair_in_tair = col in tair_splice_db
            if not can_in_ca and not col_in_co and not tair_in_tair:
                can_transcript = select_transcript(can, can_pep_dict)  # for the genes that has no isoform at all
                col_transcript = select_transcript(col, col_pep_dict)
                tair_transcript = select_transcript(tair, tair_pep_dict)
                new_line = hog + '\t' + col_transcript + '\t' + can_transcript + '\t' + tair_transcript + '\n'
                no_iso.write(new_line)
            elif can_in_ca and not col_in_co and not tair_in_tair: # only can
                can_trans_ls = generate_transcript_list(can, can_splice_db)
                col_transcript = select_transcript(col, col_pep_dict)
                tair_transcript = select_transcript(tair, tair_pep_dict)
                extract_fasta(col_transcript, col_pep_fa, name)
                add_additional_fasta(tair_transcript, tair_pep_fa, name)
                for ge in can_trans_ls:
                    add_additional_fasta(ge, can_pep_fa, name)
            elif col_in_co and not can_in_ca and not tair_in_tair:  # only col
                col_trans_ls = generate_transcript_list(col, col_splice_db)
                can_transcript = select_transcript(can, can_pep_dict)
                tair_transcript = select_transcript(tair, tair_pep_dict)
                extract_fasta(can_transcript, can_pep_fa, name)
                add_additional_fasta(tair_transcript, tair_pep_fa, name)
                for ge in col_trans_ls:
                    add_additional_fasta(ge, col_pep_fa, name)
            elif tair_in_tair and not can_in_ca and not col_in_co:  # only tair
                tair_trans_ls = generate_transcript_list(tair, tair_splice_db)
                can_transcript = select_transcript(can, can_pep_dict)
                col_transcript = select_transcript(col, col_pep_dict)
                extract_fasta(can_transcript, can_pep_fa, name)
                add_additional_fasta(col_transcript, col_pep_fa, name)
                for ge in tair_trans_ls:
                    add_additional_fasta(ge, tair_pep_fa, name)
            elif tair_in_tair and can_in_ca and not col_in_co:   # only col not
                tair_trans_ls = generate_transcript_list(tair, tair_splice_db)
                can_trans_ls = generate_transcript_list(can, can_splice_db)
                col_transcript = select_transcript(col, col_pep_dict)
                extract_fasta(col_transcript, col_pep_fa, name)
                for ge1 in can_trans_ls:
                    add_additional_fasta(ge1, can_pep_fa, name)
                for ge2 in tair_trans_ls:
                    add_additional_fasta(ge2, tair_pep_fa, name)
            elif col_in_co and can_in_ca and not tair_in_tair:  # only tair not
                col_trans_ls = generate_transcript_list(col, col_splice_db)
                can_trans_ls = generate_transcript_list(can, can_splice_db)
                tair_transcript = select_transcript(tair, tair_pep_dict)
                extract_fasta(tair_transcript, tair_pep_fa, name)
                for ge1 in can_trans_ls:
                    add_additional_fasta(ge1, can_pep_fa, name)
                for ge2 in col_trans_ls:
                    add_additional_fasta(ge2, col_pep_fa, name)
            elif tair_in_tair and col_in_co and not can_in_ca:  # only can not
                tair_trans_ls = generate_transcript_list(tair, tair_splice_db)
                col_trans_ls = generate_transcript_list(col, col_splice_db)
                can_transcript = select_transcript(can, can_pep_dict)
                extract_fasta(can_transcript, can_pep_fa, name)
                for ge1 in col_trans_ls:
                    add_additional_fasta(ge1, col_pep_fa, name)
                for ge2 in tair_trans_ls:
                    add_additional_fasta(ge2, tair_pep_fa, name)
            else:
                col_trans_ls = generate_transcript_list(col, col_splice_db)
                can_trans_ls = generate_transcript_list(can, can_splice_db)
                tair_trans_ls = generate_transcript_list(tair, tair_splice_db)
                for ge1 in col_trans_ls:
                    add_additional_fasta(ge1, col_pep_fa, name)
                for ge2 in tair_trans_ls:
                    add_additional_fasta(ge2, tair_pep_fa, name)
                for ge3 in can_trans_ls:
                    add_additional_fasta(ge3, can_pep_fa, name)




