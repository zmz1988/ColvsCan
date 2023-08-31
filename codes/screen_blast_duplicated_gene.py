import pandas as pd
import glob
import re
import os

directory = "Blast_duplicated_gene"
dir_path = "/tmp/Ziming_analysis/annotation_update/OMA_analysis_update/"
path = os.path.join(dir_path, directory)
blast_files = glob.glob(os.path.join(path, '*_out'))


def find_self_file(f):
    name = f.strip("_out") + "_*.fasta"
    file_ls = glob.glob(os.path.join(path, name))
    return len(file_ls) == 1  # file a_f is a self_blast file


def filter_result(lines):  # lines is a dataframe
    output = []
    g_lines = lines[lines['evalue'] == lines['evalue'].min()]  # select lines that has the lowest e-value
    pd.set_option('mode.chained_assignment', None)
    g_lines['diff'] = (g_lines['qlen'] - g_lines['slen']).abs()
    if len(g_lines.index) == 1:  # only one lowest e_values
        output.append(g_lines)
    else:
        further_lines = g_lines[g_lines['bitscore'] == g_lines['bitscore'].max()]  # if e-values are the same, then choose the one with the highest bitscore
        if len(further_lines.index) == 1:  # only one highest bitscore
            output.append(further_lines)
        else:
            f_further_lines = further_lines[further_lines['diff'] == further_lines['diff'].min()]  # if e-value and bitscores are the same among genes, then choose the gene length that are most similar
            # further_lines = g_lines.loc[((g_lines['diff'] == min_diff) & (g_lines['length'] == max_length))]
            output.append(f_further_lines)
    return pd.concat(output)



def tidy_non_self_result(gene_lst, df):
    output = []
    for gen in gene_lst:
        res = df[df['qseqid'] == gen]  # select one query gene blast results
        output.append(filter_result(res))
    return pd.concat(output)


def tidy_result_self(gene_lst, df, name):
    output = []
    for gen in gene_lst:
        res = df[df['qseqid'] == gen]  # select one query gene blast results
        accession = gen[:2]
        lines = res[~res.sseqid.str.startswith(accession)]  # exclude self-blast results
        if len(lines.index) == 0:
            print(name)
            print(res)  # print out genes that find no inter-accession homologs
        else:
            col_lines = lines[lines.sseqid.str.startswith("Colg")]
            can_lines = lines[lines.sseqid.str.startswith("Cang")]
            tair_lines = lines[lines.sseqid.str.startswith("AT")]
            output.append(filter_result(col_lines))
            output.append(filter_result(can_lines))
            output.append(filter_result(tair_lines))
    if len(output) != 0:
        return pd.concat(output)  # one gene one target


for f in blast_files:
    write_out_name = str(f) + "_clean"
    name = os.path.basename(f)
    if os.stat(f).st_size == 0:
        print("********************")
        print("Attention:" + name + " has not at all inter-accession homologs \n")
        print("********************")
    else:
        with open(f, 'r') as a_f:
            df = pd.read_csv(a_f, sep='\t',
                             names=["qseqid", "sseqid", "evalue", "bitscore", "length", "pident", "qstart", "qend",
                                    "qlen",
                                    "sstart", "send", "slen"])
            gene_lst = df['qseqid'].unique()  # select query genes
            if find_self_file(f):  # if is a self_blast result
                result = tidy_result_self(gene_lst, df, name)
                if result is not None:
                    result.to_csv(write_out_name, index=False, header=True, sep='\t')
                else:
                    print("********************")
                    print("Attention:" + name + " has not at all inter-accession homologs \n")
                    print("********************")
            else:
                result = tidy_non_self_result(gene_lst, df)
                if result is not None:
                    result.to_csv(write_out_name, index=False, header=True, sep='\t')

