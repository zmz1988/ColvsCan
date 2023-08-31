import pandas as pd
import os.path
import glob

dir_path = "/tmp/Ziming_analysis/annotation_update/OMA_analysis_update/Blast_transcripts_pep"
blast_files = glob.glob(os.path.join(dir_path, '*_out'))


def filter_result(lines):  # lines is a dataframe
    output = []
    g_lines = lines[lines['evalue'] == lines['evalue'].min()]  # select lines that has the lowest e-value
    pd.set_option('mode.chained_assignment', None)
    g_lines['diff'] = (g_lines['qlen'] - g_lines['slen']).abs()
    # g_lines['diff'] = g_lines['qlen'].sub(g_lines['slen'], axis=0).abs()
    if len(g_lines.index) == 1:  # only one lowest e_values
        output.append(g_lines)
    else:
        if g_lines['bitscore'].max() - g_lines['bitscore'].min() < 5:
            further_lines = g_lines[g_lines['diff'] == g_lines['diff'].min()]  # if e-values are the same,
        # then choose the one with the highest bitscore
            if len(further_lines.index) == 1:  # only one lowest diff
                output.append(further_lines)
            else:
                f_further_lines = further_lines[further_lines['bitscore'] == further_lines['bitscore'].max()]
                output.append(f_further_lines)
        else:
            further_lines = g_lines[g_lines['bitscore'] == g_lines['bitscore'].max()]
            if len(further_lines.index) == 1:  # only one highest bitscore
                output.append(further_lines)
            else:
                f_further_lines = further_lines[further_lines['diff'] == further_lines['diff'].min()]
            # if e-value and bitscores are the same among genes, then choose the gene length that are most similar
                output.append(f_further_lines)
    return pd.concat(output)


def tidy_result_self(gene_lst, df, name):
    output = []
    for gen in gene_lst:
        res = df[df['qseqid'] == gen]  # select one query gene blast results
        accession = gen[:2]
        lines = res[~res.sseqid.str.startswith(accession)]  # exclude self-blast results
        if len(lines.index) == 0:
            print(name)
            print(gen)  # print out genes that find no inter-accession homologs
        else:
            col_lines = lines[lines.sseqid.str.startswith("Colg")]
            can_lines = lines[lines.sseqid.str.startswith("Cang")]
            tair_lines = lines[lines.sseqid.str.startswith("ARATH")]
            output.append(filter_result(col_lines))
            output.append(filter_result(can_lines))
            output.append(filter_result(tair_lines))
    if len(output) != 0:  # the 1st if result in empty output
        return pd.concat(output)  # one gene one target


for f in blast_files:
    write_out_name = str(f) + "_clean"
    name = os.path.basename(f)
    with open(f, 'r') as blast_r:
        dd = pd.read_csv(blast_r, sep='\t', names=["qseqid", "sseqid", "evalue", "bitscore", "length", "pident",
                                                   "qstart", "qend", "qlen", "sstart", "send", "slen"])
        dd.sort_values(['qseqid', 'sseqid'], inplace=True)  # need to sort the results, as the order of the genes are different
        gene_lst = dd['qseqid'].unique()  # select query genes
        result = tidy_result_self(gene_lst, dd, name)
        if result is not None:
            result.to_csv(write_out_name, index=False, header=True, sep='\t')
        else:
            print("********************")
            print("Attention:" + name + " has not at all inter-accession homologs \n")
            print("********************")









