import pandas as pd
import glob
import os
import os.path
import itertools

directory = "Blast_transcripts_pep"
dir_path = "/tmp/Ziming_analysis/annotation_update/OMA_analysis_update/"
path = os.path.join(dir_path, directory)
blast_files = glob.glob(os.path.join(path, '*_clean'))


def reciprocal_merge(da):  # a reciprocal lines only df
    output = []
    for index, row in da.iterrows():
        condition = (da['qseqid'] == row['sseqid']) & (da['sseqid'] == row['qseqid'])
        row_rev = da[condition]  # reverse row
        row_rev_index = da.index[condition].to_list()[0]  # get the index of the reciprocal lines of the row
        index_list = [index, row_rev_index]
        output.append(da.iloc[[min(index_list)]])   # choose the lowest value of the index to stand for each pair
    if len(output) != 0:  # if no reciprocal records, then len() == 0
        return pd.concat(output).drop_duplicates()  # only keep the unique lines
    else:
        return pd.DataFrame()


def remove_subset_list(li):  # list is the output from the simplify_self_result
    for a, b in itertools.combinations(li, 2):
        if set(a) <= set(b):
            if a in li:
                li.remove(a)
        elif set(b) <= set(a):
            if b in li:
                li.remove(b)
        else:
            pass
    return li


def insert_hog(f_1, hogg):  # f_1 is the output of remove_subset_list
    output = []
    numb = 0
    for x in f_1:
        hogg_new = hogg + "_" + str(numb)
        x.insert(0, hogg_new)
        output.append("\t".join(x))
        numb += 1
    return "\n".join(output)+'\n'


def find_non_reciprocal_gene(df1, df2):
    df1 = set(df1.to_numpy().flatten())   # non-reciprocal
    df2 = set(df2.to_numpy().flatten())  # reciprocal
    ge_lst = list(set(df1).difference(df2))  # find the gene in non-reciprocal but not in reciprocal
    return ge_lst


def simplify_self_result(dc, hog):  # df is a dataframe read from the blast result
    output = []
    d = {'qseqid': 'sseqid', 'sseqid': 'qseqid'}
    df_recip = dc.merge(dc.rename(columns=d))  # get all the reciprocal records
    df_merge = dc.merge(df_recip, how='left', on=['qseqid', 'sseqid'], indicator=True)
    df_non_recip = df_merge[df_merge['_merge'] == 'left_only'].iloc[:, 0:1]  # get all the non-reciprocal records
    if len(df_non_recip.index) != 0:
        f_a = find_non_reciprocal_gene(df_non_recip, df_recip)
        if len(f_a) != 0:
            f_a.insert(0, hog)
            print('\t'.join(f_a))  # print out the non-reciprocal lines
    merged_recip = reciprocal_merge(df_recip)  # merge reciprocal lines into one line for each pair
    if len(merged_recip.index) == 0:
        output = []
    else:
        qgene = merged_recip['qseqid'].unique()  # query gene list (should not contain col genes)
        for gen in qgene:  # both can and arath level
            qgene_arr = merged_recip[merged_recip['qseqid'] == gen]  # get reciprocal merged lines for gen
            sgene_arr = qgene_arr['sseqid']  # get sseqid for each qgene
            for ge in sgene_arr.tolist():  # seqid at can & col level for arath query, or col level for can query
                q_sseq = merged_recip[merged_recip['qseqid'] == ge]  # find lines start with sseqid for each qgene
                if len(q_sseq) == 0:  # could be can&col only for tair query, and all the col level for can query
                    g_ls = qgene_arr[qgene_arr['sseqid'] == ge].values.tolist()
                    output.append(list(itertools.chain.from_iterable(g_ls)))
                elif len(q_sseq) >= 1:  # at this level, must be can&col&tair
                    q_sseq_sseq = q_sseq['sseqid'].tolist()
                    for g in q_sseq_sseq:
                        if ((qgene_arr['qseqid'] == gen) & (
                                qgene_arr['sseqid'] == g)).any():  # can and tair match in qgene
                            gene_set = [gen, ge, g]
                            output.append(gene_set)
                        else:  # we pass these lines here, because in the next run when having Col as the query, this line will be printed
                            pass
    return output


for f in blast_files:
    hog = os.path.basename(f).strip("_blast_out_clean")
    name = hog + "_homology"
    write_out_name = os.path.join(path, name)
    with open(f, 'r') as a_f, open(write_out_name, 'w') as b_f:
        dd = pd.read_csv(a_f, sep='\t', usecols=['qseqid', 'sseqid'])
        new = simplify_self_result(dd, hog)
        if len(new) != 0:
            remove_subset_list(new)  # remove the redundent lines
            result = insert_hog(new, hog)  # add Hog value to the first column
            b_f.write(result)

