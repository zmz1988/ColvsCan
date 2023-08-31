import pandas as pd


def find_identical_lines(df):
    v = df['Col'].value_counts()
    v_v = df[df['Col'].isin(v.index[v.gt(1)])]  # collect lines that col value appear more than once
    w = df['Can'].value_counts()
    w_w = df[df['Can'].isin(w.index[w.gt(1)])]
    z = df['Tair'].value_counts()
    z_z = df[df['Tair'].isin(z.index[z.gt(1)])]
    merged1 = v_v.merge(w_w, how='outer', on=["HOG", "Col", "Can", "Tair"], indicator=False)
    merged3 = merged1.merge(z_z, how='outer', on=["HOG", "Col", "Can", "Tair"], indicator=False)
    return merged3


input_file = "/tmp/Ziming_analysis/annotation_update/OMA_analysis_update/HOG_gene_1to1_plus_reciprocal_duplicated.txt"

with open(input_file, 'r') as a_f:
    df = pd.read_csv(a_f, sep='\t', names=["HOG", "Col", "Can", "Tair"])
    merged_trio = find_identical_lines(df)  # get the lines with identical genes
    df_merge = df.merge(merged_trio, how='left', on=["HOG", "Col", "Can", "Tair"], indicator=True)  # get the df without any identical genes
    df_anno = df_merge[df_merge['_merge'] == 'left_only'].iloc[:, 0:4]
    df_anno.to_csv("HOG_gene_1to1_plus_reciprocal_duplicated_annotation.txt", index=False, header=True, sep='\t', na_rep='NA')
    # to generate list of reciprocal duplicated(identical) genes for each accession together with their hog number
    col_gene = merged_trio[["HOG", "Col"]].drop_duplicates('Col')
    can_gene = merged_trio[["HOG", "Can"]].drop_duplicates('Can')
    tair_gene = merged_trio[["HOG", "Tair"]].drop_duplicates('Tair')
    col_gene.to_csv("reciprocal_duplicated_identical_gene_Col.txt", index=False, header=True, sep='\t')
    can_gene.to_csv("reciprocal_duplicated_identical_gene_Can.txt", index=False, header=True, sep='\t')
    tair_gene.to_csv("reciprocal_duplicated_identical_gene_Tair.txt", index=False, header=True, sep='\t')
    merged_trio.to_csv("reciprocal_duplicated_identical_gene_all.txt", na_rep='NA', index=False, header=True, sep='\t')



