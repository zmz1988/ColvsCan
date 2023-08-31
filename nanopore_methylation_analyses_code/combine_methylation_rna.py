import sys
import pandas as pd
import os

methy_file = sys.argv[1]
rna_file = sys.argv[2]

with open(methy_file, 'r') as meth_f, open(rna_file, 'r') as rna_f:
    output_name = os.path.basename(methy_file) + '_' + os.path.basename(rna_file)
    methyl_result = pd.read_csv(meth_f, sep='\t', names=["chr", "start", "end", "strand", "gene_name", "methy_persentage", "Count_locus"])
    rna_result = pd.read_csv(rna_f, sep='\t')
    rna_result.rename(columns={rna_result.columns[0]: "gene_name"}, inplace=True)
    rna_result["Can_mean"] = rna_result.filter(like='Can').mean(axis=1)
    rna_result["Col_mean"] = rna_result.filter(like='Col').mean(axis=1)
    merged_df = methyl_result.merge(rna_result, left_on='gene_name', right_on='gene_name')
    merged_df.to_csv(output_name, sep='\t', header=True, index=False)




