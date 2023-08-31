
col_trans = "/tmp/Ziming_analysis/annotation_update/OMA_analysis_update/Col_new_annotation_clean1_amended1.splice"
can_trans = "/tmp/Ziming_analysis/annotation_update/OMA_analysis_update/Can_new_annotation_clean1_amended1.splice"
col_gene = "/tmp/Ziming_analysis/annotation_update/OMA_analysis/Col_new_annotation_clean1.gene"
can_gene = "/tmp/Ziming_analysis/annotation_update/OMA_analysis/Can_new_annotation_clean1.gene"


def build_gene_lst_col_can(file):  # file is the .gene file
    output = []
    for line in file:
        if line.startswith('>'):
            gene = line.split(' ')[0]. strip('>')
            output.append(gene)
    return output


with open(col_trans, 'r') as col_splice, \
        open(can_trans, 'r') as can_splice, \
        open(col_gene, 'r') as col_g, \
        open(can_gene, 'r') as can_g, \
        open("Gene_list_no_isoform_Col.txt", 'w') as no_iso_col, \
        open("Gene_list_no_isoform_Can.txt", 'w') as no_iso_can, \
        open("Gene_list_with_isoform_Col.txt", 'w') as with_iso_col, \
        open("Gene_list_with_isoform_Can.txt", 'w') as with_iso_can:
    can_splice_db = can_splice.read()
    col_splice_db = col_splice.read()
    can_g_lst = build_gene_lst_col_can(can_g)
    col_g_lst = build_gene_lst_col_can(col_g)
    for x in can_g_lst:
        if x in can_splice_db:
            with_iso_can.write(x + "\n")
        else:
            no_iso_can.write(x + "\n")
    for y in col_g_lst:
        if y in col_splice_db:
            with_iso_col.write(y + "\n")
        else:
            no_iso_col.write(y + "\n")




