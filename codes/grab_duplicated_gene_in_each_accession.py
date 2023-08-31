import re

input_file1 = "HOG_genbeID_arabidopsis_duplication_event.txt"


with open(input_file1, 'r') as a_f, open("Col_duplicated_gene_list.txt", 'w') as f, \
        open("Can_duplicated_gene_list.txt", 'w') as f_b, open("Tair_duplicated_gene_list.txt", 'w') as f_c:
    for row in a_f:
        col = re.findall("Colg[0-9]*", row)
        can = re.findall("Cang[0-9]*", row)
        tair = re.findall("AT[0-9][0-9]*", row)
        f.write('\n'.join(col)+'\n')
        f_b.write('\n'.join(can)+'\n')
        f_c.write('\n'.join(tair)+'\n')




