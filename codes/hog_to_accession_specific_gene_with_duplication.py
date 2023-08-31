import csv
import os


def shape_data(row):
    i = 0
    output = []
    lst = row.split()
    hog = lst[0]
    for x in lst[1:]:
        hog_new = hog + "_" + str(i)
        new_line = [hog_new, x]
        output.append('\t'.join(new_line))
        i += 1
    return '\n'.join(output)+'\n'


def into_data(row):  # each row is a list
    output = []
    row.sort()
    if not row[0].startswith("AT"):
        row.insert(0, "AT:NA")
    if not row[1].startswith("Cang"):
        row.insert(1, "Cang:NA")
    if not row[2].startswith("Colg"):
        row.insert(2, "Colg:NA")
    output.append('\t'.join(row))
    return "\n".join(output)


input_file = "/tmp/Ziming_analysis/annotation_update/OMA_analysis_update/HOG_arabidopsis_accession_specific_gene_with_duplication.txt"

with open(input_file, 'r') as a_f, open("HOG_arabidopsis_accession_specific_gene_with_duplication_clean.txt", 'w') as f_b:
    for row in a_f:
        f_b.write(shape_data(row))


with open("HOG_arabidopsis_accession_specific_gene_with_duplication_clean.txt", 'r') as fi:
    reader = csv.reader(fi, delimiter="\t")
    for row in reader:
        print(into_data(row))


os.remove("HOG_arabidopsis_accession_specific_gene_with_duplication_clean.txt")



