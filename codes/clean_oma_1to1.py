import csv
import re



def into_data(row):  # each row is a list
    output = []
    if len(row) < 2:
        pass
    if len(row) >= 2:
        row.sort()
        if not row[0].startswith("AT"):
            row.insert(0, "AT:NA")
        if not row[1].startswith("Cang"):
            row.insert(1, "Cang:NA")
        if not row[2].startswith("Colg"):
            row.insert(2, "Colg:NA")
    output.append('\t'.join(row))
    return "".join(output)


with open("HOG_geneID_arabidopsis_1to1.txt", 'r') as fi:
    reader = csv.reader(fi, delimiter="\t")
    for row in reader:
        print(into_data(row))
