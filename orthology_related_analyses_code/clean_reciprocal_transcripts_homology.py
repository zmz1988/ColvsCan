import csv
import re



def into_data(row):  # each row is a list
    output = []
    if len(row) < 2:
        pass
    if len(row) >= 2:
        row.sort()
        if not row[0].startswith("Cang"):
            row.insert(0, "Cang:NA")
        if not row[1].startswith("Colg"):
            row.insert(1, "Colg:NA")
    output.append('\t'.join(row))
    return "".join(output)


with open("reciprocal_transcripts_homology.txt", 'r') as fi:
    reader = csv.reader(fi, delimiter="\t")
    for row in reader:
        print(into_data(row))
