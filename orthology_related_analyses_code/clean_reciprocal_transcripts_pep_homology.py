import csv


tair_pep_fa = "/tmp/Ziming_analysis/annotation_update/OMA_synteny/DB/ARATH.fa"


def build_dictionary_4pep_tair(file):  # file is the protein fasta file
    output = []
    output1 = []
    for line in file:
        if line.startswith(">"):
            transcript = line.split(' | ')[0].strip('>')
            output.append(transcript)
            new_transcript = line.split(' | ')[-1].split('=')[1]
            output1.append(new_transcript)
    pep_tair_dict = dict(zip(output, output1))
    return pep_tair_dict


def into_data(row):  # each row is a list
    if not row[1].startswith("ARATH"):
        row.insert(1, "NA")
    if not row[2].startswith("Cang"):
        row.insert(2, "Cang:NA")
    if not len(row) == 4:
        row.append("Colg:NA")
    return row   # row is a list


with open("reciprocal_transcripts_pep_homology.txt", 'r') as fi, open(tair_pep_fa, 'r') as fa:
    tair_dict = build_dictionary_4pep_tair(fa)
    reader = csv.reader(fi, delimiter="\t")
    for row in reader:
        new_row = into_data(row)
        arath = new_row[1]
        if arath != 'NA':
            new_arath = tair_dict[arath].strip()
            new_row[1] = new_arath
        else:
            new_row = new_row
        print('\t'.join(new_row))

