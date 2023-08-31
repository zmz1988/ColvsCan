import re


def single_species_duplication_event(row):
    output = []
    if len(re.findall("Cang", row)) > 1 or len(re.findall("Colg", row)) > 1 or len(re.findall("AT", row)) > 1:
        if len(re.findall("Colg", row)) <= 1 and len(re.findall("AT", row)) <= 1:  # only Can has duplications
            output.append(row)
        elif len(re.findall("Cang", row)) <= 1 and len(re.findall("AT", row)) <= 1:  # only Col has duplications
            output.append(row)  # generate a manual curating list
        elif len(re.findall("Colg", row)) <= 1 and len(re.findall("Cang", row)) <= 1:  # only Tair has duplications
            output.append(row)  # generate a manual curating list
    return "".join(output)


input_file = "./HOG_geneID_arabidopsis.txt"
with open(input_file, 'r') as a_file:
    with open('HOG_genbeID_arabidopsis_single_species_duplication_event.txt', 'w') as f:
        for row in a_file:
            f.write(single_species_duplication_event(row))




