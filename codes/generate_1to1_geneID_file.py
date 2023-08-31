import re


def clean_homology(row):
    output = []
    if len(re.findall("Cang", row)) > 1 or len(re.findall("Colg", row)) > 1 or len(re.findall("AT", row)) > 1:
        pass
    else:
        x = re.sub('\.[0-9]*t', '', row)
        y = re.sub('\.[0-9]*', '', x)
        output.append(y)
    return ''.join(output)


input_file = "HOG_geneID_arabidopsis.txt"
with open(input_file, 'r') as a_file:
    with open('HOG_geneID_arabidopsis_1to1.txt', 'w') as f:
        for row in a_file:
            f.write(clean_homology(row))

