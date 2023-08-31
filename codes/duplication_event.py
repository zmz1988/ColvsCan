import re
import sys


def duplication_event(row):
    output = []
    if len(re.findall("Cang", row)) > 1 or len(re.findall("Colg", row)) > 1 or len(re.findall("AT", row)) > 1:
        x = re.sub('\.[0-9]*t', '', row)
        x = re.sub('\.[0-9]*', '', x)
        output.append(x)
    return ''.join(output)



input_file = "./HOG_geneID_arabidopsis.txt"
with open(input_file, 'r') as a_file:
    with open('HOG_genbeID_arabidopsis_duplication_event.txt', 'w') as f:
        for row in a_file:
            f.write(duplication_event(row))





