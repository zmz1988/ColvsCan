import re
import pandas as pd

input_file = "/tmp/Ziming_analysis/annotation_update/OMA_analysis_update/HOG_geneID_arabidopsis.txt"


def count_can_col_tair(a_file):
    output = []
    for line in a_file.splitlines():
        if 'Cang' in line and 'Colg' in line and 'AT' in line:
            can_count = line.count("Cang")
            col_count = line.count("Colg")
            tair_count = line.count("AT")
            output.append(
                {
                    'Can': can_count,
                    'Col': col_count,
                    'Tair': tair_count
                }
            )
    df = pd.DataFrame(output)
    can_count_total = df['Can'].sum()
    col_count_total = df['Col'].sum()
    tair_count_total = df['Tair'].sum()
    print("Triple appearance: Can " + str(can_count_total) + "; Col " + str(col_count_total) + "; Tair " + str(tair_count_total))
    return count_can_col_tair


def count_can_col(a_file):
    output = []
    for line in a_file.splitlines():
        if 'Cang' in line and 'Colg' in line and 'AT' not in line:
            can_count = line.count("Cang")
            col_count = line.count("Colg")
            output.append(
                {
                    'Can': can_count,
                    'Col': col_count,
                }
            )
    df = pd.DataFrame(output)
    can_count_total = df['Can'].sum()
    col_count_total = df['Col'].sum()
    print("Pair Col&Can appearance: Can " + str(can_count_total) + "; Col " + str(col_count_total))
    return count_can_col


def count_can_tair(a_file):
    output = []
    for line in a_file.splitlines():
        if 'Cang' in line and 'AT' in line and 'Colg' not in line:
            can_count = line.count("Cang")
            tair_count = line.count("AT")
            output.append(
                {
                    'Can': can_count,
                    'Tair': tair_count
                }
            )
    df = pd.DataFrame(output)
    can_count_total = df['Can'].sum()
    tair_count_total = df['Tair'].sum()
    print("Pair Can&Tair appearance: Can " + str(can_count_total) + "; Tair " + str(tair_count_total))
    return count_can_tair


def count_col_tair(a_file):
    output = []
    for line in a_file.splitlines():
        if 'Cang' not in line and 'AT' in line and 'Colg' in line:
            col_count = line.count("Colg")
            tair_count = line.count("AT")
            output.append(
                {
                    'Col': col_count,
                    'Tair': tair_count
                }
            )
    df = pd.DataFrame(output)
    col_count_total = df['Col'].sum()
    tair_count_total = df['Tair'].sum()
    print("Pair Col&Tair appearance: Col " + str(col_count_total) + "; Tair " + str(tair_count_total))
    return count_col_tair


def count_can_specific(a_file):
    output = []
    for line in a_file.splitlines():
        if 'Cang' in line and 'AT' not in line and 'Colg' not in line:
            can_count = line.count("Cang")
            output.append({'Can': can_count})
    df = pd.DataFrame(output)
    can_count_total = df['Can'].sum()
    print("Can specific genes: " + str(can_count_total))
    return count_can_specific


def count_col_specific(a_file):
    output = []
    for line in a_file.splitlines():
        if 'Cang' not in line and 'AT' not in line and 'Colg' in line:
            col_count = line.count("Colg")
            output.append({'Col': col_count})
    df = pd.DataFrame(output)
    col_count_total = df['Col'].sum()
    print("Col specific genes: " + str(col_count_total))
    return count_col_specific


def count_tair_specific(a_file):
    output = []
    for line in a_file.splitlines():
        if 'Cang' not in line and 'Colg' not in line and 'AT' in line:
            tair_count = line.count("AT")
            output.append({'Tair': tair_count})
    df = pd.DataFrame(output)
    tair_count_total = df['Tair'].sum()
    print("Tair specific genes: " + str(tair_count_total))
    return count_tair_specific


with open(input_file, 'r') as fin:
    text = fin.read()
    can_gene_number = text.count("Cang")
    col_gene_number = text.count("Colg")
    tair_gene_number = text.count("AT")
    print("Total Can gene number: " + str(can_gene_number))
    print("Total Col gene number: " + str(col_gene_number))
    print("Total Tair gene number: " + str(tair_gene_number))
    count_can_col_tair(text)
    count_can_col(text)
    count_can_tair(text)
    count_col_tair(text)
    count_can_specific(text)
    count_col_specific(text)
    count_tair_specific(text)


