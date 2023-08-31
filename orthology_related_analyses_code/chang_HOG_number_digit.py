import re
import pandas as pd


def find_max_digit_numb(f):   # f is the input file
    dd = pd.read_csv(f, sep='\t', usecols=['HOG'])
    hog_string = dd['HOG'].to_list()
    hog_num_list = re.findall(r'\d+', ' '.join(hog_string))
    max_digit = max(list(map(int, hog_num_list)))
    return max_digit


input_file = "HOG_gene_1to1_plus_reciprocal_duplicated_annotation.txt"


with open(input_file, 'r') as a_f, open("HOG_gene_1to1_plus_reciprocal_duplicated_annotation_HOGmaxDigit.txt", 'w') as d_f:
    m_digit = len(str(find_max_digit_numb(a_f)))
    a_f.seek(0)
    for row in a_f:
        if len(re.findall(r'\d+', row)) == 0:
            d_f.write(row + '\n')
        else:
            hog = row.split()[0]
            hog_number = re.findall(r'\d+', hog)[0]
            new_hog_numb = hog_number.zfill(m_digit)
            pre_hog = hog.split(hog_number)[0]
            after_hog = hog.split(hog_number)[1]
            new_hog = pre_hog + new_hog_numb + after_hog
            new_row = [new_hog] + row.split()[1:]
            d_f.write('\t'.join(new_row) + '\n')


