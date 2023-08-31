import pyham
import re
import sys
from IPython.display import IFrame


#  Select a nwk file as a taxonomy reference
nwk_path = "/tmp/Ziming_analysis/annotation_update/OMA_synteny/Output/EstimatedSpeciesTree.nwk"
#  And extract the newick tree as a string
tree_str = pyham.utils.get_newick_string(nwk_path, type="nwk")
# Then you select your favorite orthoXML file
orthoxml_path = "/tmp/Ziming_analysis/annotation_update/OMA_synteny/Output/HierarchicalGroups.orthoxml"
# pyham.Ham is the main object that containes all information and functionalities.
ham_analysis = pyham.Ham(tree_str, orthoxml_path, use_internal_name=False)

### this section are interactive section, operating within python3
# show the ancestral genome name
for ag in ham_analysis.taxonomy.internal_nodes:
    print("\t- {}".format(ag.name))

# You can also get all extant genomes created as a list
list_genome = ham_analysis.get_list_extant_genomes()
for g in list_genome:
    print("\t- {}".format(g.name))

# get HOG
hog33935 = ham_analysis.get_hog_by_id(33935)
# get all descendant genes
desc_hog = hog33935.get_all_descendant_hogs()
print(desc_hog)

# get all descendant genes clustered by species
desc_genes_clustered = hog33935.get_all_descendant_genes_clustered_by_species()
for species, genes in desc_genes_clustered.items():
    print(species, genes)

# get all descendant level (internal node ancestral genomes)
desc_level = hog33935.get_all_descendant_hog_levels()
for genome in desc_level:
    print(genome.name)

# with a list of genes in the HOG list, we could get the ID of each gene by the below example code
gene367641 = ham_analysis.get_gene_by_id(367641)
print(gene367641.get_dict_xref())
gene76553 = ham_analysis.get_gene_by_id(76553)
print(gene76553.get_dict_xref())

#  Select the arabidopsis claid
Arabidopsis = ham_analysis.get_ancestral_genome_by_name("Can_new_annotation_clean1_amended1/Col_new_annotation_clean1_amended1/ARATH")
# the number of genes associated to Arabidopsis claid
print(Arabidopsis.get_number_genes())  # 28003 common genes/HOGs in this ancestrial genome

############### interactive section end

# get iHam tree view for HOG33935
output_name = "HOG{}.html".format(hog33935.hog_id)
ham_analysis.create_iHam(hog=hog33935, outfile=output_name)
IFrame(output_name, width=700, height=350)  # see two sub hogs among arabidopsis

# overview the whole tree
treeprofile = ham_analysis.create_tree_profile(outfile="tree_profile.html")
IFrame("tree_profile.html", width=680, height=480)

# to get all the HOG only for Arabidopsis claid, including sub HOGs (paralogs) without a HOG ID. Note: sub HOGs are not seperated from the root HOG fasta files in the OMA output folder.
# generate a file called HOG_arabidopsis.txt, containing all genes in each HOG, including sub HOGs
original_stdout = sys.stdout
with open('HOG_arabidopsis.txt', 'w') as f:
    sys.stdout = f
    for h, gs in Arabidopsis.get_ancestral_clustering().items():
        print("HOG: {} -> genes: {}".format(h, gs))
    sys.stdout = original_stdout

# extract gene ID from the dictionary project in pyHam for the HOG genes list created in last step
def find_gene_id(row, row_number):
    output = []
    gene_ls = re.findall(r"\(([0-9]+)\)", row)
    if str(re.search("(HOG\([0-9]+\))", row)) not in 'None':
        output.append("HOG" + gene_ls[0])
        for gene in gene_ls[1:len(gene_ls)]:
            gene_id = ham_analysis.get_gene_by_id(gene)
            gene_dict = gene_id.get_dict_xref()
            prot_id = list(gene_dict.values())[1]
            if prot_id.startswith('Colg') or prot_id.startswith('Cang'):
                prot_id2 = prot_id.split(' ')[0]
                output.append(prot_id2)
            if prot_id.find("transcript_id") != -1:
                prot_id2 = prot_id.split("transcript_id=")[1]
                output.append(prot_id2)
    elif str(re.search("(HOG\([0-9]+\))", row)) in 'None':
        output.append("SUBHOG" + str(row_number))
        for gene in gene_ls[0:len(gene_ls)]:
            gene_id = ham_analysis.get_gene_by_id(gene)
            gene_dict = gene_id.get_dict_xref()
            prot_id = list(gene_dict.values())[1]
            if prot_id.startswith('Colg') or prot_id.startswith('Cang'):
                prot_id2 = prot_id.split(' ')[0]
                output.append(prot_id2)
            if prot_id.find("transcript_id") != -1:
                prot_id2 = prot_id.split("transcript_id=")[1]
                output.append(prot_id2)
    return "\t".join(output)


input_file = "./HOG_arabidopsis.txt"
with open(input_file, 'r') as a_file:
    original_stdout = sys.stdout
    with open('HOG_geneID_arabidopsis.txt', 'w') as f:
        sys.stdout = f
        row_number = 1
        for row in a_file:
            print(find_gene_id(row, row_number))
            row_number += 1
        sys.stdout = original_stdout


# to clean up 'HOG_geneID_arabidopsis.txt' and generate the 1:1 homologous gene file
def clean_homology(row):
    output = []
    if len(re.findall("Cang", row)) > 1 or len(re.findall("Colg", row)) > 1 or len(re.findall("AT", row)) > 1:
        pass
    else:
        x = re.sub('.[0-9]*t', '', row)
        x = re.sub('.[0-9]*', '', x)
        output.append(x)
    return ''.join(output)

def clean_homology_trans(row):
    output = []
    if len(re.findall("Cang", row)) > 1 or len(re.findall("Colg", row)) > 1 or len(re.findall("AT", row)) > 1:
        pass
    else:
        output.append(row)
    return ''.join(output)


input_file = "./HOG_geneID_arabidopsis.txt"
with open(input_file, 'r') as a_file:
    original_stdout = sys.stdout
    with open('HOG_geneID_arabidopsis_1to1.txt', 'w') as f:
        for row in a_file:
            sys.stdout = f
            print(clean_homology(row))
            sys.stdout = original_stdout


with open(input_file, 'r') as a_file:
    original_stdout = sys.stdout
    with open('HOG_transcriptID_arabidopsis_1to1.txt', 'w') as f:
        for row in a_file:
            sys.stdout = f
            print(clean_homology_trans(row))
            sys.stdout = original_stdout

