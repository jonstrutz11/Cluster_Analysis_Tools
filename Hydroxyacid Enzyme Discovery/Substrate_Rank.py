# This script ranks enzymes based on the similarity of their native substrate
#  to a substrate of interest. Useful for ranking enzyme candidates that all
#  perform a similar function but act on different substrates (e.g.
#  promiscuous enzymes).
#
# The =fp option allows you to choose a fingerprint (or just do all of them)
#   Fingerprints/keys supported include RDK, Morgan, and MACCS (or All)
#   RDK is chosen by default
#
# Example Command Line Usage:
# python Substrate_Rank.py SmartTableFromMetaCyc.csv Ranked_Outfile.csv -fp All
#
# Updated by Jon Strutz on 4/5/2017

import argparse
import pandas as pd
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys


def read_and_format(file):
    df = pd.read_csv(file, delimiter='\t')
    # Rename InChI columns to distinguish between left and right
    df = df.rename(columns={'InChI': 'Left_InChI', 'InChI.1': 'Right_InChI'})

    for index, row in df.iterrows():
        # Extract UniProt ID from web address in field
        try:
            uniprot_string = row['UniProt']
            uniprot_id_index = uniprot_string.find('uniprot/') + 8
            uniprot_id_index_end = uniprot_string.find('\'>')
            uniprot_id = uniprot_string[uniprot_id_index:uniprot_id_index_end]
            row['UniProt'] = uniprot_id
        except AttributeError:
            print("Warning: UniProt IDs are missing for " + row['Names'])
        # Clean up substrate names and format as list within dataframe
        for side in ['Right', 'Left', 'Right_InChI', 'Left_InChI']:
            for string_to_be_deleted in ['<SUP>', '</SUP>', '<i>', '</i>',
                                         '<SUB>', '</SUB>', '\"']:
                row[side] = row[side].replace(string_to_be_deleted, '')
                row[side] = row[side].replace(string_to_be_deleted, '')
            current_reagent = ''
            reactant_list = []
            for i, char in enumerate(row[side]):
                if row[side][i: i + 2] == '//':
                    reactant_list.append(current_reagent.strip(' // '))
                    current_reagent = ''
                current_reagent += char
            reactant_list.append(current_reagent.strip(' // '))
            row[side] = reactant_list
    return df


def calculate_similarity(df, fp):
    similarity_dict = {}
    # Determine which substrate (out of the set in the reaction) is the
    # native substrate that we are looking to compare
    gene_InChI_dict = get_native_substrate(df)
    # Use rdkit to calculate the tanimoto similarity for each gene's substrate
    #   Note: maxPath is set low to prioritize smaller molecules.
    for gene in gene_InChI_dict:
        current_inchi = gene_InChI_dict[gene][0]
        plus1pathway_inchi = gene_InChI_dict[gene][1]
        ms = [Chem.MolFromInchi(current_inchi),
              Chem.MolFromInchi(plus1pathway_inchi)]
        if fp == 'All' or fp == 'all':
            similarity_dict[gene] = []
            fps = [Chem.RDKFingerprint(x, maxPath=5) for x in ms]
            similarity_dict[gene].append(DataStructs.FingerprintSimilarity(
                                         fps[0], fps[1]))
            fps = [AllChem.GetMorganFingerprint(x, 2) for x in ms]
            similarity_dict[gene].append(DataStructs.DiceSimilarity(fps[0],
                                                                    fps[1]))
            fps = [MACCSkeys.GenMACCSKeys(x) for x in ms]
            similarity_dict[gene].append(DataStructs.FingerprintSimilarity(
                                         fps[0], fps[1]))
        elif fp == 'RDK' or fp == 'rdk':
            fps = [Chem.RDKFingerprint(x, maxPath=5) for x in ms]
            similarity_dict[gene] = DataStructs.FingerprintSimilarity(fps[0],
                                                                      fps[1])
        elif fp == 'Morgan' or fp == 'morgan':
            fps = [AllChem.GetMorganFingerprint(x, 2) for x in
                   ms]
            similarity_dict[gene] = DataStructs.DiceSimilarity(fps[0], fps[1])
        elif fp == 'MACCS':
            fps = [MACCSkeys.GenMACCSKeys(x) for x in ms]
            similarity_dict[gene] = DataStructs.FingerprintSimilarity(fps[0],
                                                                      fps[1])
        else:
            print('Please input a valid fingerprint (RDK or Morgan)')
    return similarity_dict


def get_native_substrate(df):
    gene_InChI_dict = {}
    other_side_reagents = []
    same_side_reagents = []
    for index, row in df.iterrows():
        # Determine reagent to look for based on enzyme class (from EC #)
        if '4.1.1.' in row['EC-Number']:
            reagent_to_find = 'CO2'
            same_side = False
            other_side_reagents = ['H+', 'ATP']
            plus1pathway_InChI = 'InChI=1S/C9H16O5/c1-2-3-4-5-6(8(11)12)7' \
                                 '(10)9(13)14/h6-7,10H,2-5H2,1H3,(H,11,12)' \
                                 '(H,13,14)/t6-,7+/m1/s1'
        elif '1.1.1.' in row['EC-Number']:
            reagent_to_find = 'H+'
            same_side = True
            same_side_reagents = ['H+', 'NADH', 'NADPH', 'NAD(P)H', 'CO2']
            plus1pathway_InChI = 'InChI=1S/C8H14O3/c1-2-3-4-5-6-7(9)8(10)11/' \
                                 'h2-6H2,1H3,(H,10,11)'
        else:
            raise ValueError('No reagent could be found based on the given '
                             'EC Number.')
        # Get gene and its native substrate InChI key
        current_gene = row['Names']
        InChI_list = ['']
        for side in ['Left', 'Right']:
            if reagent_to_find in row[side] and same_side:
                reagent_list = row[side]
                for reagent in same_side_reagents:
                    if reagent in reagent_list:
                        InChI_list = row[side + '_InChI']
                        try:
                            del InChI_list[reagent_list.index(reagent)]
                            reagent_list.remove(reagent)
                        except IndexError:
                            print("InChI key not found for " + reagent + " in "
                                  + row['Names'])
                # Assume that native substrate is first in list after
                # deleting all the coreagents because they are occasionally
                # listed first (~90% of the time, native substrate is first,
                # rest of the time it is listed after coreagents like H+,
                # NADH, etc.)
                gene_InChI_dict[current_gene] = [InChI_list[0],
                                                 plus1pathway_InChI]
            elif reagent_to_find in row[side] and not same_side:
                # Get other side
                if side == 'Left':
                    other_side = 'Right'
                else:
                    other_side = 'Left'
                # Get InChI from other side
                InChI_list = row[other_side + '_InChI']
                for reagent in other_side_reagents:
                    if reagent in row[other_side]:
                        del InChI_list[row[other_side].index(reagent)]
                gene_InChI_dict[current_gene] = [InChI_list[0],
                                                 plus1pathway_InChI]
    return gene_InChI_dict


def generate_outfile(filename, similarity_dict, all_fps=False):
    with open(filename, 'w') as output_file:
        if all_fps:
            output_file.write('Gene' + '\t' + 'RDK_Score' + '\t' +
                              'Morgan_Score' + '\t' + 'MACCS_Score' + '\n')
        else:
            output_file.write('Gene' + '\t' + 'Score' + '\n')
        similarity_dict_length = len(similarity_dict)
        for entry in range(similarity_dict_length):
            max_gene = ''
            max_score = 0
            for gene in similarity_dict:
                # Account for the fact that multiple fingerprints exist in
                # each dictionary item if all_fps is true
                if all_fps:
                    score = similarity_dict[gene][0]
                else:
                    score = similarity_dict[gene]
                # Get maximum score
                if score > max_score:
                    max_gene = gene
                    max_score = score
            # Round to three decimal places and use format to keep trailing
            # zeroes
            if all_fps:
                morgan_score, MACCS_score = similarity_dict[max_gene][1:3]
                output_file.write(max_gene + '\t' +
                                  format(round(max_score, 3), '.3f') + '\t' +
                                  format(round(morgan_score, 3), '.3f') + '\t' +
                                  format(round(MACCS_score, 3), '.3f') + '\n')
            else:
                output_file.write(max_gene + '\t' +
                                  format(round(max_score, 3), '.3f') + '\n')
            del similarity_dict[max_gene]


if __name__ == "__main__":
    # Argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs=1, type=str,
                        help='Must be a tab delimited .csv file (i.e. exported '
                             'from MetaCyc). The fields should be the '
                             'following: Protein, UniProt, Names (gene name), '
                             'EC-Number, Left, InChI, Right, InChI')
    parser.add_argument('outfile', nargs=1, type=str,
                        help='Must be a text file.')
    parser.add_argument('-fp', nargs=1, type=str,
                        help='Defaults to RDK fingerprint. However, [Morgan] '\
                             'and [MACCS] fingerprints can be used as well.')
    args = parser.parse_args()
    # Read and format files with pandas
    infile = vars(args)['infile'][0]
    dataFrame = read_and_format(infile)
    print(dataFrame.head())
    # Calculate similarity for each gene's native substrate based on fingerprint
    try:
        fingerprint = vars(args)['fp'][0]
    except TypeError or KeyError:
        print("RDKFingerprint used by default.")
        fingerprint = 'RDK'
    gene_score_dict = calculate_similarity(dataFrame, fingerprint)
    # Rank and format output text file
    outfile = vars(args)['outfile'][0]
    if fingerprint == 'All' or fingerprint == 'all':
        generate_outfile(outfile, gene_score_dict, all_fps=True)
    else:
        generate_outfile(outfile, gene_score_dict)
