# Takes in pickaxe .tsv outputs and finds operators that match a given reaction

import pandas as pd
from rdkit import Chem


def substructure_match(smiles_string, patt):
    """Find pattern (patt) in structure identified by SMILES string,
    and return True if found, False otherwise"""
    comp_struct = Chem.MolFromSmiles(smiles_string)
    substruct = Chem.MolFromSmarts(patt)
    patt_found = comp_struct.HasSubstructMatch(substruct)
    return patt_found


def add_matches(df):
    matches_list = []
    i = 0
    for n, df_row in df.iterrows():
        current_smiles = df_row['SMILES']
        matches_list.append(substructure_match(current_smiles, pattern))
        i += 1
    df['Match?'] = matches_list
    return df


if __name__ == '__main__':
    compound_filepath = 'C:\\Users\Jon\PycharmProjects\Plus1Pathway' \
                        '\Hydroxyacid Enzyme ' \
                        'Discovery\Predicted_Compound_List.tsv '
    compound_dataframe = pd.DataFrame.from_csv(compound_filepath, sep='\t')
    pattern = '[#6h1](-[#6](=[#8])-[#8D1])(-[#8h1])(-[#6h2]-[#6h2])'
    compound_dataframe = add_matches(compound_dataframe)
    matches = 0
    num_rows = 0
    match_list = []
    for index, row in compound_dataframe.iterrows():
        num_rows += 1
        if row['Match?']:
            matches += 1
            match_list.append(index)
            print(index, row['SMILES'])
    print('Matches:', matches, 'out of', num_rows, 'compounds')
    reaction_filepath = 'C:\\Users\Jon\PycharmProjects\Plus1Pathway' \
                        '\Hydroxyacid Enzyme Discovery\Predicted_Reactions.tsv'
    starting_comps = ['C00109', 'C00902', 'C01989', 'C02123', 'C02488',
                      'C05994', 'C06255']
    reaction_dataframe = pd.DataFrame.from_csv(reaction_filepath, sep='\t')
    intermediates = []
    EC_numbers = []
    for comp_id in match_list:
        for index, row in reaction_dataframe.iterrows():
            arrow_index = row['ID Equation'].index('=>')
            if comp_id in row['ID Equation'] and row['ID Equation'].index(
                    comp_id) > arrow_index:
                one_step = False
                for starting_comp_id in starting_comps:
                    if starting_comp_id in row['ID Equation'] and \
                                    row['ID Equation'].index(starting_comp_id) \
                                    < arrow_index:
                        print('one_step: ', starting_comp_id, comp_id,
                              '\t\t\t', row['Operators'])
                        EC_numbers.append(row['Operators'])
                        one_step = True
                if not one_step:
                    current_int_index = row['ID Equation'].index('pk_cpd')
                    current_int = row['ID Equation'][
                                  current_int_index:current_int_index + 13]
                    for index2, row2 in reaction_dataframe.iterrows():
                        arrow_index = row2['ID Equation'].index('=>')
                        if current_int in row2['ID Equation'] and row2['ID ' \
                                                                       'Equation'].index(current_int) > arrow_index:
                            one_step = False
                            for starting_comp_id in starting_comps:
                                if starting_comp_id in row2['ID Equation'] and \
                                                row2['ID Equation'].index(
                                                    starting_comp_id) < arrow_index:
                                    print('two_step: ', starting_comp_id,
                                          current_int, comp_id, '\t\t',
                                          row['Operators'], '\t',
                                          row2['Operators'])
                                    EC_numbers.append((row['Operators'],
                                                       row2['Operators']))
                                    one_step = True
    print(EC_numbers)

