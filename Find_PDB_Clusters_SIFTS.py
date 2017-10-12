# Takes in a ranked list of cluster centroid uniprot ids (from
#  rank_and_plot_distmtx.py) and ranks enzymes within each cluster by number
#  of crystal structures in PDB database.
#

import os
import pandas as pd
import argparse
import shutil


def find_cluster_files(uniprot_id_list, cluster_file_directory):
    id_clusterFile_dict = {}
    for uniprot_id in uniprot_id_list:
        for file in os.listdir(cluster_file_directory):
            # Get filepath for each file
            file_path = os.path.join(cluster_file_directory, file)
            with open(file_path, 'r', encoding='ISO-8859-1') as \
                    current_cluster_infile:
                for line in current_cluster_infile:
                    if line[0] == '>' and uniprot_id in line:
                        # Add path to cluster file that uniprot ID is in
                        id_clusterFile_dict[uniprot_id] = file_path
        if uniprot_id not in id_clusterFile_dict:
            raise Exception("Error - cluster file not found for", uniprot_id)
    return id_clusterFile_dict


def extract_uniprotIDs(filepath):
    """Extracts all the uniprotIDs from a FASTA file as a list"""
    id_list = []
    with open(filepath, 'r') as fasta_infile:
        for line in fasta_infile:
            if line[0] == '>':
                # Get uniprotID from seq description line
                id_list.append(line.split('|')[1])
    return id_list


def map_uniprot_to_PDB(id_list, mapping_df):
    """Finds PDB IDs for every enzyme within a cluster (if available)"""
    PDB_id_dict = {}

    for id in id_list:
        if id in mapping_df.index:
            PDB_id_dict[id] = mapping_df.loc[id]['PDB']

    # return dictionary of form {uniprotID : PDBID}
    return PDB_id_dict


def get_centroid_uniprot_ids(file_path):
    centroid_uniprot_ids = []
    with open(file_path, 'r') as ranked_list_infile:
        # Skip first line
        ranked_list_infile.readline()
        for line in ranked_list_infile:
            # Get first id on each line (not reviewed)
            centroid_uniprot_ids.append(line.split(',')[0].strip())
    return centroid_uniprot_ids


def map_all(centroid_ids, cluster_file_dict, mapping_df):
    PDB_id_dict = {}
    for cluster_id in centroid_ids:
        cluster_filepath = cluster_file_dict[cluster_id]
        cluster_uniprotID_list = extract_uniprotIDs(cluster_filepath)
        # Map all uniprotIDs to PDB
        PDB_data = map_uniprot_to_PDB(cluster_uniprotID_list, mapping_df)
        if PDB_data:
            PDB_id_dict[cluster_id] = PDB_data
        print('Getting PDB info for', cluster_id)
    return PDB_id_dict


if __name__ == '__main__':
    # Get command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs=1, type=str,
                        help='Ranked list (centroids) text file generated from '
                             'rank_and_plot_distmtx.py.')
    parser.add_argument('outfile', nargs=1, type=str,
                        help='Ranked list of centroids with a reviewed sequence'
                             ' id (if available) inserted after the centroid '
                             'sequence id (separated by a comma). File format '
                             'is txt.')
    parser.add_argument('directory', nargs=1, type=str,
                        help='File path to a folder with cluster fasta '
                             'files (filenames should be cluster numbers with '
                             '.fasta extension). These files can be generated '
                             'using the -clusters option during usearch '
                             'clustering.')
    args = parser.parse_args()

    cluster_file_directory = vars(args)['directory'][0]
    centroid_ranked_list_file = vars(args)['infile'][0]
    outfile = vars(args)['outfile'][0]

    centroid_uniprot_ids = get_centroid_uniprot_ids(centroid_ranked_list_file)

    # Find cluster number for each id
    cluster_file_dict = find_cluster_files(centroid_uniprot_ids,
                                           cluster_file_directory)

    mapping_df = pd.read_csv('uniprot_pdb.csv')
    mapping_df = mapping_df.set_index('SP_PRIMARY')

    PDB_id_dict = map_all(centroid_uniprot_ids, cluster_file_dict, mapping_df)

    with open(outfile, 'w') as out_file:
        with open(centroid_ranked_list_file, 'r') as in_file:
            i = 0
            for line in in_file:
                if i == 0:
                    out_file.write(line.rstrip('\n') + ',PDB Structures ('
                                   'UniprotID : PDB ID list\n')
                else:
                    uniprot_id = ''
                    j = 0
                    char = line[j]
                    while char != ',' and char != '\n':
                        uniprot_id += char
                        j += 1
                        char = line[j]
                    if uniprot_id in PDB_id_dict:
                        out_file.write(line.rstrip('\n') + ',' +
                                       str(PDB_id_dict[uniprot_id]) + '\n')
                    else:
                        out_file.write(line)
                i += 1
