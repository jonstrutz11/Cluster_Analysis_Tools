# Warning: This script requires a ranked list text file to be generated from
# rank_and_plot_distmtx.py. It also requires that a folder with fasta files
# for each cluster in the distance matrix is available.
#
# This script takes a ranked list of centroid ids and finds if any reviewed
# sequences exist within the same cluster.
#
# Example:
# python .\find_reviewed.py ranked_list.txt ranked_list_with_reviewed.txt
#   'C:\\Users\\Name\\Project\\Cluster Fasta Files\\'


import os
import argparse


# Get command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('infile', nargs=1, type=str,
                    help='Ranked list (centroids) text file generated from '
                         'rank_and_plot_distmtx.py.')
parser.add_argument('outfile', nargs=1, type=str,
                    help='Ranked list of centroids with a reviewed sequence id '
                         '(if available) inserted after the centroid sequence '
                         'id (separated by a comma). File format is txt.')
parser.add_argument('directory', nargs=1, type=str,
                    help='File path to a folder with cluster fasta '
                         'files (filenames should be cluster numbers with '
                         '.fasta extension). These files can be generated '
                         'using the -clusters option during usearch '
                         'clustering.')
args = parser.parse_args()

directory = vars(args)['directory'][0]


# Get a dictionary relating cluster filenames (e.g. 123.fasta) with one
# reviewed sequence in that cluster, if available.
reviewed_dict = {}
for filename in os.listdir(directory):
    with open(directory + filename) as infile:
        cluster_num = filename[0:-6]
        for line in infile:
            # sp stands for swiss-prot which means the entry has been
            # reviewed on uniprot
            if line[1:3] == 'sp':
                end = 6
                char = line[6]
                while char != '|':
                    char = line[end]
                    end += 1
                current_id = line[4:end - 1]
                reviewed_dict[filename[0:-6]] = current_id


# Get list of ranked clusters from rank_and_plot_distmtx outfile
ranked_list = []
with open(vars(args)['infile'][0]) as infile:
    for line in infile:
        end = 6
        char = line[6]
        while char != '|':
            char = line[end]
            end += 1
        current_id = line[3:end - 1]
        ranked_list.append(current_id)


# Find the cluster that each centroid_id is found in and relate the cluster
# number to that id in ranked_dict.
ranked_dict = {}
for centroid_id in ranked_list:
    for filename in os.listdir(directory):
        # FASTA files end in cluster_num.fasta
        cluster_num = filename[0:-6]
        with open(directory + filename) as infastafile:
            for line in infastafile:
                if line[0] == '>':
                    # Shortest uniprot ids are 6 characters.
                    end = 6
                    char = line[6]
                    while char != '|':
                        char = line[end]
                        end += 1
                    current_id = line[4:end - 1]
                    if centroid_id == current_id:
                        ranked_dict[centroid_id] = cluster_num


# Create output file by writing centroids and, if the centroids cluster
# number is found in the reviewed_dict (which contains cluster numbers for
# reviewed sequences only), then write that reviewed sequence id next to the
# centroid id.
with open(vars(args)['outfile'][0], 'w') as outfile:
    outfile.write('Uniprot ID,Reviewed ID in same cluster (if '
                  'available)\n')
    for centroid in ranked_list:
        outfile.write(centroid)
        # Second condition to make sure that the reviewed sequence id does
        # not match the centroid id.
        if ranked_dict[centroid] in reviewed_dict and \
                reviewed_dict[ranked_dict[centroid]] != ranked_dict[centroid]:
            outfile.write(',' + reviewed_dict[ranked_dict[centroid]] + '\n')
        else:
            outfile.write('\n')
