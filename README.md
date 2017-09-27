# Cluster_Analysis_Tools
Tools to determine the best small subset of proteins to test for initial activity within a very large set of proteins (e.g. a protein superfamily).

Cluster Analysis Pipeline Steps:
1. Get a FASTA file with the set of sequences of interest (e.g. by searching for a group of proteins on uniprot and downloading a            FASTA file with all of them). It may be useful to filter this set of sequences by length as well (to prevent getting fragments, for        example).
2. Use sequence_cleaner.py to remove duplicate sequences. This can also filter by a minimum length as well.
3. Use usearch to cluster and get centroid sequences. I used the -cluster_fast method with an -id (sequence identity) of 0.7. Typical        ranges for this number are 0.7-0.99 depending on the application. Make sure to also set -sizeout so that the size of each cluster is      annotated. I also recommend setting a minimum cluster size via the -minsize option to reduce the number of very small clusters. Lastly,    make sure to generate a directory of fasta files for each cluster with the -fasta option (Should be of form cluster_number.fasta).
4. Use a multiple alignment tool (e.g. MUSCLE) to generate a fasta alignment 
file from the cluster centroids fasta file (calc_distmtx.py    is set to read
 fasta format).
5. Run calc_distmtx.py with this file as the infile. Also, specify an outfile.
6. Run rank_and_plot_distmtx.py with outfile from previous step as infile. 
Also, specify an outfile as well as the fasta cluster centroids    file used 
in step 4. Be sure to add a plot title with the -title option as well as the 
number of centroids to rank (if you don't want      all of them ranked) with 
the -rank option. Lastly, set a weight for cluster size (between 0 and 1). A 
higher -size_weight will tell the    algorithm to also pick larger clusters, 
even though they may not constitute the most diverse set of clusters. Values 
that seem to work    well range between 0.1 and 0.4.
7. Run find_reviewed.py with outfile from previous step as infile. Also specify outfile (this is the final output). Lastly, specify the      path to the directory containing all of the cluster fasta files (from step 3).
8. Run find_PDB_clusters_SIFTS.py with reviewed ranked list as infile to find
 PDB structures (if any) for each cluster - make sure to designate an outfile
  as well.
8. (optional) Generate FASTA file from output via Uniprot's ID Mapping Service.