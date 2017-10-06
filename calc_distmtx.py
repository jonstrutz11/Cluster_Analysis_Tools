# Uses biopython to generate a distance matrix text file from a fasta alignment
# file (e.g. calculated using the MUSCLE web service). Format is 'python
# .\calc_distmtx.py infile outfile'.

import argparse
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO

# Get arguments from command line
parser = argparse.ArgumentParser()
parser.add_argument('infile', nargs=1, type=str,
                    help='Alignment filename (.aln) used to generate distance '
                         'matrix.')
parser.add_argument('outfile', nargs=1, type=str,
                    help='Distance matrix output text filename.')
args = parser.parse_args()

# Read alignment file
alignment = AlignIO.read(vars(args)['infile'][0], 'fasta')
# Set to calculate distances based on sequence identity
calculator = DistanceCalculator('identity')
# Generate distance matrix
dm = calculator.get_distance(alignment)

with open(vars(args)['outfile'][0], 'w') as outfile:
    outfile.write(str(dm))
