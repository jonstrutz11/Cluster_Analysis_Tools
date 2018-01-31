# Note: Not developed by Tyo lab at Northwestern University. Taken from
# http://biopython.org/wiki/Sequence_Cleaner and modified slightly. Used to
# clean up sequences ( remove duplicates) before further processing. Also
# note that this was originally designed for DNA seqs (as evidenced by the
# unknown nucleotide percentage allowed - por_n), but still works for protein
# sequences without the por_n input.

import sys
from Bio import SeqIO


def sequence_cleaner(fasta_infile, fasta_outfile=None, min_length=0.0,
                     por_n=100.0):
    # Create our hash table to add the sequences
    sequences = {}

    # Using the Biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(fasta_infile, "fasta"):
        # Take the current sequence
        sequence = str(seq_record.seq).upper()
        # Check if the current sequence is according to the user parameters
        if (len(sequence) >= min_length and
           (float(sequence.count("N"))/float(len(sequence)))*100 <= por_n):
            # If the sequence passed in the test "is it clean?" and it isn't in
            # the hash table, the sequence and its id are going to be in the
            # hash
            if sequence not in sequences:
                sequences[sequence] = seq_record.id
            # If it is already in the hash table, we're just gonna
            # concatenate the ID
            # of the current sequence to another one that is already in the
            # hash table
            else:
                sequences[sequence] += "_" + seq_record.id

    # Write the clean sequences
    if not fasta_outfile:
        fasta_outfile = "clear_" + fasta_infile
    with open(fasta_outfile, "w+") as output_file:
        # Just read the hash table and write the file in fasta format
        for sequence in sequences:
            output_file.write(">" + sequences[sequence] + "\n" + sequence + "\n")

    print("CLEAN!!!\nPlease check clear_" + fasta_infile)


userParameters = sys.argv[1:]

try:
    if len(userParameters) == 1:
        sequence_cleaner(userParameters[0])
    elif len(userParameters) == 2:
        sequence_cleaner(userParameters[0], userParameters[1])
    elif len(userParameters) == 3:
        sequence_cleaner(userParameters[0], userParameters[1], float(
            userParameters[2]))
    elif len(userParameters) == 4:
        sequence_cleaner(userParameters[0], userParameters[1], float(
            userParameters[3]), float(userParameters[3]))
    else:
        print("Please provide an input fasta file")
except TypeError:
    print("There is a problem!")
