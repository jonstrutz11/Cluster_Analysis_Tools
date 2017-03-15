# Get first 10 characters of id, and then get size of cluster from FASTA file.
# Note: This assumes that the size is appended to the ID by usearch using the
# -sizeout option (e.g. with cluster_fast).


def get_size(filename):
    size = ''
    size_dict = {}
    # File must be FASTA format
    with open(filename) as infile:
        for line in infile:
            # Get line with id information
            if line[0] == '>':
                # Determine end of id (by reading until '|')
                end = 8
                char = line[end]
                while char != '|':
                    end += 1
                    char = line[end]
                seq_id = line[1:end] + '|'
                # Scan along line until we find string 'size=' and store
                # its value in a dictionary
                for char_num in range(len(line)):
                    if line[char_num:char_num + 5] == 'size=':
                        char_num += 5
                        while line[char_num] != ';':
                            size += line[char_num]
                            char_num += 1
                        size_int = int(size)
                        size_dict[seq_id] = size_int
                        size = ''
    return size_dict
