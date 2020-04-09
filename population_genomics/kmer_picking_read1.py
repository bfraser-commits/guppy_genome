import os

## Load the male_kmers dictionary
male_kmers = {}

# Define the path to the kmer matrix file
in_file = "$PATH/female_filtered_kmer_matrix.tsv"

# Obtain a file handle                                                                                                                                                                                         
fh = open(in_file, 'r')

# Iterate over each line in the file                                                                                                                                                                           
for line in fh:
    # Remove newline characters                                                                                                                                                                                       
    line = line.strip('\n')

    # Split the TSV file                                                                                                                                                                                   
    parts = line.split('\t')
    # the kmers are in the first column of the tsv file, so the first array
    kmer = parts[0]
    ## males (females already filtered)
    GL_males = int(parts[7])
    GH_males = int(parts[8])
    UM_males = int(parts[9])
    LM_males = int(parts[10])
    UQ_males = int(parts[11])
    LO_males = int(parts[12])

## Print kmers that are only >= mean and <= 95% confidence limits males (from female_filtered_matrix)
    if (GL_males >= 6 and GL_males <= 18 and GH_males >= 7 and GH_males <= 24 and ML_males >= 5 and ML_males <= 18 and \
    	MH_males >= 6 and MH_males <= 23 and OL_males >=3 and OL_males <=12 and OH_males >=10 and OH_males <=36):
        male_kmers[kmer] = 1
        
fh.close() 

kmer_len = 31

## File path for the population specific fasta file
file_path = "$PATH/read1.fasta"

## Output file name for population specific read1_ymers
out_file = "$PATH/kmer_picking_outputs/read_1_ymers.fasta"

# Iterate over each line in the file                                                                                                                                                                           
fh = open(file_path, 'r')
fh1 = open(out_file, "w")

for id in fh:

# Remove newline character from id line (starts with >)           
    id = id.strip('\n')

# If ID is > then the sequence is the next line (as fasta file), strip newline character from sequence
    if id == ">" + id[1:]:
        seq = fh.next()
        seq = seq.strip('\n')

# Iterate over the reads looking for num_kmers (if they appear in the male kmers dictionary)

        for i in range(0, num_kmers):
            j = i + kmer_len
            kmer = seq[i:j]
            if kmer in male_kmers:
                guanapo_males = id + "\n" + seq + "\n"
                fh1.write(guanapo_males)

## Close the file                                                                                                                                                                                               
fh.close()
fh1.close()
