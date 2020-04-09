import os

## Make the read 1 dictionary
read_1 = {}

# Define the path to the list of reads for read 1 (as a tsv file)

in_file = "$PATH/kmer_picking_outputs/read_1_ymers.tsv"

# Obtain a file handle                                                                                                                                                                                     
fh = open(in_file, 'r')

# Iterate over each line in the file                                                                                                                                                
for line in fh:
    # Remove newline characters                                                                                                                                                                            
           
    line = line.strip('\n')

    parts = line.split('\t')

    reads = parts[0]
    read_1[reads] = 1

fh.close() 

in_file = "$PATH/clean_reads/sample_read2.fasta"

out_file = "$PATH/kmer_picking_outputs/read_2_ymers.fasta"

# Iterate over each line in the file                                                                                                                                                                       
    
fh = open(in_file, 'r')
fh1 = open(out_file, "w")

for id in fh:

# # Remove newline character from id line (starts with >)           
     id = id.strip('\n')

# # If Id is > then the sequence is the next line (as fasta file), strip newline character from sequence
     if id == ">" + id[1:]:
         seq = fh.next()
         seq = seq.strip('\n')

         if id in read_1:
            read_2 = id + "\n" + seq + "\n"
            fh1.write(read_2)
                                                                                                                                                 
fh.close()
fh1.close()
