num_targets = todo
genome_index_name = 'genome.fasta.index'

for i in range(1, num_targets+1):
    i_str = str(i)
    print('python kmer_generator_from_target.py target'+i_str+'.fasta')
    print('guidescan enumerate -o target' + i_str + '_guidescan_out -f tmp ' + genome_index_name)
