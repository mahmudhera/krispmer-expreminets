start = 1
end = 105
genome_name = 'staphylococcus-aureus'
genome_name = 'rhodobacter-sphaeroides'
jf_reads_filename = 'krispmer_temp/jf_binary_file.jf'
histo_filename = 'krispmer_temp/k_spectrum_histo_data'

f = open('script.sh', 'w')
for i in range(start, end+1):
    num = str(i)
    str_ = 'python main.py ../runs-on-genomes/' + genome_name + '/mer_counts.jf ../runs-on-genomes/' + genome_name + '/target' + num + '.fasta ../runs-on-genomes/' + genome_name + '/scores' + num + ' 3 ' + genome_name + '/out_for_target'+num + ' -J ./' + jf_reads_filename + ' -H ./' + histo_filename
    f.write(str_ + '\n')
    str_ = "echo done for " + num
    f.write(str_ + '\n')
f.close()
