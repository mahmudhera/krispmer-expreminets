from Bio import SeqIO
import random

fasta_filename = 'NC_010079.fasta'
num_targets_each_chr = 15
target_len_low = 3000
target_len_high = 4000
reads_filename = 'reads.fastq'
hd = 3

bash_file = open('script.sh', 'w')
target_id = 1
for seq_record in SeqIO.parse(fasta_filename, "fasta"):
	num_bases = len(seq_record)
	start_positions = [random.randint(1,num_bases-target_len_high-1) for i in range(num_targets_each_chr)]
	end_positions = [start_positions[i] + random.randint(target_len_low, target_len_high) for i in range(num_targets_each_chr)]
	for pos in list(zip(start_positions, end_positions)):
		#print(str(seq_record.id)+":"+str(pos[0])+"-"+str(pos[1]))
		f = open('target' + str(target_id) + '.fasta', 'w')
		f.write ('> target' + str(target_id) + '\n')
		f.write (str(seq_record.seq[pos[0]:pos[1]]) + '\n')
		f.close()
		bash_file.write('krispmer ' + reads_filename + ' ' + 'target' + str(target_id) + '.fasta' + ' ' + 'scores' + str(target_id) + ' ' + str(hd) + ' -n\n')
		bash_file.write('echo Done for ' + str(target_id) + '\n')
		target_id += 1
bash_file.close()
