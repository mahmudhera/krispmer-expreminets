from Bio import SeqIO
import random

fasta_filename = 'one_assembly.fasta'
num_targets_each_chr = 15
target_len_low = 2000
target_len_high = 3000

target_id = 1
for seq_record in SeqIO.parse(fasta_filename, "fasta"):
	num_bases = len(seq_record)
	start_positions = [random.randint(1,num_bases-target_len_high-1) for i in range(num_targets_each_chr)]
	end_positions = [start_positions[i] + random.randint(target_len_low, target_len_high) for i in range(num_targets_each_chr)]
	for pos in list(zip(start_positions, end_positions)):
		#print(str(seq_record.id)+":"+str(pos[0])+"-"+str(pos[1]))
		f = open('target' + str(target_id), 'w')
		f.write ('> target' + str(target_id) + '\n')
		target_id += 1
		f.write (str(seq_record.seq[pos[0]:pos[1]]) + '\n')
		f.close()
