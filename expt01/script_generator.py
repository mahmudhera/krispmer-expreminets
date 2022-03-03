start = 6
end = 15

f = open('script.sh', 'w')
for i in range(start, end+1):
    num = str(i)
    str_ = 'python main.py ../runs-on-genomes/staphylococcus-aureus/mer_counts.jf ../runs-on-genomes/staphylococcus-aureus/target' + num + '.fasta ../runs-on-genomes/staphylococcus-aureus/scores' + num + ' 3 staphylococcus-aureus/out_for_target'+num
    f.write(str_ + '\n')
f.close()
