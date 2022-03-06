start = 6
end = 15

f = open('script.sh', 'w')
for i in range(start, end+1):
    num = str(i)
    str_ = 'python main.py ../runs-on-genomes/rhodobacter-sphaeroides/mer_counts.jf ../runs-on-genomes/rhodobacter-sphaeroides/target' + num + '.fasta ../runs-on-genomes/rhodobacter-sphaeroides/scores' + num + ' 3 rhodobacter-sphaeroides/out_for_target'+num
    f.write(str_ + '\n')
f.close()
