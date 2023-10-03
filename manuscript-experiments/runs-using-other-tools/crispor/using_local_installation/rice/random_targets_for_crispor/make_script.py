f = open('script.sh', 'w')
for i in range(100):
    f.write(f'python /var/www/html/crispor.py pz9Osativa target_{i}.fasta target_{i}_results.tsv\n')
    f.flush()
f.close()
