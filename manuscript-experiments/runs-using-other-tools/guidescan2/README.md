## Steps to implement in script

1. For every target, prepare a kmer-file
`python generate_kmers.py --kmer-length 20 target3.fasta --min-chr-length 1 > tmp`
1. For every target, enumerate the grnas
`guidescan enumerate -m 3 -f tmp test_index --output tmp2`
1. Read the tmp2 and extract only the grnas
`tail -n +3 tmp2 | cut -d$'\t' -f10 > tmp3`

## Steps to perform to run using guidescan
1. Copy all targets in directory "inputs"
1. Copy the genome in directory
1. Make a database using:
`guidescan index <FASTA-GENOME>`
1. Copy the `generate_kmers.py` script in directory
1. Make sure guidescan is installed
1. Run script
