import pandas as pd

num_targets = 50
target_length = 150
random_seed = 1
fasta_filename = 'genome.fasta'

if __name__ == "__main__":
    df = pd.read_csv('yeast_gRNAs_multiple.sam', delimiter='\t', header=None)
    print(df)
