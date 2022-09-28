import os

def main():
    list_grna_scores = []
    for file in os.listdir('.'):
        if file.endswith('.fasta') and not file.startswith('scores'):
            print(file)

if __name__ == "__main__":
    main()
