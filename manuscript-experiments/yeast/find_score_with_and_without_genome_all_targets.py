import os

def main():
    list_grna_scores = []
    for file in os.listdir('.'):
        if file.endswith('.fasta'):
            print(file)

if __name__ == "__main__":
    main()
