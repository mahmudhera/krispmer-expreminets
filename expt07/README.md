## Overview
In this experiment, we will determine sequences from the Human genome (GrCh38 or hg38) which have potential gRNA target sites that appear only once in the reference, but multiple times in the CHM1 reference. The idea is that we will run krispmer using the sequencing reads of CHM1, and a reference based tool using hg38. krispmer will hopefully identify differences between multiple target where this gRNA occurs, and exclude this to avoid off-target effects, whereas, the reference based tools may choose this as a high specificity gRNA. Rank of this gRNA would be particularly interesting to investigate.

### Step 1: determining the gRNA target sites which appear once in hg38 and multiple times in CHM1
1. Counted all 23-mers ending in `GG` (or beginning with `CC`) in the hg38 and CHM1 assembly (use `jellyfish count -m 23` to generate the binary file containing k-mer counts)
1. Listed them in readable files (use `jellyfish dump`)
1. Generated all 23-mers that have count=1 in hg38, >1 in CHM1 (TODO: upload these files)

The number of such 23-mers is very large. Therefore, we next only extract those whichmap to the exons of hg38

### Step 2: determining all exons in hg38
1. Collected the hg38 gene annotations from here: http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/refGene.txt.gz
1. Used the following code to list all exon coordinates in bed file
`zcat refGene.txt.gz|sort -u -k13,13|cut -f3,10,11|awk 'BEGIN{OFS="\t"}{split($2,start,",");split($3,end,","); for(i=1;i<length(start);++i){print $1,start[i],end[i]}}' > exons.bed
`

The bed file is available in this repository.