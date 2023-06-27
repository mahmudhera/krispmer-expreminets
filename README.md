# krispmer-expreminets
Codes of all experiments to prepare results are in this repository, along with comments. This parent readme will also store directions and guidelines as well to recreate the experiments. Will also contain meeting notes.


## Log in info
experiments are run in this machine:
131.215.78.41
ajker jonno

## Installation tips
Had to do the following:
```
export PYTHONPATH=/home/atif/krispmer/soft/krispmer-github-repo/kRISP-mER/kRISP-meR_source:$PYTHONPATH
export PATH=/home/atif/krispmer/soft/bowtie2-2.4.5-linux-x86_64:$PATH
export PATH=/home/atif/guidescan2/guidescan-cli-master/build/bin:$PATH
```

## Todo for the tool
1. Add on-tgt activity scores
1. Fix package error


Time to run on human genome:
----------------------------
Hamming distance = 3, time = ~14 sec, mem = 11.7 GB
Hamming distance = 2, time = ~7  sec, mem = 0.64 GB
Hamming distance = 1, time = ~6  sec, mem = 0.14 GB
