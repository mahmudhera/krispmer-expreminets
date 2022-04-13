Steps to run:

1. Create directory
1. Copy your genome there
1. `guidescan index genome_name`
1. Copy target files in directory
1. Copy kmer_generator script there
1. Copy all_target_handler.py there
1. Edit the file with number of targets and genome index name
1. `python all_target_handler.py > handle_all_targets.sh`
1. Finally, `bash handle_all_targets.sh`
