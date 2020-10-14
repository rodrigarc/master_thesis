#!/bin/bash -l

for f in *.fasta;
do FastTree -topm 1.0 -close 0.75 -refresh 0.8 -2nd -log ../trees/"${f%.fasta}_logfile.txt" -nt -sprlength 10 -nni 10 -spr 2 -cat 20 -gtr "$f" > ../trees/"${f%.fasta}_tree_file.tre"
done

