#!/bin/bash -l
 
for f in *.gz; 
do igdiscover clonoquery --cdr3-core 2:-2 --mismatches 0.2 --aa "$f" 2020-08-10_public_clonotypes_summary_for_clonoquery_added_FR4-0.txt > "${f%.gz}_clonoquery_full.txt" 
done

