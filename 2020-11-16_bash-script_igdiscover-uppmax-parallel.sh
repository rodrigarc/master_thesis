find . -type f -name "*.yaml" -exec sed -i 's/barcode_length_3prime: 0/barcode_length_3prime: 21/g' {} \;
 
#!/bin/sh -l
  
 #SBATCH -A snic2020-15-190
 #SBATCH -p node
 #SBATCH -n 9
 #SBATCH -t 0-12:00:00
 #SBATCH -J ig-cirelli
 
conda activate igdiscover
 
foo () {
igdiscover init --db path/to/my/database/ --reads1 "${f%.gz}-B1-D14_S2_L001_R1_001.fastq.gz" "/proj/snic2020-16-131/results/2020-11-16/${f%.gz}" 

cd "$f" 
vim igdiscover.yaml


igdiscover init --db path/to/my/database/ --reads1 mylibrary.1.fastq.gz myexperiment


igdiscover clonoquery --cdr3-core 2:-2 --mismatches 0.2 --aa --summary "${f%.gz}_clonoquery_summary.txt" "$f" 2020-08-10_public_clonotypes_summary_for_clonoquery_added_FR4-0.txt
    }
 
for f in *.gz; do foo "$f" & done


  sbatch -A snic2020-15-190 -n 2 -t 00:30:00 -J test01 --wrap "igdiscover run -j 2 iteration-01/filtered.tab"