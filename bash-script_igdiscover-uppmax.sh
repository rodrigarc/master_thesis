 #!/bin/sh -l
  
 #SBATCH -A snic-2020-15-190
 #SBATCH -p core
 #SBATCH -n 1
 #SBATCH -t 00:15:00
 #SBATCH -J igdiscovertest
 
 #module load conda
 #conda activate igdiscover
 
igdiscover run -j 1
