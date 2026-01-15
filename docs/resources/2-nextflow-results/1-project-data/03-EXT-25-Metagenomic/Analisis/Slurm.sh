#!/bin/sh
#SBATCH --mail-type=ALL              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=lauteca@ipb.csic.es # Where to send mail 
#SBATCH --nodes=1                    # Use one node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=20            # Number of CPU cores per task
#SBATCH --output=Analisis_%j.out     # Standard output and error log
#SBATCH --get-user-env 
. /etc/profile.d/modules.sh 

/home/eandres/bin/miARma1.7/miARma /mnt/beegfs/Bioinformatica/EXT/18-EXT-24-Laura_Lozano/Analisis/Quality_human.ini
