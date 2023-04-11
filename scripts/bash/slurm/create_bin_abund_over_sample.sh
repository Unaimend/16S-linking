#!/bin/bash
#SBATCH -J LB_hr 
#SBATCH --qos=normal
#SBATCH --time=4:30:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH -e /zfshome/sukem127/error_hr_%A_%a.log
#SBATCH -o /zfshome/sukem127/output_hr_%A_%a.log

args=("$@")  # Access all arguments as an array
echo "Current Task ID is ${SLURM_ARRAY_TASK_ID}"
echo "Currently working on ${args[$SLURM_ARRAY_TASK_ID]}"
#this is easier todo with $(basename <file> suffix)
sample3=$(basename ${args[$SLURM_ARRAY_TASK_ID]})
sample2=$(basename ${sample3})
sample="${sample2%.*}"
echo "Starting with ${sample}"
echo "Folder is ${args[0]}"

bbmap.sh -Xmx20g unpigz=t threads=2 minid=0.90 \
statsfile=bin-statsfiles/${sample}.statsfile \
scafstats=bin-scafstats-statsfiles/${sample}.scafstats \
covstats=bin-cov-statsfiles/${sample}.covstat \
rpkm=bin-rpkm-statsfiles/${sample}.rpkm \
sortscafs=f nzo=f ambiguous=all \
in="${args[0]}/${sample}/filter/${sample}_r1_nohost.fastq.gz" \
in2="${args[0]}/${sample}/filter/${sample}_r2_nohost.fastq.gz"



exit 0
