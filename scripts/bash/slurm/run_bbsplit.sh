#!/bin/bash
#SBATCH --array 0-51
#SBATCH -J LB_hr 
#SBATCH --qos=normal
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH -e /zfshome/sukem127/error_hr_%A_%a.log
#SBATCH -o /zfshome/sukem127/output_hr_%A_%a.log



#set -euo pipefail
echo "Started script"
_sample_path=$1
_read_folder=$2
_out_path=$3
#rm -r  slurm-* || true 
#ret=$?
#if [ $ret= 0 ]; then 
#	echo "Sucessfully deleted  old logs"; 
#else
#	echo "Could not delete old logs"; 
#fi

#rm -rf  slurm_out/reads_to_bins/ 
#ret=$?
#if [ $ret = 0 ]; then 
#	echo "Sucessfully deleted old folder"; 
#else
#	echo "Could not deleted  old fodler"; 
#fi
#mkdir -p slurm_out/reads_to_bins/  &&
#ret=$?
#if [ $ret = 0 ]; then 
#	echo "Sucessfully created log folder"; 
#else
#	echo "Could no create log folder"; 
#fi
source ~/.bashrc
mamba activate bbmap

echo "Loading samples from ${_sample_path}"
echo "Loading reads from ${_read_folder}"
echo "Saving to from ${_out_path}"

IFS=$'\n' read -d '' -r -a lines < "${_sample_path}"
_my_sample=${lines[ $SLURM_ARRAY_TASK_ID ]}
echo "Working on ${_my_sample}"

folder=${_read_folder}
line=${_my_sample}
bbsplit.sh -Xmx20g usejni=t unpigz=t threads=2 minid=0.90 ambiguous2=toss  \
path=${_out_path}/bin_index/ \
in="${folder}/${line}/filter/${line}_r1_nohost.fastq.gz" \
in2="${folder}/${line}/filter/${line}_r2_nohost.fastq.gz" \
basename=./splited-read-${line}/%.out.fastq.gz 

echo "finished sample on ${_my_sample}"

exit 0

