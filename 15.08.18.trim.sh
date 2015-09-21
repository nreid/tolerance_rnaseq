#!/bin/bash -l
#SBATCH -J array_job
#SBATCH -o array_job_out_%A_%a.txt
#SBATCH -e array_job_err_%A_%a.txt
#SBATCH --array=1-728
#SBATCH -p hi
#SBATCH --cpus-per-task=4



trimmo=java\ -jar\ /home/nreid/bin/trimmomatic-0.33/trimmomatic-0.33.jar

#fq root
fq1=$(find ~/tolerance_rnaseq/ -type f -name "*.gz" | \
grep R1 | \
grep -v "ut.fastq" | \
grep -v "pt.fastq" | \
sed -n $(echo $SLURM_ARRAY_TASK_ID)p)

echo $fq1
fq2=$(echo $fq1 | sed 's/_R1_/_R2_/')
echo $fq2

ofq1=$(echo $fq1 | sed 's/fastq/pt\.fastq/')
ofq1u=$(echo $fq1 | sed 's/fastq/ut\.fastq/')
ofq2=$(echo $fq2 | sed 's/fastq/pt\.fastq/')
ofq2u=$(echo $fq2 | sed 's/fastq/ut\.fastq/')
echo $ofq1
echo $ofq1u
echo $ofq2
echo $ofq2u

$trimmo PE -threads 4 -phred33 \
$fq1 $fq2 $ofq1 $ofq1u $ofq2 $ofq2u \
ILLUMINACLIP:NEBnextAdapt.fa:2:30:10 \
LEADING:5 \
TRAILING:5 \
MINLEN:50
