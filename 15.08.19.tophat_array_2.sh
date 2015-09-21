#!/bin/bash -l
#SBATCH -J array_job
#SBATCH -o array_job_out_%A_%a.txt
#SBATCH -e array_job_err_%A_%a.txt
#SBATCH --array=4-75
#SBATCH -p hi
#SBATCH --cpus-per-task=8
###### number of nodes
###### number of processors
####SBATCH -n 1

module load samtools/1.2 
module load boost/1.55.0
module load bowtie2/2.2.5 
module load tophat/2.1.0

#tell bowtie where the indexed genome is
export BOWTIE2_INDEXES=/home/nreid/popgen/kfish3

#various relevant directories and files
genomebase=/home/nreid/popgen/kfish3/killifish20130322asm
genome=/home/nreid/popgen/kfish3/killifish20130322asm.fa
gff=/home/nreid/popgen/kfish3/kfish2rae5g.main.pub.gff.gz

root=$(find /home/nreid/tolerance_rnaseq/AWBIG9001/ -type f -name "*.fastq.gz"  | grep -o "Sample[^/]*" | sort | uniq | sed -n $(echo $SLURM_ARRAY_TASK_ID)p)
outdir=/home/nreid/tolerance_rnaseq/bowtie/$root.out

#first read
fq1=$(find /home/nreid/tolerance_rnaseq/ -type f -name "*.fastq.gz" | grep $root | grep "pt.fastq" | grep "R1_" | sort | tr "\n" "," | sed 's/,$//')
#second read
fq2=$(find /home/nreid/tolerance_rnaseq/ -type f -name "*.fastq.gz" | grep $root | grep "pt.fastq" | grep "R2_" | sort | tr "\n" ",")
#unpaired
fq3=$(find /home/nreid/tolerance_rnaseq/ -type f -name "*.fastq.gz" | grep $root | grep "ut.fastq" | tr "\n" "," | sed 's/,$//')

fq2=$fq2$fq3

echo $fq1
echo $fq2
echo $fq3
echo $outdir
echo $genomebase

tophat -p 8 -o $outdir $genomebase $fq1 $fq2
