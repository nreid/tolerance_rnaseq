#!/bin/bash
#SBATCH -p hi
#SBATCH --cpus-per-task=24
#SBATCH --mem=16000
###### number of nodes
###### number of processors

module load samtools/1.2 
module load boost/1.55.0

#various relevant directories and files
genomebase=/home/nreid/popgen/kfish3/killifish20130322asm
genome=/home/nreid/popgen/kfish3/killifish20130322asm.fa

fc=/home/nreid/bin/subread-1.4.6-p4-source/bin/featureCounts
saf=/home/nreid/popgen/kfish3/kfish2rae5g.main.pub.corrected.exons.saf
inbam=$(cat bowtie.bams.list | tr '\n' ',' | sed 's/,/ /g')

root=tolerance.counts
outdir=/home/nreid/tolerance_rnaseq/featurecounts/$root.fc
outdir2=/home/nreid/tolerance_rnaseq/featurecounts/$root.meta.fc

#$fc -Q 10 -T 24 -a $saf -F SAF -o $outdir $inbam
$fc -Q 10 -T 24 -f -p -a $saf -F SAF -o $outdir2 $inbam
