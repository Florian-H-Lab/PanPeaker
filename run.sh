#!/bin/bash
#$ -N PanPeaker
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=4G
#$ -o /scratch/bi03/heylf/Snakemake_Cluster_Logs/
#$ -e /scratch/bi03/heylf/Snakemake_Cluster_Logs/

source activate panpeaker
python PanPeaker.py \
--piranha "\"-s -p 0.05 -b 50 -i 50 -u 10\"" \
--pureclip "\"-iv 'chr1;'\"" \
-i /scratch/bi03/heylf/PanPeaker/test/RBFOX2_rep1.bam /scratch/bi03/heylf/PanPeaker/test/RBFOX2_rep2.bam \
-b /scratch/bi03/heylf/PanPeaker/test/RBFOX2_rep1.bai /scratch/bi03/heylf/PanPeaker/test/RBFOX2_rep2.bai \
-c /scratch/bi03/heylf/PanPeaker/test/Input.bam \
-k /scratch/bi03/heylf/PanPeaker/test/Input.bai \
-nt 1 \
-g ../genomes/hg19/GRCh37.p13.genome.fa -o /scratch/bi03/heylf/PanPeaker \
--chr_sizes ../genomes/hg19/hg19_chr_sizes.txt \
--seed 123 --refinement
source deactivate panpeaker


