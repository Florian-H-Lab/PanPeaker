#!/bin/bash
#$ -N PanPeaker
#$ -cwd
#$ -pe smp 10
#$ -l h_vmem=4G
#$ -o /scratch/bi03/heylf/Snakemake_Cluster_Logs/
#$ -e /scratch/bi03/heylf/Snakemake_Cluster_Logs/

python PanPeaker.py --piranha "\"-s -p 0.05 -b 50 -i 50 -u 10\"" --pureclip "\"-iv 'chr1;chr2;chr3;'\"" \
-i /scratch/bi03/heylf/PanPeaker/test/RBFOX2_rep1.bam -b /scratch/bi03/heylf/PanPeaker/test/RBFOX2_rep1.bai \
-c /scratch/bi03/heylf/PanPeaker/test/Input.bam -k /scratch/bi03/heylf/PanPeaker/test/Input.bai -nt 10 \
-g ../genomes/hg19/GRCh37.p13.genome.fa -o /scratch/bi03/heylf/PanPeaker



