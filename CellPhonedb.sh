#!/bin/bash
#SBATCH -t 4-18:00:00
#SBATCH --mem=40000M
#SBATCH -J Cellphonedb
#SBATCH -p himem
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -o %x-%j.out


source ~/.bashrc

source activate scRNAseq

cd /cluster/projects/macparland/TA/LiverMap2.0/RawData/Analysis/EmptyDrops/Subcluster/ManualAnnotation

SAMPLE=$1

META="${SAMPLE}_meta.txt"
COUNTS="${SAMPLE}_counts.txt"
OUT="${SAMPLE}_cellphonedb"

mkdir -p $OUT
cellphonedb method statistical_analysis $META $COUNTS --output-path=$OUT --threads=1 --debug-seed 42 --database ~/.cpdb/releases/v2.0.0/cellphone.db --counts-data hgnc_symbol --iterations 5


