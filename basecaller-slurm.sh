#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --mail-user=jomojaco@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --exclude=phoenix-09
#SBATCH --mem=100gb
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-task=25
#SBATCH --output=logs/%x.%j.log
#SBATCH --time=06:00:00

set -o pipefail
# cause a bash script to exit immediately when a command fails
set -e
# cause the bash shell to treat unset variables as an error and exit immediately
set -u
# echo each line of the script to stdout so we can see what is happening
# to turn off echo do 'set +o xtrace'
set -o xtrace

POD5DIR=$1
OUTDIR=$2
REF=$3
#/private/home/jomojaco/dorado-0.7.3-linux-x64/bin/dorado basecaller hac $POD5DIR  --device cuda:all --output-dir $OUTDIR
#/private/home/jomojaco/dorado-0.7.3-linux-x64/bin/dorado basecaller demux $OUTDIR -o $OUTDIR

/private/home/jomojaco/dorado-0.7.3-linux-x64/bin/dorado basecaller hac $POD5DIR  --device cuda:all --reference $REF --kit-name SQK-NBD114-24 > $OUTDIR"/aligned_basecalled.bam" 
/private/home/jomojaco/dorado-0.7.3-linux-x64/bin/dorado demux $OUTDIR -o $OUTDIR"/demux" --kit-name SQK-NBD114-24 
 