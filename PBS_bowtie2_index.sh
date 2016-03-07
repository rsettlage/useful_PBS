#!/bin/bash

## qsub /home/rsettlag/useful_PBS/PBS_bowtie2_index.sh -v FILE=/work/newriver/rsettlag/genomes/Bos_taurus.UMD3.1.dna.toplevel.fa

###########################################################################
## environment & variable setup

####### job customization
#PBS -N BT2_idx
## merge stdout and stderr
#PBS -j oe
#PBS -m a -M rsettlag@vt.edu
# ARC needs:
#PBS -l walltime=0:40:00
#PBS -l nodes=1:ppn=1
#PBS -W group_list=newriver
#PBS -q normal_q
#PBS -A arc_test
####### end of job customization
# end of environment & variable setup
###########################################################################
module purge
module load gcc/5.2.0
###print PBS script
PBS_script="/home/rsettlag/useful_PBS/PBS_bowtie2_index.sh"
echo '#############################################################################'
more $PBS_script
echo '#############################################################################'
###########################################################################
echo start:
date

###$TMPFS

database_NAME=${FILE##*/}
database_BASE=${database_NAME%.*}
database_DIR=${FILE%/*}

echo $database_NAME
echo $database_BASE

cd $database_DIR

###########################################################################
## begin execution stage

## actual command and its options to be executed by this job
echo bowtie2-build -f $database_NAME $database_BASE
EXEC_CMD="/home/rsettlag/bowtie2/bowtie2-2.2.7/bowtie2-build -f $database_NAME $database_BASE"
$EXEC_CMD

## end execution stage
###########################################################################

echo end:
date

exit;
