#!/bin/bash

## qsub /groups/DAC/useful_PBS/PBS_bowtie2_index.sh -v FILE=/groups/DAC/Mattison_March2014/working_data/Pecan_Trinity_NSS/Pecan_Trinity_NSS.fasta

###########################################################################
## environment & variable setup

####### job customization
## name our process
#PBS -N BT2_idx
## merge stdout and stderr
#PBS -j oe
## e-mail us when the job: '-M b' = begins, '-M a' = aborts, '-M e' = ends
#PBS -m a -M rsettlage@vbi.vt.edu
# Access group and queue
#PBS -W group_list=grl
#PBS -q grl_q
# Set the walltime, which is the maximum time your job can run in HH:MM:SS
#PBS -lwalltime=5:00:00
# Set the number of nodes, and the number of processors per node (up to 12)
#PBS -lnodes=1:ppn=20:grlhighmem
#######PBS -lnodes=1:ppn=16:smp ##not used for the sfx queue
####### end of job customization
# end of environment & variable setup
###########################################################################
export MODULEPATH=/apps/packages/bio/modulefiles:$MODULEPATH
export MODULEPATH=/apps/modulefiles/stats:$MODULEPATH
module load bio/bowtie2/2.1.0
###print PBS script
PBS_script="/groups/DAC/useful_PBS/PBS_bowtie2_index.sh"
echo '#############################################################################'
more $PBS_script
echo '#############################################################################'
###########################################################################
echo start:
date

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
EXEC_CMD="bowtie2-build -f $database_NAME $database_BASE"
$EXEC_CMD

## end execution stage
###########################################################################

echo end:
date

exit;
