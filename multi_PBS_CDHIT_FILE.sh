#! /bin/bash

## all PBS/torque specific commands are prepended with '#PBS'
####qsub /groups/DAC/useful_PBS/multi_PBS_CDHIT_FILE.sh -v Identity=0.95,Query=/groups/DAC/vector_contaminants_file/adaptors_contaminants.fa
####

###########################################################################
## PBS stuff
####### job customization
#PBS -N CDH_Trinity
## merge stdout and stderr
#PBS -j oe
## e-mail us when the job: '-M b' = begins, '-M a' = aborts, '-M e' = ends
#PBS -m a -M rsettlage@vbi.vt.edu
## request queue sfx_q, sfxsmp with group list as sfx or sfxsmp, need to add a :smp to the nodes section
#PBS -W group_list=sfx
#PBS -q sandybridge_q
#PBS -lnodes=1:ppn=12
####PBS -l naccesspolicy=singlejob or #PBS -W x="NACCESPOLICY:SINGLEJOB"
######## specify resource allocation
## ask for 1 hour of wall time
#PBS -l walltime=120:00:00
######## additional sample resource reservations
## ask for 500MB of memory
###PBS -l mem=500mb
# end of PBS stuff
###########################################################################
## environment & variable setup
export MODULEPATH=/apps/packages/bio/modulefiles:$MODULEPATH
export PATH=/apps/packages/bio/cd-hit-auxtools/current/:$PATH
module load bio/cd-hit/4.6.1
module list
## end of loading modules
## end of environment & variable setup
###########################################################################
echo start:
date
echo jobid is $PBS_JOBID
JOBID=${PBS_JOBID%%.*}
echo jobnumber is $JOBID
###########################################################################
echo
SCRIPT=/groups/DAC/useful_PBS/multi_PBS_CDHIT_FILE.sh
more $SCRIPT
echo
echo
module list
echo
###########################################################################
## begin execution stage

work_DIR=$PBS_O_WORKDIR
cd $work_DIR
local_DIR=/localscratch/$JOBID
mkdir $local_DIR

echo originating directory is $work_DIR
echo node directory is $local_DIR
echo
echo


query_FILE=${Query##*/}
query_DIR=${Query%/*}
echo query is $query_FILE

query_FILE_rep_EXT="."$Identity
query_FILE_rep=$query_FILE$query_FILE_rep_EXT
echo results will go in file $query_FILE_rep

results=_CDHIT_$Identity
result_DIR=$query_DIR/$query_FILE$results
echo result directory is $result_DIR
mkdir $result_DIR
cd $result_DIR
pwd
ln -s $query_DIR/$query_FILE

echo

###cd $work_DIR

cp -v ./* $local_DIR/
cd $local_DIR
echo directory we are working in is...
pwd
echo


ls -lah

echo
echo
pwd
echo directory listing:
ls -la
echo ready to run cd-hit
### tgicl -F $query_FILE -c 6

echo cd-hit-est -i $query_FILE -o $query_FILE_rep -c $Identity -d 0 -n 10 -p 1 -M 0 -T 12

`cd-hit-est -i $query_FILE -o $query_FILE_rep -c $Identity -d 0 -n 10 -p 1 -M 0 -T 12 >$query_FILE.out_log`

echo 
echo files on node:
ls -lah
echo copying result over using cp -v $result_FILE $query_DIR/$result_FILE

cp -v -u ./* $result_DIR

## end execution stage
###########################################################################
echo end:
date

rm -r $local_DIR/

echo $query_FILE
echo finished
exit