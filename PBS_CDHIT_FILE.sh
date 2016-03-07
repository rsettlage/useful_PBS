#! /bin/bash

## all PBS/torque specific commands are prepended with '#PBS'
####qsub /groups/DAC/useful_PBS/PBS_CDHIT_FILE.sh -v Query=/groups/DAC/vector_contaminants_file/adaptors_contaminants.fa
####

###########################################################################
## PBS stuff
####### job customization
#PBS -N hCDH_T_red
## merge stdout and stderr
#PBS -j oe
## e-mail us when the job: '-M b' = begins, '-M a' = aborts, '-M e' = ends
#PBS -m a -M rsettlage@vbi.vt.edu
## request queue sfx_q, sfxsmp with group list as sfx or sfxsmp, need to add a :smp to the nodes section
#PBS -W group_list=grl
#PBS -q grl_q
#PBS -l nodes=1:ppn=32:grlsmp
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


query_FILE_95_EXT=".95"
query_FILE_90_EXT=".90"
query_FILE_80_EXT=".80"
query_FILE_60_EXT=".60"
query_FILE_60_EXT=".50"
query_FILE_95=$query_FILE$query_FILE_95_EXT
query_FILE_90=$query_FILE$query_FILE_90_EXT
query_FILE_80=$query_FILE$query_FILE_80_EXT
query_FILE_60=$query_FILE$query_FILE_60_EXT
query_FILE_50=$query_FILE$query_FILE_50_EXT
echo results will go in file $query_FILE_rep

results="_hCDHIT_red"
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

echo cd-hit -i $query_FILE -o $query_FILE_95 -c 0.95 -n 5 -g 1 -G 0 -aS 0.8  -d 0 -p 1 -T 16 -M 0 

`cd-hit -i $query_FILE -o $query_FILE_95 -c 0.95 -n 5 -p 1 -T 31 -M 0 -B 1 > db_95.log`
`cd-hit -i $query_FILE_95 -o $query_FILE_90 -c 0.90 -n 5 -p 1 -T 31 -M 0 -B 1 > db_90.log`
`cd-hit -i $query_FILE_90 -o $query_FILE_80 -c 0.80 -n 5 -p 1 -T 31 -M 0 -B 1 > db_80.log`
`cd-hit -i $query_FILE_80 -o $query_FILE_60 -c 0.60 -n 4 -p 1 -T 31 -M 0 -B 1 > db_60.log`
`cd-hit -i $query_FILE_60 -o $query_FILE_50 -c 0.50 -n 3 -p 1 -T 31 -M 0 -B 1 > db_50.log`

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