#!/bin/bash

###usage qsub /groups/DAC/useful_PBS/PBS_TranscriptAssembly_3_velvetg_oases.sh -v FILE=/groups/DAC/Rami_RNAseq/denovo_assemblies/TT_wholeBird1/21/Graph2
###use find with absolute paths
#PBS -N TTL_fish
#PBS -j oe
# Set the walltime, which is the maximum time your job can run in HH:MM:SS
#PBS -lwalltime=30:00:00
# Set the number of nodes, and the number of processors per node (up to 12)
# Access group and queue - Do not change
#PBS -W group_list=grl
#PBS -q grl_q
#PBS -lnodes=1:ppn=12
######add modules
export MODULEPATH=/apps/packages/bio/modulefiles:$MODULEPATH
export MODULEPATH=/apps/modulefiles/stats:$MODULEPATH
module add velvet
module add oases
module list
###print PBS script
PBS_script="/groups/DAC/useful_PBS/PBS_TranscriptAssembly_3_velvetg_oases.sh"
echo '#############################################################################'
more $PBS_script
echo '#############################################################################'
###########################################################################
echo start:
date
echo jobid is $PBS_JOBID
JOBID=${PBS_JOBID%%.*}
echo jobnumber is $JOBID
###########################################################################
################set up the directories on node and for results#########
work_DIR=$PBS_O_WORKDIR
results_DIR=$work_DIR/$JOBID
#####mkdir $results_DIR ####<------good thing to do
node_DIR=/localscratch/$JOBID
mkdir $node_DIR
echo originating directory is $work_DIR
echo node directory is $node_DIR
echo results directory is $results_DIR
echo
cd $work_DIR
pwd
##########################################################################
# Below here enter the commands to start your job

###pass it the location of the Graph2 file, it will cd into that directory and start velvetg and oases

echo moving to directory containing Graph2 file
graph2_DIR=${FILE%/*}
echo $graph2_DIR
cd $graph2_DIR
pwd
echo copy over files
cp -v ./* $node_DIR/

cd $node_DIR/
pwd
velvetg ./ -read_trkg yes -unused_reads yes 
echo ########
echo starting oases
echo ########
oases ./ -min_trans_lgth 200 -insert_length 350

###summarize results
perl /groups/DAC/useful_perl/count_loci_hash.pl transcripts.fa

###########cleanup and exit
cp -u -v ./* $graph2_DIR/
rm -r $node_DIR
cd $graph2_DIR
pwd
ls -lah

echo $FILE
echo finished!
exit