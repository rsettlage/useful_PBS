#!/bin/bash

###usage qsub /groups/DAC/useful_PBS/PBS_TranscriptAssembly_NoConvey_2_velveth_velvetg_oases.sh -v KMER=33
###for i in 41 39 37 35 33 31 29 27 25 23; do qsub /groups/DAC/useful_PBS/PBS_TranscriptAssembly_NoConvey_2_velveth_velvetg_oases.sh -v KMER=$i; done
###for i in 81 79 77 75 73 71 69; do qsub /groups/DAC/useful_PBS/PBS_TranscriptAssembly_NoConvey_2_velveth_velvetg_oases.sh -v KMER=$i; done
###pass in kmer
#PBS -N Mattison_GRL4553

# Set the walltime, which is the maximum time your job can run in HH:MM:SS
#PBS -lwalltime=50:00:00
#PBS -lnodes=1:ppn=12
#PBS -W group_list=grl
#PBS -q grl_q
######PBS -l mem=256GB
####only specify convey-ex01 for large genome/transcriptomes needing more memory
###print PBS script
PBS_script="/groups/DAC/useful_PBS/PBS_TranscriptAssembly_NoConvey_2_velveth_velvetg_oases.sh"
######add modules
export MODULEPATH=/apps/packages/bio/modulefiles:$MODULEPATH
export MODULEPATH=/apps/modulefiles/stats:$MODULEPATH
module add velvet
module add oases
module list
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
results_DIR=$work_DIR/$KMER
mkdir $results_DIR ####<------good thing to do
node_DIR=/localscratch/$JOBID
mkdir $node_DIR 
df -h /localscratch/
echo originating directory is $work_DIR
echo node directory is $node_DIR
echo results directory is $results_DIR
echo
cd $work_DIR
pwd
##########################################################################
# Below here enter the commands to start your job

#print the starting time of the job
echo start:
date

pwd
echo kmer is $KMER
cd $results_DIR
echo now in results directory...
pwd
ln -s $work_DIR/CnyUnifiedSeq
##ln -s $work_DIR/CnyUnifiedSeq or Sequences <-------------------------------------choose for previous line ------######
ls -la --block-size=MB
cp -v ./* $node_DIR
cd $node_DIR
echo now in node directory...
pwd
ls -lah

echo ########
echo starting velveth, takes about 10 min
velveth ./ $KMER -reuse_binary
##velveth ./ $KMER -reuse_binary <---------------------change if able to use binary otherwise reuse_Sequences -----##########
cp -u -v -r $node_DIR/* $results_DIR
echo ########
echo starting velvetg, takes about 3 hours
velvetg ./ -read_trkg yes -unused_reads no -exp_cov auto -ins_length 188
cp -u -v -r $node_DIR/* $results_DIR
echo starting oases, takes about 2 hours
echo ######## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<---------------------------------insert size is fragment size with reads...NOT adaptors
###oases ./ -min_trans_lgth 200 -ins_length 160
oases ./ -min_trans_lgth 200 -ins_length 240

###summarize results
echo finished, do summaries
perl /groups/DAC/useful_perl/count_loci_hash.pl transcripts.fa
tail -n 3 transcripts.fa.stats


###clean up and exit
ls -lah
cp -u -v -r $node_DIR/* $results_DIR
rm $node_DIR -r
cd $results_DIR
pwd
ls -lah

#print the end time of the job
echo end:
date
echo $PBS_JOBID.out

df -h
hostname
echo $node_DIR
echo $KMER
echo finished!
exit
