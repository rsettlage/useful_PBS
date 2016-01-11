  #!/bin/bash

###usage qsub /groups/DAC/useful_PBS/PBS_TranscriptAssembly_4_velvet_oases_merge.sh
#PBS -N Sub_S1
#PBS -m a -M rsettlage@vbi.vt.edu
# Set the walltime, which is the maximum time your job can run in HH:MM:SS
#PBS -lwalltime=40:00:00
# Set the number of nodes, and the number of processors per node (up to 12)
# Access group and queue - Do not change
#PBS -W group_list=grl
#PBS -q grl_q
#PBS -lnodes=1:ppn=10
######add modules
export MODULEPATH=/apps/packages/bio/modulefiles:$MODULEPATH
export MODULEPATH=/apps/modulefiles/stats:$MODULEPATH
module load velvet
module load oases
module list
###print PBS script
PBS_script="/groups/DAC/useful_PBS/PBS_TranscriptAssembly_4_velvet_oases_merge.sh"
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

kmer_merge=27
merge_DIR=merge_S1_R3

velveth $merge_DIR $kmer_merge -fasta -long /groups/DAC/Subbiah_Aug2014/working_data/zz_test/S1_R3_cc2*/contigs.fa 

# move to node dir
cd $merge_DIR
cp -v ./* $node_DIR
cd $node_DIR

echo #############
echo velveth portion done, starting velvetg
velvetg ./ -read_trkg yes -conserveLong yes
cp -v -u ./* $PBS_O_WORKDIR/$merge_DIR/
echo #############
echo velvet protion done, starting oases
echo #############
oases ./ -merge yes -min_trans_lgth 200
cp -v -u ./* $PBS_O_WORKDIR/$merge_DIR


###summarize results
perl /groups/DAC/useful_perl/getContigSummary.pl contigs.fa >contigs.fa.stats
perl /groups/DAC/useful_perl/count_loci_hash.pl transcripts.fa
tail -n 3 transcripts.fa.stats

###########cleanup and exit
cp -v -u ./* $PBS_O_WORKDIR/$merge_DIR
rm -r $node_DIR

echo end:
date
echo finished
exit
