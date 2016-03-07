#!/bin/bash

###usage qsub /groups/DAC/useful_PBS/PBS_TranscriptAssembly_Trinity.sh
#PBS -N Trinity_M1Q
#PBS -j oe
# Set the walltime, which is the maximum time your job can run in HH:MM:SS
#PBS -lwalltime=050:00:00
# Access group and queue - Do not change
#PBS -W group_list=sfxsmp
#PBS -q sfxsmp_q
#PBS -lnodes=1:ppn=40:smp
######add modules
export MODULEPATH=/apps/packages/bio/modulefiles:$MODULEPATH
export MODULEPATH=/apps/modulefiles/stats:$MODULEPATH
module add bio/trinity/2014-07-17
module rm devel/java/1.8.0_31
module load devel/java/1.7.0_51
module list
###print PBS script
PBS_script="/groups/DAC/useful_PBS/PBS_TranscriptAssembly_Trinity.sh"
echo '#############################################################################'
more $PBS_script
echo '#############################################################################'
###########################################################################
echo start:
date
echo jobid is $PBS_JOBID
JOBID=${PBS_JOBID%%.*}
echo jobnumber is $JOBID
hostname
###########################################################################
################set up the directories on node and for results#########
work_DIR=$PBS_O_WORKDIR
echo originating directory is $work_DIR
echo
cd $work_DIR
pwd
##########################################################################
# Below here enter the commands to start your job

#left="/groups/DAC/Liu_Jan2014/working_data/M2F_R1.fastq /groups/DAC/Liu_Jan2014/working_data/M2H_R1.fastq /groups/DAC/Liu_Jan2014/working_data/M2P_R1.fastq"
#left="/groups/DAC/Dervisis_Apr2014/working_data/Reads_R1.fastq"
#right="/groups/DAC/Liu_Jan2014/working_data/M2F_R2.fastq /groups/DAC/Liu_Jan2014/working_data/M2H_R2.fastq /groups/DAC/Liu_Jan2014/working_data/M2P_R2.fastq"
#right="/groups/DAC/Dervisis_Apr2014/working_data/Reads_R2.fastq"

left="/groups/DAC/Liu_Jan2014/working_data/MIQ_R1.fastq.gz"
right="/groups/DAC/Liu_Jan2014/working_data/MIQ_R2.fastq.gz"
#right="/groups/DAC/Liu_Jan2014/working_data/MIQ_R2.fastq /groups/DAC/Liu_Jan2014/working_data/MID_R2.fastq /groups/DAC/Liu_Jan2014/working_data/MIN_R2.fastq"

mkdir /localscratch/$JOBID
cd /localscratch/$JOBID

cp -v $left ./
cp -v $right ./

pwd
ls -lah

assembly_name=zz_barnacle_MIQ

/apps/packages/bio/trinity/2014-07-17/Trinity --seqType fq --JM 400G --left $left --right $right --CPU 38 --SS_lib_type RF --min_kmer_cov 2 --output $assembly_name

cp -v -R /localscratch/$JOBID/$assembly_name/ $PBS_O_WORKDIR/

rm -rv /localscratch/$JOBID/

###########cleanup and exit

echo finished!
exit
