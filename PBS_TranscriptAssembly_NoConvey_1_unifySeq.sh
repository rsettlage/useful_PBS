#!/bin/bash

###usage qsub /groups/DAC/useful_PBS/PBS_TranscriptAssembly_NoConvey_1_unifySeq.sh
#PBS -N Grabau_CombinedWlines
#PBS -j oe
# Set the walltime, which is the maximum time your job can run in HH:MM:SS
#PBS -lwalltime=040:00:00
#PBS -lnodes=1:ppn=4
# Access group and queue - Do not change
#PBS -W group_list=sfx
#PBS -q sfx_q
######add modules
export MODULEPATH=/apps/packages/bio/modulefiles:$MODULEPATH
export MODULEPATH=/apps/modulefiles/stats:$MODULEPATH
module add velvet
module add oases
module list
###print PBS script
PBS_script="/groups/DAC/useful_PBS/PBS_TranscriptAssembly_NoConvey_1_unifySeq.sh"
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
echo originating directory is $work_DIR
echo
cd $work_DIR
pwd
##########################################################################
# Below here enter the commands to start your job

#short="/groups/DAC/Grabau_Oct2013/BACblast/N70/*gz"
#/groups/DAC/Oscillated_turkey/Illx/assembly/L1_R1_AT_QT_run2.paired_matched.fastq.gz /groups/DAC/Oscillated_turkey/Illx/assembly/L1_R2_AT_QT_run2.paired_matched.fastq.gz _
#/groups/DAC/Oscillated_turkey/Illx/assembly/L2_R1_AT_QT_run2.paired_matched.fastq.gz /groups/DAC/Oscillated_turkey/Illx/assembly/L2_R2_AT_QT_run2.paired_matched.fastq.gz
#"

shortPaired="/groups/DAC/Mattison_March2014/working_data/GRL4553_9-20-13-BR2_CTTGTA_L008_R1.fastq /groups/DAC/Mattison_March2014/working_data/GRL4553_9-20-13-BR2_CTTGTA_L008_R2.fastq"
#long_singles="/groups/DAC/Jorge_June2013/2013_04_26_SergioMarshall_PS-LF89-MP-8_reads.fasta.fastq"

###no reason to specify the kmer, but the program still needs it
kmer=21
echo $file_MASK
assembly_name=GRL4553

velveth ./ 71 -fmtAuto -separate -shortPaired1 /groups/DAC/Grabau_Oct2013/BACblast/Wlines/combined_w/W73_R1.fastq.gz /groups/DAC/Grabau_Oct2013/BACblast/Wlines/combined_w/W73_R2.fastq.gz -shortPaired2 /groups/DAC/Grabau_Oct2013/BACblast/Wlines/combined_w/W171_R1.fastq.gz /groups/DAC/Grabau_Oct2013/BACblast/Wlines/combined_w/W171_R2.fastq.gz -create_binary -noHash
###velveth ./$assembly_name/ $kmer -fmtAuto -separate -shortPaired $shortPaired -create_binary -noHash

###########cleanup and exit

echo finished!
exit
