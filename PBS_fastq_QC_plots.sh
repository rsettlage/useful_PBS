#!/bin/bash
#usage qsub /groups/DAC/useful_PBS/PBS_fastq_QC_plots.sh -v instrument=HiSeq,FILE='{}' \;
###
###for HiSeq:
#find ./ -iname '*001.fastq.gz' -exec qsub /groups/DAC/useful_PBS/PBS_fastq_QC_plots.sh -v instrument=HiSeq,FILE='{}' \;
###for HiSeq, if doing data that has been run through the trimming functions, use instrument = HiSeq-trimmed
#find ./ -iname '*001_AT_QT.fastq.gz' -exec qsub /groups/DAC/useful_PBS/PBS_fastq_QC_plots.sh -v instrument=HiSeq-trimmed,FILE='{}' \;
###for HiSeq, if doing data that has been run through the trimming functions and is paired, use instrument = HiSeq-trimmed-paired
#find ./ -iname '*001_AT_QT.paired_matched.fastq.gz' -exec qsub /groups/DAC/useful_PBS/PBS_fastq_QC_plots.sh -v instrument=HiSeq-trimmed-paired,FILE='{}' \;
###
###for MiSeq:
#find ./ -iname '*.fastq.gz' -exec qsub /groups/DAC/useful_PBS/PBS_fastq_QC_plots.sh -v instrument=MiSeq,FILE='{}' \;
###
#launch in the subdirectory with all the files, need to change this to take a file mask
####### job customization
## name our process
##        <----------------------------------------------set the next line
#PBS -N dataQC
## merge stdout and stderr
#PBS -j oe
## e-mail us when the job: '-M b' = begins, '-M a' = aborts, '-M e' = ends
#PBS -m a -M rsettlage@vbi.vt.edu   
#PBS -lwalltime=2:00:00
################## Access group and queue, use one or the other############################
###PBS -W group_list=sfx
###PBS -q sfx_q
#PBS -W group_list=sfx
#PBS -q sandybridge_q
#PBS -lnodes=1:ppn=1
########  OR  ##############
#####'PBS -W group_list=sfxsmp
#####'PBS -q sfxsmp_q
#####'PBS -lnodes=1:ppn=16:smp
####### end of job customization
export MODULEPATH=/apps/packages/bio/modulefiles:$MODULEPATH
export MODULEPATH=/apps/modulefiles/stats:$MODULEPATH
module load bio/R/3.1.0
###print PBS script
PBS_script="/groups/DAC/useful_PBS/PBS_fastq_QC_plots.sh"
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
## begin execution stage # Below here enter the commands to start your job

if [ ${instrument} == 'HiSeq' ];then
	echo this looks like raw HiSeq data, do:
	CMD="R CMD BATCH -fastq.gz -001.fastq.gz -$FILE /groups/DAC/useful_R/illumina_HiSeq_QC_args.R"
	echo $CMD
	${CMD}
fi

if [ ${instrument} == 'HiSeq-trimmed' ];then
	echo this looks like trimmed HiSeq data, do:
	CMD="R CMD BATCH -_AT_QT.fastq.gz -001_AT_QT.fastq.gz -$FILE /groups/DAC/useful_R/illumina_HiSeq_QC_args.R"
	echo $CMD
	${CMD}
fi

if [ ${instrument} == 'HiSeq-trimmed-paired' ];then
	echo this looks like trimmed HiSeq data, do:
	CMD="R CMD BATCH -_AT_QT.paired_matched.fastq.gz -001_AT_QT.paired_matched.fastq.gz -$FILE /groups/DAC/useful_R/illumina_HiSeq_QC_args.R"
	echo $CMD
	${CMD}
fi

if [ ${instrument} == 'MiSeq' ];then
	echo this looks like trimmed MiSeq data, do:
	CMD="R CMD BATCH -$FILE /groups/DAC/useful_R/illumina_MiSeq_QC_args.R"
	echo $CMD
	${CMD}
fi

exit;