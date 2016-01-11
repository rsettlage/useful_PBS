#!/bin/bash

###usage find `pwd` -iname '*accepted_hits.bam' -exec qsub /groups/DAC/useful_PBS/PBS_cufflinks.sh -v EXT=,BAM='{}',STRANDED=second,GENOME=,GTF= \;
####### job customization
## name our process
#PBS -N EW_CLO
## merge stdout and stderr
#PBS -j oe
## e-mail us when the job: '-M b' = begins, '-M a' = aborts, '-M e' = ends
#PBS -m a -M rsettlage@vbi.vt.edu
# Access group and queue
#PBS -W group_list=sfx
#PBS -q sandybridge_q
# Set the walltime, which is the maximum time your job can run in HH:MM:SS
#PBS -lwalltime=200:00:00
# Set the number of nodes, and the number of processors per node (up to 12)
#PBS -lnodes=1:ppn=6
###########################################################################
#### add modules:
. /etc/profile.d/modules.sh
export MODULEPATH=/apps/packages/bio/modulefiles:$MODULEPATH
module load bio/cufflinks/2.2.1
module list
#end of add modules
###########################################################################
###print PBS script
PBS_script="/groups/DAC/useful_PBS/PBS_cufflinks.sh"
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
echo changing to working directory
cd $work_DIR
pwd
echo hostname
hostname
hostname >$JOBID.node_track
df -h >>$JOBID.node_track
hostname >/groups/DAC/job_history/$JOBID.node_track
df -h >>/groups/DAC/job_history/$JOBID.node_track
##########################################################################
## begin execution stage # Below here enter the commands to start your job
#print the starting time of the job
echo start:
date

echo running Cufflinks on $BAM using $GTF

BAM_NAME=${BAM##*/}
BAM_DIR=${BAM%/*}
GTF_FILE=${GTF##*/}
GTF_DIR=${GTF%/*}
GENOME_FILE=${GENOME##*/}
GENOME_DIR=${GENOME%/*}


echo getting bam file from $BAM_DIR/$BAM_NAME
echo working in $node_DIR

echo cp $BAM_DIR/$BAM_NAME $node_DIR
cp -v $BAM_DIR/$BAM_NAME $node_DIR
echo cp $GTF_DIR/$GTF_FILE $node_DIR
cp -v $GTF_DIR/$GTF_FILE $node_DIR
echo cp $GENOME_DIR/$GENOME_FILE $node_DIR
cp -v $GENOME_DIR/$GENOME_FILE $node_DIR


cd $node_DIR
echo files in node_DIR
pwd
ls -lah
####run cufflinks

CLO_DIR=$BAM_NAME$EXT.CLO

lib_TYPE="fr-unstranded"
if [ $STRANDED = "second" ]; then
	echo library is stranded and set to 'second', correct for GRL, 
	lib_TYPE="fr-secondstrand"
	CLO_DIR=$CLO_DIR.SS
fi
if [ $STRANDED = "reverse" ]; then
	echo library is stranded and set to reverse, correct for Illumina UDP
	lib_TYPE="fr-firststrand"
	CLO_DIR=$CLO_DIR.FS
fi

pwd
ls -lah
echo
echo

##use -g to get novels, use -G to use GTF file as is
echo running cufflinks -p 6 --overlap-radius 1 -q --library-type $lib_TYPE -b $GENOME_FILE -N -g $GTF_FILE -o $CLO_DIR $BAM_NAME
cufflinks -p 6 -q --library-type $lib_TYPE -b $GENOME_FILE -N -g $GTF_FILE -o $CLO_DIR $BAM_NAME 
##rm accepted_hits.sam
cp -v -R $CLO_DIR/ $work_DIR/

#print the end time of the job
echo end:
date

echo ./$PBS_JOBID.out
rm -r $node_DIR
rm $work_DIR/$JOBID.node_track
rm /groups/DAC/job_history/$JOBID.node_track
echo $FILE
echo finished!
exit
