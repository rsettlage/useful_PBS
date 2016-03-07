#!/bin/bash

###use find `pwd` -maxdepth 1 -iname '*accepted_hits.bam' -exec qsub /groups/DAC/useful_PBS/PBS_HTSeqEXONS_FILE.sh -v STRANDED=yes,GFF=,EXT=,BAM='{}' \;
###note, send in the merged bam, script will sort and convert to sam with the headers <-----might need to work on the headers part for more complex cases
###
###NOTE ----------need to fix up the gf file with dexseq_prepare_annotation.py in the R package for DEXSeq??
###NOTE2----------use python /groups/DAC/useful_programs/dexseq_prepare_annotation.py and dexseq_count.py
###########################################################################
## environment & variable setup
####### job customization
## name our process
#PBS -N GV_DEX
## merge stdout and stderr
#PBS -j oe
## e-mail us when the job: '-M b' = begins, '-M a' = aborts, '-M e' = ends
#PBS -m a -M rsettlage@vbi.vt.edu
# Access group and queue
#PBS -W group_list=sfx
#PBS -q sandybridge_q
# Set the walltime, which is the maximum time your job can run in HH:MM:SS
#PBS -lwalltime=5:00:00
# Set the number of nodes, and the number of processors per node (up to 12)
#PBS -lnodes=1:ppn=3
####### end of job customization
# end of environment & variable setup
###########################################################################
#### add modules:
. /etc/profile.d/modules.sh
export MODULEPATH=/apps/packages/bio/modulefiles:$MODULEPATH
module load bio/samtools/0.1.19
module load bio/HTSeq/0.6.1p1
module list
#end of add modules
###########################################################################
###print PBS script
PBS_script="/groups/DAC/useful_PBS/PBS_HTSeqEXONS_FILE.sh"
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
hostname >$JOBID.node_track
df -h >>$JOBID.node_track
hostname >/groups/DAC/job_history/$JOBID.node_track
df -h >>/groups/DAC/job_history/$JOBID.node_track
##########################################################################
## begin execution stage # Below here enter the commands to start your job
#print the starting time of the job
echo start:
date

echo running HTseq on $BAM using $GFF

BAM_NAME=${BAM##*/}
BAM_DIR=${BAM%/*}
GFF_FILE=${GFF##*/}
GFF_DIR=${GFF%/*}


echo getting bam file from $BAM_DIR/$BAM_NAME
echo working in $node_DIR

echo cp $BAM_DIR/$BAM_NAME $node_DIR
cp -v $BAM_DIR/$BAM_NAME $node_DIR
echo cp $GFF_DIR/$GFF_FILE $node_DIR
cp -v $GFF_DIR/$GFF_FILE $node_DIR

cd $node_DIR
echo files in node_DIR
pwd
ls -lah

lib_TYPE="no"
LT="NS"
if [ $STRANDED = "yes" ]; then
	echo library is stranded and set to 'yes', correct for GRL
	lib_TYPE="yes"
	LT="YS"
fi
if [ $STRANDED = "reverse" ]; then
	echo library is stranded and set to reverse, correct for Illumina UDP
	lib_TYPE="reverse"
	LT="RS"
fi

#dont need to convert to sam anymore, default from TH is name sort
#BAM_EXT=${BAM_NAME##*.}
#if [ ${BAM_EXT} == 'bam' ];then
#	echo this is a bam, need to do a few extra things
#	####first sort bam, then convert it making sure the headers are included
#	###samtools sort -n -@6 889_merged.bam 889_merged.sorted
#	echo sorting bamfile using samtools sort -n -@3 $BAM_NAME $BAM_NAME.sorted
#	samtools sort -@3 $BAM_NAME $BAM_NAME.sorted
#	###samtools view -h 889_merged.sorted.bam -o 889_merged.sorted.bam.sam
#	echo converting bam to sam using samtools view -h $BAM_NAME.sorted.bam -o $BAM_NAME.sorted.bam.sam
#	samtools view -h $BAM_NAME.sorted.bam -o $BAM_NAME.sorted.bam.sam
#	SAM_NAME=$BAM_NAME.sorted.bam.sam
#	OUT_NAME=$BAM_NAME.sorted.bam.sam.EXONS.$LT.$EXT.out
#fi

#if [ ${BAM_EXT} == 'sam' ];then
#	echo this is a sam, need to get the name right
#	echo assuming sam is already sorted and has header, might need to check this...
#	SAM_NAME=$BAM_NAME
#	OUT_NAME=$BAM_NAME.EXONS.$LT.$EXT.out
#fi

OUT_NAME=$BAM_NAME.EXONS.$LT.$EXT.out

###need to first create a gff file using the scripts provided in the DEXSeq R package dexseq_prepare_annotation.py <in.gtf> <out.gff>
###python ~/R/library/DEXSeq/python_scripts/dexseq_count.py --paired yes --stranded reverse /groups/DAC/blastdbs/GRCh37/Annotation/GR37_genes_DEXseq.gff 889_merged.bam.sam 889_merged.bam.sam.out 
###use -s reverse for UDP method, use yes, no, reverse
echo running DEXSeq count script using python ~/R/library/DEXSeq/python_scripts/dexseq_count.py --paired=yes --stranded=$lib_TYPE $GFF_FILE $BAM_NAME $OUT_NAME

python /groups/DAC/useful_programs/dexseq_count.py -f bam -r name --paired=no --stranded=$lib_TYPE $GFF_FILE $BAM_NAME $OUT_NAME

#####get the data back...

echo files in local dir
ls -l
cp -v ./$OUT_NAME $BAM_DIR/

rm -r $node_DIR

#print the end time of the job
echo end:
date

echo ./$PBS_JOBID.out

echo $FILE
echo finished!

rm $work_DIR/$JOBID.node_track
rm /groups/DAC/job_history/$JOBID.node_track
exit;
