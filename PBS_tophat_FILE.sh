#!/bin/bash

###use find `pwd` -iname '*_R1_*fastq.gz' -exec qsub /groups/DAC/useful_PBS/PBS_tophat_FILE.sh -v EXT=,PAIRED=Y,STRANDED=second,FILE1='{}',GTF=,INDEX= \;
###note, for index, just pass it the base name and it will grab them all
###for Sai's stranded method, use second strand for TH/Bowtie and stranded=first for HTSeq
###########################################################################
## environment & variable setup
####### job customization
## name our process
#PBS -N aln_TH
## merge stdout and stderr
#PBS -j oe
## e-mail us when the job: '-M b' = begins, '-M a' = aborts, '-M e' = ends
#PBS -m a -M rsettlage@vbi.vt.edu
# Access group and queue
#PBS -W group_list=sfx
#PBS -q sandybridge_q
# Set the walltime, which is the maximum time your job can run in HH:MM:SS
#PBS -lwalltime=10:00:00
# Set the number of nodes, and the number of processors per node (up to 12)
#PBS -lnodes=1:ppn=6
####### end of job customization
# end of environment & variable setup
###########################################################################
#### add modules:
. /etc/profile.d/modules.sh
export MODULEPATH=/apps/packages/bio/modulefiles:$MODULEPATH
module load bio/tophat2/2.0.13
module list
#end of add modules
###########################################################################
###print PBS script
PBS_script="/groups/DAC/useful_PBS/PBS_tophat_FILE.sh"
echo '#############################################################################'
more $PBS_script
echo '#############################################################################'
###########################################################################
echo start:
date
echo jobid is $PBS_JOBID
JOBID=${PBS_JOBID%%.*}
echo jobnumber is $JOBID
echo hostname is:
hostname
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
##########################################################################
## begin execution stage # Below here enter the commands to start your job
#print the starting time of the job
echo start:
date

echo running tophat on $FILE1 using $INDEX

FILE_NAME1=${FILE1##*/}
FILE_DIR1=${FILE1%/*}
INDEX_BASE=${INDEX##*/}
INDEX_DIR=${INDEX%/*}
GTF_FILE=${GTF##*/}
GTF_DIR=${GTF%/*}


echo getting read file from $FILE_DIR1/$FILE_NAME1
echo working in $node_DIR

echo cp $FILE_DIR1/$FILE_NAME1 $node_DIR
cp -v $FILE_DIR1/$FILE_NAME1 $node_DIR
echo cp $INDEX_DIR/$INDEX_BASE* $node_DIR
cp -v $INDEX_DIR/$INDEX_BASE* $node_DIR
echo cp $GTF_DIR/$GTF_FILE $node_DIR
cp -v $GTF_DIR/$GTF_FILE $node_DIR

align_FILES=$FILE_NAME1
if [ $PAIRED = "Y" ]; then
	echo is paired, so get file2
	READ1_indicator="_R1"
	READ2_indicator="_R2"
	FILE_NAME2=${FILE_NAME1/"$READ1_indicator"/"$READ2_indicator"}
	echo cp $FILE_DIR1/$FILE_NAME2 $node_DIR
	cp -v $FILE_DIR1/$FILE_NAME2 $node_DIR
	align_FILES=$FILE_NAME1" "$FILE_NAME2
fi

lib_TYPE="fr-unstranded"
out_dir_suffix="_US_THO2"
if [ $STRANDED = "second" ]; then
	echo library is stranded and set to 'second', correct for GRL, 
	lib_TYPE="fr-secondstrand"
	out_dir_suffix="_SS_THO2"
fi
if [ $STRANDED = "first" ]; then
	echo library is stranded and set to 'first', correct for Illumina UDP
	lib_TYPE="fr-firststrand"
	out_dir_suffix="_FS_THO2"
fi

cd $node_DIR
echo files in node_DIR
pwd
ls -lah

####first look to see if it is zipped, if so, unzip it and set the extension approprately
###can take zipped files now, should test using pigz as the unzipper as TH can take alternative unzippers

out_DIR=$FILE_NAME1$EXT$out_dir_suffix
echo output to $out_DIR

## for Ina: tophat2 -p 6 -r 100 --library-type fr-firststrand -G $GTF_FILE -o $out_DIR $INDEX_BASE $align_FILES  
##tophat2 -p 6 -r 100 --library-type fr-firststrand -G $GTF_FILE -o $out_DIR $INDEX_BASE $align_FILES  
##tophat2 -p 6 --library-type fr-unstranded -G $GTF_FILE -o $out_DIR $INDEX_BASE $align_FILES  
###some confusion here on first vs second strand, for SS sequencing with UDP, seems to work best as fr-firststrand and then use reverse in the HTseq part
###according to Illumina, this is correct.  For normal libs, use fr-unstranded
##for just alignment, not considering the gtf, use:
echo tophat2 -p 6 --library-type $lib_TYPE -o ./ -G $GTF_FILE $INDEX_BASE $align_FILES
#tophat2 -p 6 --library-type $lib_TYPE -o ./ -G $GTF_FILE $INDEX_BASE $align_FILES

##no GTF file
tophat2 -p 6 --library-type $lib_TYPE -o ./ $INDEX_BASE $align_FILES
###trying to map to the transcriptome
##do quick more of the align_summary.txt file
more .//align_summary.txt


#####get the data back...

echo files in local dir
ls -l
cp -v ./unmapped.bam $FILE_DIR1/$out_DIR."unmapped.bam"
cp -v ./accepted_hits.bam $FILE_DIR1/$out_DIR."accepted_hits.bam"

rm -r $node_DIR

#print the end time of the job
echo end:
date

echo ./$PBS_JOBID.out

echo $FILE
echo finished!

rm $work_DIR/$JOBID.node_track
exit;
