#!/bin/bash

###use find `pwd` -iname '*fastq.gz' -exec qsub /groups/DAC/useful_PBS/PBS_bowtie2_miRNA_FILE.sh -v EXT=EC81,INDEX=/groups/DAC/blastdbs/GRCh37/Sequence/Bowtie2Index/genome,FILE1='{}' \;
###note, for index, just pass it the base name and it will grab them all
###########################################################################
## environment & variable setup
####### job customization
## name our process
#PBS -N bt_mir
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
module load bio/bowtie2/2.2.4
module load samtools
module list
#end of add modules
###########################################################################
###print PBS script
PBS_script="/groups/DAC/useful_PBS/PBS_bowtie2_miRNA_FILE.sh"
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
##if [ $PAIRED = "Y" ]; then   <-----------------------------------paired reads doesn't make sense
##	echo is paired, so get file2
##	READ1_indicator="_R1_"
##	READ2_indicator="_R2_"
##	FILE_NAME2=${FILE_NAME1/"$READ1_indicator"/"$READ2_indicator"}
##	echo cp $FILE_DIR1/$FILE_NAME2 $node_DIR
##	cp -v $FILE_DIR1/$FILE_NAME2 $node_DIR
##	align_FILES=$FILE_NAME1" "$FILE_NAME2
##fi

cd $node_DIR
echo files in node_DIR
pwd
ls -lah

echo aligning using bowtie2
echo bowtie2 -L 14 -p 6 -x $INDEX_BASE -U $align_FILES -S $align_FILES.sam
###optimized for short sequences???  
###Need to check, but did bowtie2 -L 14 -p 6 -x $INDEX_BASE -U $align_FILES -S $align_FILES.sam for Ina
bowtie2 -L 14 -p 6 -x $INDEX_BASE -U $align_FILES -S $align_FILES.sam

echo
echo doing some sam manipulation
`awk '{ if($12~/AS:i:/ || NF<5) print }' $FILE_NAME1.sam > $FILE_NAME1.sam_aligned`
`awk '{ if($5!="0" || NF<5) print }' $FILE_NAME1.sam_aligned > $FILE_NAME1.sam_nomulti`

##convert to bam
echo samtools view -@ 6 -bS $align_FILES.sam \>$align_FILES.$EXT.sam.bam
#echo samtools view -@ 6 -bS $FILE_NAME1.sam_aligned \>$FILE_NAME1.$EXT.sam.aligned.bam
echo samtools view -@ 6 -bS $FILE_NAME1.sam_nomulti \>$FILE_NAME1.$EXT.sam.nomulti.bam

samtools view -@ 6 -bS $align_FILES.sam >$align_FILES.$EXT.sam.bam
#samtools view -@ 6 -bS $FILE_NAME1.sam_aligned >$FILE_NAME1.$EXT.sam.aligned.bam
samtools view -@ 6 -bS $FILE_NAME1.sam_nomulti >$FILE_NAME1.$EXT.sam.nomulti.bam

#####get the data back...
echo 
echo files in local dir
ls -l
cp -v ./*.bam $FILE_DIR1/
#cp -v ./$align_FILES.* $FILE_DIR1/

rm -r $node_DIR

#print the end time of the job
echo end:
date

echo ./$PBS_JOBID.out

echo $FILE
echo finished!

rm $work_DIR/$JOBID.node_track
exit;
