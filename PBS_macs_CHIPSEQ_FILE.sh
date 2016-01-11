#!/bin/bash

###use find `pwd` -maxdepth 1 -iname '*.bam' -exec qsub /groups/DAC/useful_PBS/PBS_macs_CHIPSEQ_FILE.sh -v PAIRED=N,TYPE=narrow,FILT=Y,NAME=zz_macs2_test,FILE1='{}' \;
###note, for index, just pass it the base name and it will grab them all
###########################################################################
## environment & variable setup
####### job customization
## name our process
#PBS -N macs
## merge stdout and stderr
#PBS -j oe
## e-mail us when the job: '-M b' = begins, '-M a' = aborts, '-M e' = ends
#PBS -m a -M rsettlage@vbi.vt.edu
# Access group and queue
#PBS -W group_list=sfx
#PBS -q sandybridge_q
# Set the walltime, which is the maximum time your job can run in HH:MM:SS
#PBS -lwalltime=30:00:00
# Set the number of nodes, and the number of processors per node (up to 12)
#PBS -lnodes=1:ppn=12
####### end of job customization
# end of environment & variable setup
###########################################################################
#### add modules:
. /etc/profile.d/modules.sh
export MODULEPATH=/apps/packages/bio/modulefiles:$MODULEPATH
module load site/shadowfax/easybuild/setup
module load MACS2
module load samtools
module list
#end of add modules
###########################################################################
###print PBS script
PBS_script="/groups/DAC/useful_PBS/PBS_macs_CHIPSEQ_FILE.sh"
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

echo running macs2 on $FILE1 

FILE_NAME1=${FILE1##*/}
FILE_DIR1=${FILE1%/*}
outputName=$NAME


echo getting read file from $FILE_DIR1/$FILE_NAME1
echo working in $node_DIR

echo get $FILE_NAME1
echo cp $FILE_DIR1/$FILE_NAME1 $node_DIR
cp -v $FILE_DIR1/$FILE_NAME1 $node_DIR


align_FILE=$FILE_NAME1
if [ $PAIRED = "Y" ]; then
	echo is paired, so get file2
	READ1_indicator="_R1_"
	READ2_indicator="_R2_"
	FILE_NAME2=${FILE_NAME1/"$READ1_indicator"/"$READ2_indicator"}
	echo cp $FILE_DIR1/$FILE_NAME2 $node_DIR
	cp -v $FILE_DIR1/$FILE_NAME2 $node_DIR
	align_FILES="-1 "$FILE_NAME1" -2 "$FILE_NAME2
fi

cd $node_DIR
echo files in node_DIR
pwd
ls -lah

###filter out unmapped
macs_FILE=$align_FILE
if [ $FILT = "Y" ]; then
	echo samtools view -hb -@ 12 -q 10 $align_FILE -o $align_FILE.q10.bam
	samtools view -hb -@ 12 -q 10 $align_FILE -o $align_FILE.q10.bam
	macs_FILE=$align_FILE.q10.bam
fi
###

###run macs2
if [ $TYPE = "narrow" ]; then
	macs2 --version
	echo macs2 callpeak -t $macs_FILE --name=$macs_FILE."_".$outputName --gsize 2.7e9 --keep-dup all 
	macs2 callpeak -t $macs_FILE --name=$macs_FILE."_".$outputName --gsize 2.7e9 --keep-dup all
fi

if [ $TYPE = "broad" ]; then
	macs2 --version
	echo macs2 callpeak -t $macs_FILE --name=$macs_FILE."_".$outputName --gsize 2.7e9 --keep-dup all --nomodel --broad --nolambda
	macs2 callpeak -t $macs_FILE --name=$macs_FILE."_".$outputName --gsize 2.7e9 --keep-dup all --nomodel --broad --nolambda
fi

#####get the data back...

echo files in local dir
ls -l
rm *sam
rm *bam
cp -v *$outputName* $FILE_DIR1/

rm -r $node_DIR

#print the end time of the job
echo end:
date

echo ./$PBS_JOBID.out

echo $FILE
echo finished!

rm $work_DIR/$JOBID.node_track
exit;
