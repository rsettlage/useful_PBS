#!/bin/bash

###use find ./ \( -name '*fastq' \) -exec qsub /groups/DAC/useful_PBS/PBS_STAR_FILE.sh -v PAIRED=Y,genome_DIR=/groups/DAC/blastdbs/GRCh37/Sequence/WholeGenomeFasta/STAR_db/,FILE1='{}' \;
###note, for index, just pass it the base name and it will grab them all
###########################################################################
## environment & variable setup
####### job customization
## name our process
#PBS -N Ina_STAR
## merge stdout and stderr
#PBS -j oe
## e-mail us when the job: '-M b' = begins, '-M a' = aborts, '-M e' = ends
#PBS -m a -M rsettlage@vbi.vt.edu
# Access group and queue
#PBS -W group_list=sfx
#PBS -q sandybridge_q
# Set the walltime, which is the maximum time your job can run in HH:MM:SS
#PBS -lwalltime=1:00:00
# Set the number of nodes, and the number of processors per node (up to 12)
#PBS -lnodes=1:ppn=12
####### end of job customization
# end of environment & variable setup
###########################################################################
#### add modules:
. /etc/profile.d/modules.sh
export MODULEPATH=/apps/packages/bio/modulefiles:$MODULEPATH
module load test/bio/STAR/2.3.0
module load samtools
module load pigz
module list
#end of add modules
###########################################################################
###print PBS script
PBS_script="/groups/DAC/useful_PBS/PBS_STAR_FILE.sh"
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
hostname >>host_list.txt
##########################################################################
## begin execution stage # Below here enter the commands to start your job
#print the starting time of the job
echo start:
date

echo running STAR on $FILE1 using $INDEX

FILE_NAME1=${FILE1##*/}
FILE_DIR1=${FILE1%/*}
##INDEX_BASE=${INDEX##*/}
##genome_DIR=${INDEX%/*}

echo getting read file from $FILE_DIR1/$FILE_NAME1
echo working in $node_DIR
cd $node_DIR

##echo cp $FILE_DIR1/$FILE_NAME1 $node_DIR
##cp -v $FILE_DIR1/$FILE_NAME1 $node_DIR
##echo cp $INDEX_DIR/$INDEX_BASE* $node_DIR
##cp -v $INDEX_DIR/$INDEX_BASE* $node_DIR
##echo cp $GTF_DIR/$GTF_FILE $node_DIR
##cp -v $GTF_DIR/$GTF_FILE $node_DIR

align_FILES=$FILE_DIR1/$FILE_NAME1
query_EXT=${FILE_NAME1##*.}
if [ ${query_EXT} == 'gz' ];then
	new_FILE_NAME1=${FILE_NAME1%.*}
	echo new query file is $new_FILE_NAME1
	echo unzipping $FILE_NAME1 to local node $new_FILE_NAME1
	pigz -p 12 -v -c -d $FILE_DIR1/$FILE_NAME1 >$new_FILE_NAME1
	align_FILES=$new_FILE_NAME1
fi


if [ $PAIRED = "Y" ]; then
	echo is paired, so get file2
	READ1_indicator="_R1_"
	READ2_indicator="_R2_"
	FILE_NAME2=${FILE_NAME1/"$READ1_indicator"/"$READ2_indicator"}
	new_FILE_NAME2=${FILE_NAME2%.*}
	echo new query file is $new_FILE_NAME2
	##echo cp $FILE_DIR1/$FILE_NAME2 $node_DIR
	##cp -v $FILE_DIR1/$FILE_NAME2 $node_DIR
	align_FILES=$FILE_DIR1/$FILE_NAME1" "$FILE_DIR1/$FILE_NAME2
	if [ ${query_EXT} == 'gz' ];then
		echo unzipping $FILE_NAME2 to local node $new_FILE_NAME2
		pigz -p 12 -v -c -d $FILE_DIR1/$FILE_NAME2 >$new_FILE_NAME2
		align_FILES=$new_FILE_NAME1" "$new_FILE_NAME2
	fi
fi


echo files in node_DIR
pwd
ls -lah

####first look to see if it is zipped, if so, unzip it and set the extension approprately
###can take zipped files now, should test using pigz as the unzipper as TH can take alternative unzippers

out_dir_suffix="_STAR"
out_DIR=$FILE_NAME1$out_dir_suffix
echo output to $out_DIR

##run star

echo running with STAR --genomeDir $genome_DIR --runThreadN 12 --genomeLoad NoSharedMemory --outSAMmode Full --outFilterMismatchNmax 3 --outFilterMultimapNmax 5 --readFilesIn $align_FILES
STAR --genomeDir $genome_DIR --runThreadN 12 --genomeLoad NoSharedMemory --outFilterMismatchNmax 3 --outFilterMultimapNmax 5 --readFilesIn $align_FILES

##convert samt to bam to save space
samtools view -Sb -@12 Aligned.out.sam -o Aligned.out.sam.bam

#####get the data back...

echo files in local dir
ls -l
cp -v ./Aligned.out.sam.bam $FILE_DIR1/$FILE_NAME1.STAR_aligned.sam.bam
cp -v ./SJ.out.tab $FILE_DIR1/$FILE_NAME1.STAR_SJ.tab
cp -v ./Log.final.out $FILE_DIR1/$FILE_NAME1.STAR_stats.txt
cp -v ./Log.out $FILE_DIR1/$FILE_NAME1.STAR.log

rm -r $node_DIR

#print the end time of the job
echo end:
date

echo ./$PBS_JOBID.out

echo $FILE
echo finished!

rm $work_DIR/$JOBID.node_track
exit;
