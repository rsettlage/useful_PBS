#!/bin/bash

###############################################################################
###      @author: Bob Settlage                                                                                          ###
###        Data Analysis Core @ the Virginia Bioinformatics Institute                                                   ###
###        December 2011                                                                                                ###
###Launch in target directory                                                                                           ###
###need reference                                                                                                       ###
###need file to map                                                                                                     ###
###align paired reads independently                    ###
###need to 
###        use absolute path for files                                                                                  ###
#### usage for SE/PE:  qsub /groups/DAC/useful_PBS/PBS_Popoolation_2_align.sh -v Results=_good,Reference=transcripts.fa,Reads=reads.fastq                               ###
####
#### usage find `pwd` -maxdepth 1 -iname '*fastq_[12]' -exec qsub /groups/DAC/useful_PBS/PBS_Popoolation_2_align.sh -v Results=_indi,Reference=/groups/DAC/Igor_Jan2014/Anopheles-stephensi-Indian_SCAFFOLDS_AsteI2.fa,Reads='{}' \;    ###
###                                                                                                                     ###
###############################################################################

####### job customization
## name our process
##        <----------------------------------------------set the next line
#PBS -N Pop_paki_cat
## merge stdout and stderr
#PBS -j oe
## e-mail us when the job: '-M b' = begins, '-M a' = aborts, '-M e' = ends
#PBS -m a -M rsettlage@vbi.vt.edu   
#PBS -lwalltime=15:00:00
################## Access group and queue, use one or the other#######max per node sfx=12/, smps are 40/#####################
#PBS -W group_list=sfx
#PBS -q sandybridge_q
#PBS -lnodes=1:ppn=12
####### end of job customization
export MODULEPATH=/apps/packages/bio/modulefiles:$MODULEPATH
export MODULEPATH=/apps/modulefiles/stats:$MODULEPATH
module load bio/bwa/0.7.4
module load pigz
###print PBS script
PBS_script="/groups/DAC/useful_PBS/PBS_Popoolation_2_align.sh"
echo '#############################################################################'
more $PBS_script
echo '#############################################################################'
###########################################################################
echo start:
date
echo jobid is $PBS_JOBID
JOBID=${PBS_JOBID%%.*}
echo jobnumber is $JOBID
echo job was launched using:
echo Type set to $Type 
echo Reference set to $Reference 
echo Reads set to $Reads
###########################################################################
################set up the directories on node and for results#########

reference_FILE=${Reference##*/}
reference_DIR=${Reference%/*}
reads1_FILE=${Reads##*/}
reads1_DIR=${Reads%/*}
reads2_FILE=${Reads##*/}

work_DIR=$PBS_O_WORKDIR
results_DIR=$reads1_DIR
node_DIR=/localscratch/$JOBID
mkdir $node_DIR
echo originating directory is $work_DIR
echo node directory is $node_DIR
echo results directory is $results_DIR
echo
cd $work_DIR
pwd
hostname >$JOBID.txt
df -h >>$JOBID.txt
hostname >/groups/DAC/job_history/$JOBID.txt
df -h >/groups/DAC/job_history/$JOBID.txt
##########################################################################
## begin execution stage # Below here enter the commands to start your job

echo starting processing
echo reference file is $reference_FILE

echo reads 1 file is $reads1_FILE
if [ -e $reads1_DIR/$reads1_FILE ]; then
	echo found read1 file, proceed with script
	cp -v $reads1_DIR/$reads1_FILE $node_DIR/
	reads1_EXT=${reads1_FILE##*.}
	if [ ${reads1_EXT} == 'gz' ];then
		echo unzipping $node_DIR/$reads1_FILE
		pigz -v -d -p 12 $node_DIR/$reads1_FILE
		new_reads1_FILE=${reads1_FILE%.*}
		echo new reads1 file is $new_reads1_FILE
		reads1_FILE=$new_reads1_FILE
	fi

else
	echo did not find read1 file reads1_DIR/$reads1_FILE so aborting script
	exit
fi

cd $node_DIR
echo
echo
ls -la
echo 
echo


### if needed, index reference    ##################################index################################

suffix=".rbwt"
indexed_reference_FILE=$reference_FILE$suffix
echo checking for presence of $indexed_reference_FILE

if [ -e $reference_DIR/$indexed_reference_FILE ]; then #only index the reference file if it is not present
	echo reference file is already indexed, moving files to local node and going to align stage
	cp -v $reference_DIR/$reference_FILE* $node_DIR/
else
	echo moving reference over and indexing reference file
	cp -v $reference_DIR/$reference_FILE $node_DIR/
	echo bwa index $reference_FILE
	bwa index $reference_FILE
fi

### do bwa aln  -uses suffix array##################################align#################################

ls -lah
echo time to start aligning....
echo created trimmed files via Popoolation trim, use them
echo reads1_FILE is $reads1_FILE
output_suffix=".sai"
reads1_out_FILE=$reads1_FILE$Results$output_suffix
echo creating suffix align array for $reads1_FILE and putting it into $reads1_out_FILE
echo bwa aln -t 11 -l 100 -o 2 -d 12 -e 12 -n 0.01 $reference_FILE $reads1_FILE -f $reads1_out_FILE
bwa aln -t 11 -l 100 -o 2 -d 12 -e 12 -n 0.01 $reference_FILE $reads1_FILE -f $reads1_out_FILE
echo finished aln stage for reads1
echo

###########################copy the results back########################
echo cp $node_DIR/$reads1_out_FILE $work_DIR/  ####<---change this to go to results dir
cd $node_DIR/
ls -lah
cp -v $node_DIR/$reads1_out_FILE $work_DIR/
rm $node_DIR/ -r

echo finished
hostname
rm /groups/DAC/job_history/$JOBID.txt
exit
