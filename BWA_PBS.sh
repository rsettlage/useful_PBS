#!/bin/bash

###############################################################################
###      @author: Bob Settlage                                                                                          ###
###        Data Analysis Core @ the Virginia Bioinformatics Institute                                                   ###
###        December 2011                                                                                                ###
###Launch in target directory                                                                                           ###
###need reference                                                                                                       ###
###need file to map                                                                                                     ###
###need SE flag, 1 = SE, 2 = PE  --for PE, expecting read indicator in file name to be '_R1_"                    ###
###need to 
###        use absolute path for files                                                                                  ###
#### usage for SE:  qsub /groups/DAC/useful_PBS/BWA_PBS.sh -v Type=SE,Results=_good,Reference=transcripts.fa,Reads=reads.fastq                               ###
#### usage for PE:  qsub /groups/DAC/useful_PBS/BWA_PBS.sh -v Type=PE,Results=_M_23,Reference=/groups/DAC/Fogelgren/mm_ref_GRCm38_chr17.fa,Reads=/groups/DAC/Fogelgren/working_data/BrBr2-2_R1.fastq.gz                 ###
#### usage find `pwd` -maxdepth 1 -iname '*R1*.fastq.gz' -exec qsub /groups/DAC/useful_PBS/BWA_PBS.sh -v Type=PE,Reference=transcripts.fa,Reads='{}' \;    ###
###                                                                                                                     ###
###############################################################################

####### job customization
## name our process
##        <----------------------------------------------set the next line
#PBS -N BWA_ND
## merge stdout and stderr
#PBS -j oe
## e-mail us when the job: '-M b' = begins, '-M a' = aborts, '-M e' = ends
#PBS -m a -M rsettlage@vbi.vt.edu   
#PBS -lwalltime=10:00:00
################## Access group and queue, use one or the other#######max per node sfx=12/, smps are 40/#####################
#PBS -W group_list=sfx
#PBS -q sandybridge_q
#PBS -lnodes=1:ppn=3
####### end of job customization
export MODULEPATH=/apps/packages/bio/modulefiles:$MODULEPATH
export MODULEPATH=/apps/modulefiles/stats:$MODULEPATH
module load bio/bwa/0.7.6a
module load pigz
###print PBS script
PBS_script="/groups/DAC/useful_PBS/BWA_PBS.sh"
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
read_type=$Type
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
		pigz -p 3 -v -d $node_DIR/$reads1_FILE
		new_reads1_FILE=${reads1_FILE%.*}
		echo new reads1 file is $new_reads1_FILE
		reads1_FILE=$new_reads1_FILE
	fi

else
	echo did not find read1 file reads1_DIR/$reads1_FILE so aborting script
	exit
fi
if [ $read_type = "PE" ]; then
	echo this is for paired end data
	echo reads 1 file is $reads1_FILE
	READ1_indicator="_R1"
	READ2_indicator="_R2"
	COMBINED_indicator="_R12"
	echo going to do substitution
	echo read1/2 indicators are $READ1_indicator $READ2_indicator
	echo file name we are changing $reads2_FILE
	reads2_FILE=${reads2_FILE/"$READ1_indicator"/"$READ2_indicator"}
	reads12_out_sam_FILE=${reads1_FILE/"$READ1_indicator"/"$COMBINED_indicator"}
	echo reads 2 file is $reads2_FILE 
	echo sam file is $reads12_out_sam_FILE
	if [ -e $reads1_DIR/$reads2_FILE ]; then
		echo found read2 file as $reads2_FILE, proceed with script
		cp -v $reads1_DIR/$reads2_FILE $node_DIR/
		reads2_EXT=${reads2_FILE##*.}
		if [ ${reads2_EXT} = 'gz' ];then
			echo unzipping $node_DIR/$reads2_FILE
			pigz -p 3 -v -d $node_DIR/$reads2_FILE
			new_reads2_FILE=${reads2_FILE%.*}
			echo new reads2 file is $new_reads2_FILE
			reads2_FILE=$new_reads2_FILE
		fi
	else
		echo did not find read2 file $reads1_DIR/$reads2_FILE so aborting script
		echo files on node
		ls -lah $node_DIR/
		echo 
		echo
		echo files in read dir
		ls -lah $reads1_DIR/
		exit
	fi
fi

echo
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
	cd $node_DIR/
	echo bwa index -p index $reference_FILE
	bwa index -p index $reference_FILE
fi

### do bwa mem  -uses suffix array##################################align#################################

cd $node_DIR/
ls -lah
echo time to start aligning....
echo
output_suffix=".sam"

RG_string=`echo $reads1_FILE | cut -d"_" -f1,2`
echo $RG_string

RG_header="@RG\tID:"$RG_string"\tLB:"$RG_string"\tSM:"$RG_string"\tPL:ILLUMINA"
echo $RG_header

if [ $read_type = "PE" ]; then #if paired end, need to align mate as well
	echo reads1_FILE is $reads1_FILE
	echo reads2_FILE is $reads2_FILE
	echo doing alignment using $reference_FILE
	reads2_out_FILE=$reads2_FILE$output_suffix
	`echo bwa mem -t 2 index $reads1_FILE $reads2_FILE > $reads2_out_FILE`
	`bwa mem -t 2 index $reads1_FILE $reads2_FILE > $reads2_out_FILE`
else
	echo reads1_FILE is $reads1_FILE
	reads1_out_FILE=$reads1_FILE$output_suffix
	echo doing alignment using $reference_FILE
	`echo bwa mem -t 3 -R $RG_header index $reads1_FILE > $reads1_out_FILE`
	`bwa mem -t 3 -R $RG_header index $reads1_FILE > $reads1_out_FILE`
	echo finished aln stage for reads1
fi

### generate SAM file            #################################SAM#####################################
 ##output_suffix=".sam"
##if [ $read_type = "SE" ]; then #use either the samse or sampe branch for single or paired end data as appropriate
##	reads1_out_sam_FILE=$reads1_FILE$Results$output_suffix
##	echo generating SAM file for $reads1_FILE
##	echo bwa samse -f $reads1_out_sam_FILE $reference_FILE $reads1_out_FILE $reads1_FILE
##	bwa samse -f $reads1_out_sam_FILE $reference_FILE $reads1_out_FILE $reads1_FILE
##else
##	reads12_out_sam_FILE=$reads12_out_sam_FILE$Results$output_suffix
##	echo generating SAM file for $reads1_FILE and $reads2_FILE for ouptut in $reads12_out_sam_FILE
##	echo bwa sampe -n 30 -f $reads12_out_sam_FILE $reference_FILE $reads1_out_FILE $reads2_out_FILE $reads1_FILE $reads2_FILE
##	bwa sampe -n 30 -P -f $reads12_out_sam_FILE $reference_FILE $reads1_out_FILE $reads2_out_FILE $reads1_FILE $reads2_FILE
##fi

###########################copy the results back########################
echo cp $node_DIR/* $work_DIR/  ####<---change this to go to results dir
cd $node_DIR/
ls -lah
cp -v $node_DIR/*sam $work_DIR/
rm $node_DIR/ -r

rm $work_DIR/$JOBID.txt
rm /groups/DAC/job_history/$JOBID.txt

echo finished
exit
