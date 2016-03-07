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
#### usage for SE:  qsub /groups/DAC/useful_PBS/PBS_Popoolation_1_trim.sh -v Type=SE,Reads=reads.fastq                               ###
#### usage for PE:  qsub /groups/DAC/useful_PBS/PBS_Popoolation_1_trim.sh -v Type=PE,Reads=/groups/DAC/Fogelgren/working_data/BrBr2-2_R1.fastq.gz                 ###
#### usage find `pwd` -maxdepth 1 -iname '*_1.fastq.gz' -exec qsub /groups/DAC/useful_PBS/PBS_Popoolation_1_trim.sh -v Type=PE,Reads='{}' \;    ###
###                                                                                                                     ###
###############################################################################

####### job customization
## name our process
##        <----------------------------------------------set the next line
#PBS -N Pop_trim
## merge stdout and stderr
#PBS -j oe
## e-mail us when the job: '-M b' = begins, '-M a' = aborts, '-M e' = ends
#PBS -m a -M rsettlage@vbi.vt.edu   
#PBS -lwalltime=100:00:00
################## Access group and queue, use one or the other#######max per node sfx=12/, smps are 40/#####################
## this process doesn't parallelize well, but due to the disk IO, keep it at about 6 ppn, change the pigz if ppn changes
#PBS -W group_list=sfx
#PBS -q sandybridge_q
#PBS -lnodes=1:ppn=6
####### end of job customization
export MODULEPATH=/apps/packages/bio/modulefiles:$MODULEPATH
export MODULEPATH=/apps/modulefiles/stats:$MODULEPATH
module load bwa
bio/samtools/0.1.18
bio/popoolation/1.2.2
module load pigz
###print PBS script
PBS_script="/groups/DAC/useful_PBS/PBS_Popoolation_1_trim.sh"
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
		pigz -v -d -p 6 $node_DIR/$reads1_FILE
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
	READ1_indicator="_1"
	READ2_indicator="_2"
	COMBINED_indicator="_R12_"
	echo read1/2 indicators are $READ1_indicator $READ2_indicator
	reads2_FILE=${reads2_FILE/"$READ1_indicator"/"$READ2_indicator"}
	reads_out_FILE=${reads1_FILE/"$READ1_indicator"/"$COMBINED_indicator"}
	echo reads 2 file is $reads2_FILE 
	echo sam file is $reads12_out_sam_FILE
	if [ -e $reads1_DIR/$reads2_FILE ]; then
		echo found read2 file as $reads2_FILE, proceed with script
		cp -v $reads1_DIR/$reads2_FILE $node_DIR/
		reads2_EXT=${reads2_FILE##*.}
		if [ ${reads2_EXT} = 'gz' ];then
			echo unzipping $node_DIR/$reads2_FILE
			pigz -v -d -p 6 $node_DIR/$reads2_FILE
			new_reads2_FILE=${reads2_FILE%.*}
			echo new reads2 file is $new_reads2_FILE
			reads2_FILE=$new_reads2_FILE
		fi
	else
		echo did not find read2 file $reads1_DIR/$reads2_FILE so aborting script
		exit
	fi
fi

echo
echo
ls -la
echo 
echo

### do trimming using the Popoolation trim script
cd $node_DIR/

trimmed_out=".trimmed"
reads_out_FILE=$reads1_FILE$trimmed_out

CMD="perl /apps/packages/bio/popoolation/1.2.2/basic-pipeline/trim-fastq.pl --input1 $reads1_FILE --input2 $reads2_FILE --output $reads_out_FILE --quality-threshold 20 --min-length 50 --fastq-type sanger"
echo $CMD
${CMD}
ls -lah

if [ ${reads2_EXT} = 'gz' ];then
	echo zipping output files(s) -- dont worry about error if this wasn't PE
	pigz -v -p 6 $node_DIR/$reads_out_FILE*
fi

###########################copy the results back########################
echo cp $node_DIR/*_[1/2] $work_DIR/  
cd $node_DIR/
ls -lah
cp -v $node_DIR/$reads_out_FILE* $work_DIR/
rm $node_DIR/ -r

echo finished
exit
