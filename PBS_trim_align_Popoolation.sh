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
#### usage for SE:  qsub /groups/DAC/useful_PBS/PBS_trim_align_Popoolation.sh -v Type=SE,Results=_good,Reference=transcripts.fa,Reads=reads.fastq                               ###
#### usage for PE:  qsub /groups/DAC/useful_PBS/PBS_trim_align_Popoolation.sh -v Type=PE,Results=_M_23,Reference=/groups/DAC/Fogelgren/mm_ref_GRCm38_chr17.fa,Reads=/groups/DAC/Fogelgren/working_data/BrBr2-2_R1.fastq.gz                 ###
#### usage find `pwd` -maxdepth 1 -iname '*R1*.fastq.gz' -exec qsub /groups/DAC/useful_PBS/PBS_trim_align_Popoolation.sh -v Type=PE,Results=indi,Reference=/groups/DAC/Igor_Jan2014/Anopheles-stephensi-Indian_SCAFFOLDS_AsteI2.fa,Reads='{}' \;    ###
###                                                                                                                     ###
###############################################################################

####### job customization
## name our process
##        <----------------------------------------------set the next line
#PBS -N Pop_indi
## merge stdout and stderr
#PBS -j oe
## e-mail us when the job: '-M b' = begins, '-M a' = aborts, '-M e' = ends
#PBS -m a -M rsettlage@vbi.vt.edu   
#PBS -lwalltime=100:00:00
################## Access group and queue, use one or the other#######max per node sfx=12/, smps are 40/#####################
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
PBS_script="/groups/DAC/useful_PBS/PBS_trim_align_Popoolation.sh"
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
	READ1_indicator="_R1_"
	READ2_indicator="_R2_"
	COMBINED_indicator="_R12_"
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
reads3_FILE=$reads1_FILE$trimmed_out

CMD="perl /apps/packages/bio/popoolation/1.2.2/basic-pipeline/trim-fastq.pl --input1 $reads1_FILE --input2 $reads2_FILE --output $reads3_FILE --quality-threshold 20 --min-length 50 --fastq-type sanger"
echo $CMD
${CMD}
ls -lah
ext_1="_1"
ext_2="_2"
reads1_FILE_trimmed=$reads3_FILE$ext_1
reads2_FILE_trimmed=$reads3_FILE$ext_2

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
echo reads1_FILE is $reads1_FILE_trimmed
output_suffix=".sai"
reads1_out_FILE=$reads1_FILE_trimmed$output_suffix
echo creating suffix align array for $reads1_FILE_trimmed
echo bwa aln -t 6 -l 100 -o 2 -d 12 -e 12 -n 0.01 $reference_FILE $reads1_FILE_trimmed -f $reads1_out_FILE
bwa aln -t 6 -l 100 -o 2 -d 12 -e 12 -n 0.01 $reference_FILE $reads1_FILE_trimmed -f $reads1_out_FILE
echo finished aln stage for reads1
echo

if [ $read_type = "PE" ]; then #if paired end, need to align mate as well
	echo creating suffix align array for $reads2_FILE_trimmed
	reads2_out_FILE=$reads2_FILE_trimmed$output_suffix
	echo bwa aln -t 6 -l 100 -o 2 -d 12 -e 12 -n 0.01 $reference_FILE $reads2_FILE_trimmed -f $reads2_out_FILE
	bwa aln -t 6 -l 100 -o 2 -d 12 -e 12 -n 0.01 $reference_FILE $reads2_FILE_trimmed -f $reads2_out_FILE
fi

### generate SAM file            #################################SAM#####################################
output_suffix=".sam"
if [ $read_type = "SE" ]; then #use either the samse or sampe branch for single or paired end data as appropriate
	reads1_out_sam_FILE=$reads1_FILE$Results$output_suffix
	echo generating SAM file for $reads1_FILE
	echo bwa samse -f $reads1_out_sam_FILE $reference_FILE $reads1_out_FILE $reads1_FILE
	bwa samse -f $reads1_out_sam_FILE $reference_FILE $reads1_out_FILE $reads1_FILE
else
	reads12_out_sam_FILE=$reads12_out_sam_FILE$Results$output_suffix
	echo generating SAM file for $reads1_FILE and $reads2_FILE for ouptut in $reads12_out_sam_FILE
	echo bwa sampe -f $reads12_out_sam_FILE $reference_FILE $reads1_out_FILE $reads2_out_FILE $reads1_FILE $reads2_FILE
	bwa sampe -f $reads12_out_sam_FILE $reference_FILE $reads1_out_FILE $reads2_out_FILE $reads1_FILE $reads2_FILE
fi

## Extract reads with a mapping quality of at least 20 (unambiguously mapped reads), create a sorted bam file and pileup
output_suffix=".mapped.sorted"
if [ $read_type = "SE" ]; then #use either the samse or sampe branch for single or paired end data as appropriate
	reads1_out_sorted_bam_FILE=$reads1_out_sam_FILE$output_suffix
	echo generating sorted BAM file for $reads1_FILE_trimmed
	echo samtools view -q 20 -bS $reads1_out_sam_FILE | samtools sort - $reads1_out_sorted_bam_FILE$Results
	samtools view -q 20 -bS $reads1_out_sam_FILE | samtools sort - $reads1_out_sorted_bam_FILE$Results
	echo samtools pileup $reads1_out_sorted_bam_FILE.bam >$reads1_out_sorted_bam_FILE$Results.bam.pileup
	samtools pileup $reads1_out_sorted_bam_FILE.bam >$reads1_out_sorted_bam_FILE$Results.bam.pileup
else
	reads12_out_sorted_bam_FILE=$reads12_out_sam_FILE$Results$output_suffix
	echo generating sorted BAM file for $reads1_FILE_trimmed and $reads2_FILE_trimmed for ouptut in $reads12_out_sam_FILE
	echo bwa sampe -n 30 -f $reads12_out_sam_FILE $reference_FILE $reads1_out_FILE $reads2_out_FILE $reads1_FILE $reads2_FILE
	echo samtools view -q 20 -bS $reads12_out_sam_FILE | samtools sort - reads12_out_sorted_bam_FILE$Results
	samtools view -q 20 -bS $reads12_out_sam_FILE | samtools sort - $reads12_out_sorted_bam_FILE$Results
	echo samtools pileup $reads12_out_sorted_bam_FILE.bam >$reads12_out_sorted_bam_FILE$Results.bam.pileup
	samtools pileup $reads12_out_sorted_bam_FILE.bam >$reads12_out_sorted_bam_FILE$Results.bam.pileup
fi

###########################copy the results back########################
echo cp $node_DIR/*bam $work_DIR/  ####<---change this to go to results dir
cd $node_DIR/
ls -lah
cp -v $node_DIR/*bam $work_DIR/
cp -v $node_DIR/*pileup $work_DIR/
rm $node_DIR/ -r

echo finished
exit
