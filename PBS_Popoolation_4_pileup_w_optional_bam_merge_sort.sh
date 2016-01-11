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
#### usage for SE:  qsub /groups/DAC/useful_PBS/PBS_Popoolation_4_pileup_w_optional_bam_merge_sort.sh -v Type=SE,Results=_good,Reference=transcripts.fa,BAM=reads.fastq                               ###
#### 
#### usage find `pwd` -maxdepth 1 -iname '*bam' -exec qsub /groups/DAC/useful_PBS/PBS_Popoolation_4_pileup_w_optional_bam_merge_sort.sh -v Reference=transcripts.fa,BAM='{}' \;    ###
###                                                                                                                     ###
###############################################################################

####### job customization
## name our process
##        <----------------------------------------------set the next line
#PBS -N Pop_pile_paki
## merge stdout and stderr
#PBS -j oe
## e-mail us when the job: '-M b' = begins, '-M a' = aborts, '-M e' = ends
#PBS -m a -M rsettlage@vbi.vt.edu   
#PBS -lwalltime=20:00:00
################## Access group and queue, use one or the other#######max per node sfx=12/, smps are 40/#####################
#PBS -W group_list=sfx
#PBS -q sandybridge_q
#PBS -lnodes=1:ppn=6
####### end of job customization
export MODULEPATH=/apps/packages/bio/modulefiles:$MODULEPATH
export MODULEPATH=/apps/modulefiles/stats:$MODULEPATH
bio/popoolation/1.2.2
module load pigz
###print PBS script
PBS_script="/groups/DAC/useful_PBS/PBS_Popoolation_4_pileup_w_optional_bam_merge_sort.sh"
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
bam_FILE=${BAM##*/}
bam_DIR=${BAM%/*}

work_DIR=$PBS_O_WORKDIR
results_DIR=$reads1_DIR
node_DIR=/localscratch/$JOBID
mkdir $node_DIR
echo originating directory is $work_DIR
echo node directory is $node_DIR
echo results directory is $results_DIR
echo
cd $node_DIR
pwd
hostname >$JOBID.txt
df -h >>$JOBID.txt
hostname >/groups/DAC/job_history/$JOBID.txt
df -h >/groups/DAC/job_history/$JOBID.txt
##########################################################################
## begin execution stage # Below here enter the commands to start your job

echo starting processing
echo reference file is $reference_FILE
###get reference    ##################################index################################

if [ -e $reference_DIR/$reference_FILE ]; then #foudn it
	echo found reference, cpying it over
	cp -v $reference_DIR/$reference_FILE $node_DIR/
else
	echo didnt find reference, abort
	exit
fi

### now get the data    ##############################################################################
echo bam file is $bam_FILE
if [ -e $bam_DIR/$bam_FILE ]; then
	echo found bam file, proceed with script
	cp -v $bam_DIR/$bam_FILE $node_DIR/
else
	echo did not find bam file $bam_DIR/$bam_FILE so aborting script
	exit
fi


echo
echo
ls -la
echo 
echo


## Extract reads with a mapping quality of at least 20 (unambiguously mapped reads), create a sorted bam file and pileup
output_suffix=".mapped.sorted"

	reads12_out_sorted_bam_FILE=$bam_FILE$Results$output_suffix
	echo generating sorted BAM file for $bam_FILE and $reads2_FILE for ouptut in $reads12_out_sorted_bam_FILE
	echo /groups/DAC/Novozymes_pipeline/samtools_0.1.16 view -q 20 -b $bam_FILE \| samtools sort - $reads12_out_sorted_bam_FILE
	`/groups/DAC/Novozymes_pipeline/samtools_0.1.16 view -q 20 -b $bam_FILE | /groups/DAC/Novozymes_pipeline/samtools_0.1.16 sort - $reads12_out_sorted_bam_FILE`

###make pileups

echo /groups/DAC/Novozymes_pipeline/samtools_0.1.16 pileup -f $reference_FILE $reads12_out_sorted_bam_FILE.bam >$reads12_out_sorted_bam_FILE.bam.pileup
`/groups/DAC/Novozymes_pipeline/samtools_0.1.16 pileup -f $reference_FILE $reads12_out_sorted_bam_FILE.bam >$reads12_out_sorted_bam_FILE.bam.pileup`


###########################copy the results back########################
echo cp $node_DIR/*bam $work_DIR/  ####<---change this to go to results dir
cd $node_DIR/
ls -lah
cp -v $node_DIR/*bam $work_DIR/
cp -v $node_DIR/*pileup $work_DIR/
rm $node_DIR/ -r

echo finished
exit
