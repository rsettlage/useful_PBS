#!/bin/bash
#usage qsub /groups/DAC/useful_PBS/generic_PBS.sh -v command=/groups/DAC/useful_perl/blah.pl,FILE1='{}',FILE2='{}'

####### job customization
## name our process
##        <----------------------------------------------set the next line
#PBS -N sam_sort
## merge stdout and stderr
#PBS -j oe
## e-mail us when the job: '-M b' = begins, '-M a' = aborts, '-M e' = ends
#PBS -m a -M rsettlage@vbi.vt.edu   
#PBS -lwalltime=5:00:00
################## Access group and queue, use one or the other############################
#PBS -W group_list=sfx
#PBS -q sandybridge_q
#PBS -lnodes=1:ppn=6
####### end of job customization
export MODULEPATH=/apps/packages/bio/modulefiles:$MODULEPATH
export MODULEPATH=/apps/modulefiles/stats:$MODULEPATH
module load bio/samtools/0.1.19
module load R
module load pigz
module list
###print PBS script
PBS_script="/groups/DAC/useful_PBS/generic_PBS.sh"
echo '#############################################################################'
more $PBS_script
more $command
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
results_DIR=$work_DIR/$JOBID
#####mkdir $results_DIR ####<------good thing to do
node_DIR=/localscratch/$JOBID
mkdir $node_DIR
echo originating directory is $work_DIR
echo node directory is $node_DIR
echo results directory is $results_DIR
echo
cd $work_DIR
pwd
##########################################################################
## begin execution stage # Below here enter the commands to start your job

##CMD="samtools view -bS $FILE1 -o $FILE1.bam"
CMD="samtools sort -@ 6 $FILE1 $FILE1.sorted_coords"
##CMD="R CMD BATCH -_AT_QT.fastq.gz -001_AT_QT.fastq.gz -$FILE /groups/DAC/useful_R/illumina_HiSeq_QC_args.R"
##CMD="R CMD BATCH -.fastq.gz -001.fastq.gz -$FILE /groups/DAC/useful_R/illumina_HiSeq_QC_args.R"
##CMD="perl /groups/DAC/useful_perl/match_pend_fastq_RES3.pl $FILE1 $FILE2"
##CMD="samtools view -@ 4 -q 20 -bS $FILE1 -o $FILE1.bam"

##CMD="perl /apps/packages/bio/popoolation/1.2.2/basic-pipeline/trim-fastq.pl --input1 $FILE1 --input2 $FILE2 --output 
##$FILE3 --quality-threshold 20 --min-length 50 --fastq-type sanger"
echo $CMD
${CMD}

###########################copy the results back########################
##echo cp $node_DIR/* $work_DIR/  ####<---change this to go to results dir
##cd $node_DIR/
##ls >test.txt
##cp -v $node_DIR/* $work_DIR/
rm $node_DIR/ -r

exit;
