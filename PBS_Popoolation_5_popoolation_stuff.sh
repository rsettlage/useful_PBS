#!/bin/bash

###############################################################################
###      @author: Bob Settlage                                                                                          ###
###        Data Analysis Core @ the Virginia Bioinformatics Institute                                                   ###
###        December 2011                                                                                                ###
###Launch in target directory                                                                                           ###
###need reference                                                                                                       ###
###need file to map                                                                                                     ###
###
###need to 
###        use absolute path for files                                                                                  ###
#### usage for SE/PE:  qsub /groups/DAC/useful_PBS/PBS_Popoolation_5_popoolation_stuff.sh -v GTF=transcripts.fa,Pileup=reads.fastq                               ###
#### 
#### usage find `pwd` -maxdepth 1 -iname '*pileup' -exec qsub /groups/DAC/useful_PBS/PBS_Popoolation_5_popoolation_stuff.sh -v GTF=/groups/DAC/Igor_Jan2014/Anopheles-stephensi-,Pileup='{}' \;    ###
###                                                                                                                     ###
###############################################################################

####### job customization
## name our process
##        <----------------------------------------------set the next line
#PBS -N Pop_paki
## merge stdout and stderr
#PBS -j oe
## e-mail us when the job: '-M b' = begins, '-M a' = aborts, '-M e' = ends
#PBS -m a -M rsettlage@vbi.vt.edu   
#PBS -lwalltime=19:00:00
################## Access group and queue, use one or the other#######max per node sfx=12/, smps are 40/#####################
#PBS -W group_list=sfx
#PBS -q sandybridge_q
#PBS -lnodes=1:ppn=6
####### end of job customization
export MODULEPATH=/apps/packages/bio/modulefiles:$MODULEPATH
export MODULEPATH=/apps/modulefiles/stats:$MODULEPATH
module load bio/popoolation/1.2.2
module load pigz
###print PBS script
PBS_script="/groups/DAC/useful_PBS/PBS_Popoolation_5_popoolation_stuff.sh"
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
reference_FILE=${GTF##*/}
reference_DIR=${GTF%/*}
pileup_FILE=${Pileup##*/}
pileup_DIR=${Pileup%/*}

work_DIR=$PBS_O_WORKDIR
results_DIR=$pileup_DIR
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
### get reference GTF, if it is a GFF, convert it    ##################################Reference GTF################################

echo getting/converting GTF/GFF

if [ -e $reference_DIR/$reference_FILE ]; then #make sure it is there, then get it
	cp -v $reference_DIR/$reference_FILE $node_DIR/
else
	echo did not find read1 file pileup_DIR/$pileup_FILE so not going to do variation by position
	echo give a gtf file next time
	#exit
fi

reference_EXT=${reference_FILE##*.}
if [ ${reference_EXT} == 'gff' ];then
	echo converting gff to gtf
	`cat $reference_FILE | awk '$2=="FlyBase" && $3=="exon"'| perl -pe 's/ID=([^:;]+)([^;]+)?;.*/gene_id "$1"; transcript_id "$1:1";/'> $reference_FILE.gtf`
	new_reference_FILE=${reference_FILE%.*}.GTF
	echo new reference file is $new_reference_FILE
	reference_FILE=$new_reference_FILE
else
	echo extracting exons from gtf file
	`awk '$3=="exon"' $reference_FILE >$reference_FILE.exons`
	new_reference_FILE=$reference_FILE.exons
	reference_FILE=$new_reference_FILE
fi


### now get the data, ie Pileup file  ##############################################################################
echo reads 1 file is $pileup_FILE
if [ -e $pileup_DIR/$pileup_FILE ]; then
	echo found read1 file, proceed with script
	cp -v $pileup_DIR/$pileup_FILE $node_DIR/
	pileup_EXT=${pileup_FILE##*.}
else
	echo did not find read1 file $pileup_DIR/$pileup_FILE so aborting script
	exit
fi


echo
echo
ls -la
echo 
echo


##now do Popoolation stuff

##Variance-sliding.pl
pi=$pileup_FILE.sliding.pi1
D=$pileup_FILE.sliding.D1
echo perl /apps/packages/bio/popoolation/1.2.2/Variance-sliding.pl --input $pileup_FILE --output $pi --measure pi --window-size 10000 --step-size 10000 --min-count 2 --min-coverage 4 --max-coverage 400 --min-qual 20 --pool-size 20 --min-covered-fraction 0.2
perl /apps/packages/bio/popoolation/1.2.2/Variance-sliding.pl --input $pileup_FILE --output $pi --measure pi --fastq-type sanger --window-size 10000 --step-size 10000 --min-count 2 --min-coverage 4 --max-coverage 400 --min-qual 20 --pool-size 20 --min-covered-fraction 0.2
cp -v $node_DIR/*pi $work_DIR/

echo perl /apps/packages/bio/popoolation/1.2.2/Variance-sliding.pl --input $pileup_FILE --output $D --measure D --window-size 10000 --step-size 10000 --min-count 2 --min-coverage 4 --max-coverage 400 --min-qual 20 --pool-size 20 --min-covered-fraction 0.2
perl /apps/packages/bio/popoolation/1.2.2/Variance-sliding.pl --input $pileup_FILE --output $D --measure D --fastq-type sanger --window-size 10000 --step-size 10000 --min-count 2 --min-coverage 4 --max-coverage 400 --min-qual 20 --pool-size 20 --min-covered-fraction 0.2
cp -v $node_DIR/*D $work_DIR/

##Visualise output of Variance-sliding.pl  -- need to provide a chr or scaffold to do this for
perl /apps/packages/bio/popoolation/1.2.2/Visualise-output.pl --input $pi --output $pi.pdf --ylab pi --chromosomes "Indian_2L Indian_2R Indian_3L Indian_3R Indian_X Indian_unknown"
perl /apps/packages/bio/popoolation/1.2.2/Visualise-output.pl --input $D --output $D.pdf --ylab D --chromosomes "Indian_2L Indian_2R Indian_3L Indian_3R Indian_X Indian_unknown"

cp -v $node_DIR/*pdf $work_DIR/

perl /apps/packages/bio/popoolation/1.2.2/VarSliding2Wiggle.pl --input $pi --trackname "dmel Pi" --output $pi.wig
perl /apps/packages/bio/popoolation/1.2.2/VarSliding2Wiggle.pl --input $D --trackname "dmel D" --output $D.wig
cp -v $node_DIR/*wig $work_DIR/

##Variance-at-position
pi=$pileup_FILE.position.pi
D=$pileup_FILE.position.D
echo perl /apps/packages/bio/popoolation/1.2.2/Variance-at-position.pl --pool-size 500 --min-qual 20 --min-coverage 4 --min-count 2 --max-coverage 4 --pileup $pileup_FILE --gtf $reference_FILE --output $pi --measure pi --min-covered-fraction 0.2
perl /apps/packages/bio/popoolation/1.2.2/Variance-at-position.pl --pool-size 20 --min-qual 20 --min-coverage 4 --min-count 2 --max-coverage 4 --pileup $pileup_FILE --gtf $reference_FILE --output $pi --measure pi --min-covered-fraction 0.2 --fastq-type sanger

echo perl /apps/packages/bio/popoolation/1.2.2/Variance-at-position.pl --pool-size 500 --min-qual 20 --min-coverage 4 --min-count 2 --max-coverage 4 --pileup $pileup_FILE --gtf $reference_FILE --output $D --measure D --min-covered-fraction 0.2
perl /apps/packages/bio/popoolation/1.2.2/Variance-at-position.pl --pool-size 20 --min-qual 20 --min-coverage 4 --min-count 2 --max-coverage 4 --pileup $pileup_FILE --gtf $reference_FILE --output $D --measure D --min-covered-fraction 0.2 --fastq-type sanger



###########################copy the results back########################
cd $node_DIR/
ls -lah
cp -v $node_DIR/$pi $work_DIR/
cp -v $node_DIR/$D $work_DIR/
cp -v $node_DIR/*params $work_DIR/
cp -v $node_DIR/*pdf $work_DIR/
cp -v $node_DIR/*wig $work_DIR/
rm $node_DIR/ -r

echo finished
exit
