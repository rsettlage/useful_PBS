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
#### usage for SE/PE:  qsub /groups/DAC/useful_PBS/PBS_Popoolation_6_popoolation2_stuff.sh -v EXT=all_indi_cat,GTF=/groups/DAC/Igor_Jan2014/Anopheles-stephensi-SDA-500_BASEFEATURES_AsteS1.0.gtf                              ###
#### 
#### MUST SUPPLY READS IN EXPERIMENT DESIGN FORMAT
###                                                                                                                     ###
###############################################################################

####### job customization
## name our process
##        <----------------------------------------------set the next line
#PBS -N Pop2_indi
## merge stdout and stderr
#PBS -j oe
## e-mail us when the job: '-M b' = begins, '-M a' = aborts, '-M e' = ends
#PBS -m a -M rsettlage@vbi.vt.edu   
#PBS -lwalltime=100:00:00
################## Access group and queue, use one or the other#######max per node sfx=12/, smps are 40/#####################
#PBS -W group_list=sfx
#PBS -q sandybridge_q
#PBS -lnodes=1:ppn=16
####### end of job customization
export MODULEPATH=/apps/packages/bio/modulefiles:$MODULEPATH
export MODULEPATH=/apps/modulefiles/stats:$MODULEPATH
module load bio/samtools/0.1.19
module load bio/popoolation2/1.201
module load pigz
###print PBS script
PBS_script="/groups/DAC/useful_PBS/PBS_Popoolation_6_popoolation_stuff.sh"
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
reference_FILE=${GTF##*/}
reference_DIR=${GTF%/*}


###must specify all files pakividually, copy all pakividually, etc
###assuming these to be sorted!!!
bam_DIR="/groups/DAC/Igor_Jan2014/working_data_Pop_trimmed/zz_split/zz_indi_split"
bam_FILE1="BAN.paki_cat_merged.bam.mapped.sorted.bam"
bam_FILE2="CHB.paki_cat_merged.bam.mapped.sorted.bam"
bam_FILE3="IRN.paki_cat_merged.bam.mapped.sorted.bam"
bam_FILE4="KAZ.paki_cat_merged.bam.mapped.sorted.bam"
bam_FILE5="WBAN.paki_cat_merged.bam.mapped.sorted.bam"

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
	echo did not find read1 file pileup_DIR/$pileup_FILE so aborting script
	exit
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


### now get the data, ie bam files  ##############################################################################
echo reads 1 file is $pileup_FILE
if [ -e $bam_DIR/$bam_FILE1 ]; then
	echo found read1 file, proceed with script
	cp -v $bam_DIR/$bam_FILE1 $node_DIR/
	cp -v $bam_DIR/$bam_FILE2 $node_DIR/
	cp -v $bam_DIR/$bam_FILE3 $node_DIR/
	cp -v $bam_DIR/$bam_FILE4 $node_DIR/
	cp -v $bam_DIR/$bam_FILE5 $node_DIR/
else
	echo did not find read1 file $bam_DIR/$bam_FILE1 so aborting script
	exit
fi


echo
echo
ls -lah
echo 
echo

###create mpileup with all data
samtools mpileup -B $bam_FILE1 $bam_FILE2 $bam_FILE3 $bam_FILE4 $bam_FILE5 >p1_p2_p3_p4_p5$EXT.mpileup

echo starting Popoolation2 stuff
ls -lah
##now do Popoolation2 stuff
echo java -ea -Xmx7g -jar  /apps/packages/bio/popoolation2/1.201/mpileup2sync.jar --input p1_p2_p3_p4_p5$EXT.mpileup --output p1_p2_p3_p4_p5$EXT_java.sync --fastq-type sanger --min-qual 20 --threads 16
java -ea -Xmx7g -jar  /apps/packages/bio/popoolation2/1.201/mpileup2sync.jar --input p1_p2_p3_p4_p5$EXT.mpileup --output p1_p2_p3_p4_p5$EXT.java.sync --fastq-type sanger --min-qual 20 --threads 16
##/apps/packages/bio/popoolation2/1.201/create-genewise-sync.pl --input p1_p2_p3_p4_p5$EXT.java.sync --gtf $reference_FILE --output pp1_p2_p3_p4_p5$EXT_genes.sync

##allele frequency difference
/apps/packages/bio/popoolation2/1.201/snp-frequency-diff.pl --input p1_p2_p3_p4_p5$EXT.java.sync --output-prefix pp1_p2_p3_p4_p5$EXT --min-count 6 --min-coverage 10 --max-coverage 200

##fst values
/apps/packages/bio/popoolation2/1.201/fst-sliding.pl --input p1_p2_p3_p4_p5$EXT.java.sync --output pp1_p2_p3_p4_p5$EXT.fst --suppress-noninformative --min-count 6 --min-coverage 50 --max-coverage 200 --min-covered-fraction 1 --window-size 1 --step-size 1 --pool-size 500
/apps/packages/bio/popoolation2/1.201/fst-sliding.pl --input p1_p2_p3_p4_p5$EXT.java.sync --output p1_p2_p3_p4_p5$EXT.w10000.fst --min-count 6 --min-coverage 50 --max-coverage 200 --min-covered-fraction 1 --window-size 10000 --step-size 10000 --pool-size 50
##/apps/packages/bio/popoolation2/1.201/fst-sliding.pl --min-count 6 --min-coverage 50 --max-coverage 200 --pool-size 500 --min-covered-fraction 0.0 --window-size 1000000 --step-size 1000000 --input p1_p2_p3_p4_p5$EXT.genes.sync --output p1_p2_p3_p4_p5$EXT.genewise.fst
/apps/packages/bio/popoolation2/1.201/export/pwc2igv.pl --input p1_p2_p3_p4_p5$EXT.fst --output p1_p2_p3_p4_p5$EXT.igv

##fishers exact test
/apps/packages/bio/popoolation2/1.201/fisher-test.pl --input p1_p2_p3_p4_p5$EXT.java.sync --output p1_p2_p3_p4_p5$EXT.fet --min-count 6 --min-coverage 50 --max-coverage 200 --suppress-noninformative
/apps/packages/bio/popoolation2/1.201/export/pwc2igv.pl --input p1_p2_p3_p4_p5$EXT.fst --output p1_p2_p3_p4_p5$EXT.fet.igv

###########################copy the results back########################
echo cp $node_DIR/*bam $work_DIR/  ####<---change this to go to results dir
cd $node_DIR/
ls -lah
cp -v $node_DIR/*sync $work_DIR/
cp -v $node_DIR/*fst $work_DIR/
cp -v $node_DIR/*rc $work_DIR/
cp -v $node_DIR/*pwc $work_DIR
cp -v $node_DIR/*fet $work_DIR
cp -v $node_DIR/*igv $work_DIR

rm $node_DIR/ -r
rm /groups/DAC/job_history/$JOBID.txt
date
echo finished
exit
