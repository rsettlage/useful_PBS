###############################################################################
###      @author: Bob Settlage                                                                                          ###
###        Data Analysis Core @ the Virginia Bioinformatics Institute                                                   ###
###        December 2011                                                                                                ###
###Launch in target directory                                                                                           ###
###combining adapter trimming, qual trimming and if needed re-pairing pairs
###
###use ---absolute paths---- in anticipation of making this into a PBS script
###find `pwd` -iname '*gz' -exec qsub /groups/DAC/useful_PBS/PBS_FASTQ_QC_control_EA_UTILS.sh -v Type=s,FILE1='{}' \;
###find `pwd` -iname '*_R1_*gz' -exec qsub /groups/DAC/useful_PBS/PBS_FASTQ_QC_control_EA_UTILS.sh -v Type=p,FILE1='{}' \;
###qsub /groups/DAC/useful_PBS/PBS_FASTQ_QC_control_EA_UTILS.sh -v Type=p,SHUFFLE=no,FILE1=/path/file1.fastq
###qsub /groups/DAC/useful_PBS/PBS_FASTQ_QC_control_EA_UTILS.sh -v Type=p,SHUFFLE=no,FILE1='{}' \;
###                                                                                                                     ###
###############################################################################

####### job customization
## name our process
##        <----------------------------------------------set the next line
#PBS -N Trim_mgDNA
## merge stdout and stderr
#PBS -j oe
## e-mail us when the job: '-M b' = begins, '-M a' = aborts, '-M e' = ends
#PBS -m a -M rsettlage@vbi.vt.edu   
#PBS -lwalltime=40:00:00
################## Access group and queue, use one or the other#######max per node sfx=12/, smps are 40/#####################
#PBS -W group_list=sfx
#PBS -q sandybridge_q
#PBS -lnodes=1:ppn=16 
####### end of job customization
export MODULEPATH=/apps/packages/bio/modulefiles:$MODULEPATH
export MODULEPATH=/apps/modulefiles/stats:$MODULEPATH
module load pigz
module load bio/ea-utils
###print PBS script
PBS_script="/groups/DAC/useful_PBS/PBS_FASTQ_QC_control_EA_UTILS.sh"
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
hostname >$JOBID.node_track
df -h >>$JOBID.node_track
##########################################################################
## begin execution stage # Below here enter the commands to start your job

quality_EXT="_QT"
adaptor_EXT="_AT"
summary_EXT=".sum.txt"
ADAPTER_FILE="/groups/DAC/vector_contaminant_fasta/Illumina_adapters.txt"
fastq_EXT=".fastq"

read_FILE1=${FILE1##*/}
read_DIR1=${FILE1%/*}
cp -v $read_DIR1/$read_FILE1 $node_DIR
cd $node_DIR
pwd
ls -lah

read_FILE1_EXT=${read_FILE1##*.}
if [ ${read_FILE1_EXT} == 'gz' ];then
	echo unzipping $read_FILE1
	pigz -p 16 -v -d $read_FILE1
	new_read_FILE1=${read_FILE1%.*}
	echo new read file is $new_read_FILE1
	read_FILE1=$new_read_FILE1
fi


temp1=${read_FILE1%.*}
adaptor_trimmed_FILE1=$temp1$adaptor_EXT$fastq_EXT
quality_trimmed_FILE1=$temp1$adaptor_EXT$quality_EXT$fastq_EXT
quality_summary_FILE1=$temp1$summary_EXT

echo this is $Type data
echo starting processing on Read 1
echo putting output in
echo $adaptor_trimmed_FILE1
echo $quality_trimmed_FILE1
echo $quality_summary_FILE1
echo $final_out1
echo

##first remove adaptor, then trim by qual

	echo running adapter trimming
	fastq-mcf -l 50 -U -t 0 -C 1000000 $ADAPTER_FILE $node_DIR/$read_FILE1 -o $node_DIR/$adaptor_trimmed_FILE1
	echo starting the quality trimming
	code="/groups/DAC/useful_programs/Btrim/Btrim64 -t $node_DIR/$adaptor_trimmed_FILE1 -o $node_DIR/$quality_trimmed_FILE1 -s $node_DIR/$quality_summary_FILE1 -q -l 15 -w 5 -a 15 -S"
	${code}
	rm $quality_summary_FILE1
	perl /groups/DAC/useful_perl/removeOptHead.pl $node_DIR/$quality_trimmed_FILE1
	noOpt_EXT="_noOpt.fastq"
	mv -v $node_DIR/$quality_trimmed_FILE1$noOpt_EXT $node_DIR/$quality_trimmed_FILE1

if [ $Type = "p" ]; then
	echo starting R2 processing
	READ1_indicator="_R1_"
	READ2_indicator="_R2_"
	read_FILE2=${read_FILE1/"$READ1_indicator"/"$READ2_indicator"}
	echo altering $read_FILE1 to give $read_FILE2
	read_DIR2=$read_DIR1
	echo read file2 is $read_FILE2
	if [ ${read_FILE1_EXT} == 'gz' ];then
		cp -v $read_DIR2/$read_FILE2.gz $node_DIR/
		echo unzipping $read_FILE2
		pigz -v -p 16 -d $node_DIR/$read_FILE2
		echo unzipped read file is $read_FILE2
		wasZipped="1"
	else
		cp -v $read_DIR2/$read_FILE2 $node_DIR/
	fi
	
	pwd
	ls -lah
	temp2=${read_FILE2%.*}
	echo basefile is $temp2
	adaptor_trimmed_FILE2=$temp2$adaptor_EXT$fastq_EXT
	quality_trimmed_FILE2=$temp2$adaptor_EXT$quality_EXT$fastq_EXT
	quality_summary_FILE2=$temp2$summary_EXT
	echo
	echo putting read2 output in
	echo $adaptor_trimmed_FILE2
	echo $quality_trimmed_FILE2
	echo $quality_summary_FILE2
	echo
	fastq-mcf -l 50 -U -t 0 -C 1000000 $ADAPTER_FILE $node_DIR/$read_FILE2 -o $node_DIR/$adaptor_trimmed_FILE2
	code="/groups/DAC/useful_programs/Btrim/Btrim64 -t $node_DIR/$adaptor_trimmed_FILE2 -o $node_DIR/$quality_trimmed_FILE2 -s $node_DIR/$quality_summary_FILE2 -q -l 15 -w 5 -a 15 -S >test2.txt"
	echo $code
	${code}
	rm $quality_summary_FILE2
	perl /groups/DAC/useful_perl/removeOptHead.pl $node_DIR/$quality_trimmed_FILE2
	noOpt_EXT="_noOpt.fastq"
	mv -v $node_DIR/$quality_trimmed_FILE2$noOpt_EXT $node_DIR/$quality_trimmed_FILE2
	code="perl /groups/DAC/useful_perl/match_pend_fastq_RES3.pl $node_DIR/$quality_trimmed_FILE1 $node_DIR/$quality_trimmed_FILE2"
	echo $code
	${code}	
#	if [ $SHUFFLE = "yes" ]; then
#		matched_EXT=".paired_matched.fastq"
#		shuffled_EXT=".paired_matched.shuffled.fastq"
#		paired_FILE1=${quality_trimmed_FILE1%.*}
#		paired_FILE1=$paired_FILE1$matched_EXT
#		paired_FILE2=${quality_trimmed_FILE2%.*}
#		paired_FILE2=$paired_FILE2$matched_EXT
#		shuffled_FILE=${quality_trimmed_FILE1%.*}
#		shuffled_FILE=$shuffled_FILE$shuffled_EXT
#		code="perl /groups/DAC/useful_perl/shuffleSequences_fastq.pl $node_DIR/$paired_FILE1 $node_DIR/$paired_FILE2 $node_DIR/$shuffled_FILE"
#		echo $code
#		${code}
#	fi
fi

#if [ $wasZipped = "1" ]; then
	ls -lah
	echo zipping $quality_trimmed_FILE1
	wc -l *AT_QT*fastq
	pigz -v -p 16 $node_DIR/*AT_QT*fastq
	cp -v $node_DIR/$quality_trimmed_FILE1.gz $read_DIR1
	###cp -v $node_DIR/* $read_DIR1
	ls -l $read_DIR1/$quality_trimmed_FILE1.gz
	if [ $Type = "p" ]; then
		cp -v $node_DIR/$quality_trimmed_FILE2.gz $read_DIR1
		##ls -l $read_DIR1/$quality_trimmed_FILE2.gz
		cp -v $node_DIR/*paired_matched* $read_DIR1
		cp -v $node_DIR/*shuffled* $read_DIR1
	fi
#fi

rm -rf $node_DIR
rm $work_DIR/$JOBID.node_track
echo finished!
exit 1
