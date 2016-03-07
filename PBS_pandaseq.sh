###############################################################################
###      @author: Bob Settlage                                                                                          ###
###        Data Analysis Core @ the Virginia Bioinformatics Institute                                                   ###
###        December 2011                                                                                                ###
###Launch in target directory                                                                                           ###
###run pandaseq, need to supply read 1 file, it does rest
###
###use ---absolute paths---- in anticipation of making this into a PBS script
###find `pwd` -iname '*R1*gz' -exec qsub /groups/DAC/useful_PBS/PBS_pandaseq.sh -v FILE1='{}' \;
###                                                                                                                     ###
###############################################################################

####### job customization
## name our process
##        <----------------------------------------------set the next line
#PBS -N PS_PM
## merge stdout and stderr
#PBS -j oe
## e-mail us when the job: '-M b' = begins, '-M a' = aborts, '-M e' = ends
#PBS -m a -M rsettlage@vbi.vt.edu   
#PBS -lwalltime=1:00:00
################## Access group and queue, use one or the other#######max per node sfx=12/, smps are 40/#####################
#PBS -W group_list=sfx
#PBS -q sandybridge_q
#PBS -lnodes=1:ppn=6
####### end of job customization
export MODULEPATH=/apps/packages/bio/modulefiles:$MODULEPATH
export MODULEPATH=/apps/modulefiles/stats:$MODULEPATH
module load bio/pandaseq/2013-05-07
###print PBS script
PBS_script="/groups/DAC/useful_PBS/PBS_pandaseq.sh"
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
echo
cd $work_DIR
pwd
hostname >$JOBID.node_track
df -h >>$JOBID.node_track
##########################################################################
## begin execution stage # Below here enter the commands to start your job

overlapped_EXT=".overlapped.fasta"
singles_EXT=".not_overlapped.fasta"


fastq_EXT=".fastq"

temp1=$FILE1
read_FILE1=${temp1##*/}
read_DIR1=${temp1%/*}
cp -v $read_DIR1/$read_FILE1 $node_DIR

READ1_indicator="_R1_"
READ2_indicator="_R2_"
read_FILE2=${read_FILE1/"$READ1_indicator"/"$READ2_indicator"}
echo altering $read_FILE1 to give $read_FILE2
read_DIR2=$read_DIR1
cp -v $read_DIR2/$read_FILE2 $node_DIR

cd $node_DIR
pwd
ls -lah

read_EXT=${read_FILE1##*.}

temp1=${read_FILE1%.*}
overlapped_FILE1=$temp1$overlapped_EXT
singles_FILE1=$temp1$singles_EXT


echo starting processing 
echo putting output in
echo $overlapped_FILE1
echo $singles_FILE1
echo

##output is fasta
##pandaseq -f $read_FILE1 -r $read_FILE2 -T 6 -o 10 -B >$overlapped_FILE1 -u $singles_FILE1
##output is fastq, use -o 10 to require 10 base overlap
pandaseq -f $read_FILE1 -r $read_FILE2 -T 6 -B -F >$overlapped_FILE1 -u $singles_FILE1

cp -v $node_DIR/$overlapped_FILE1 $read_DIR1
cp -v $node_DIR/$singles_FILE1 $read_DIR1

rm -rf $node_DIR
rm $work_DIR/$JOBID.node_track
echo finished!
exit 1
