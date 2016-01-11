#!/bin/bash


### qsub /groups/DAC/useful_PBS/PBS_RSEM_trinity_assembly_FILE.sh -v PAIRED=Y,db_FILE=/groups/DAC/Liu_Jan2014/assemblies_trinity/Trinity_all.fasta,FILE1=/groups/DAC/Liu_Jan2014/working_data/M2F_R1.fastq
####qsub /groups/DAC/useful_PBS/PBS_RSEM_trinity_assembly_FILE.sh -v PAIRED=Y,db_FILE=/groups/DAC/Mattison_March2014/working_data/Pecan_Trinity_SS/Pecan_Trinity_SS.fasta,FILE1=/groups/DAC/Mattison_March2014/working_data/GRL4548_8-11-12-BR3_CAGATC_L008_R1.fastq
####### job customization
## name our process
##        <----------------------------------------------set the next line
#PBS -N RSEM_GRL3884_MID
## merge stdout and stderr
#PBS -j oe
## e-mail us when the job: '-M b' = begins, '-M a' = aborts, '-M e' = ends
#PBS -m a -M rsettlage@vbi.vt.edu
# Access group and queue
#PBS -W group_list=sfx
#PBS -q sfx_q
# Set the walltime, which is the maximum time your job can run in HH:MM:SS
#PBS -lwalltime=100:00:00
# Set the number of nodes, and the number of processors per node (up to 12)
#PBS -lnodes=1:ppn=12
####### end of job customization
export MODULEPATH=/apps/packages/bio/modulefiles:$MODULEPATH
export MODULEPATH=/apps/modulefiles/stats:$MODULEPATH
module load bio/RSEM/1.2.7
module load bio/trinity/2013-11-10
###print PBS script
PBS_script="/groups/DAC/useful_PBS/PBS_RSEM_trinity_assembly_FILE.sh"
echo #############################################################################
more $PBS_script
echo #############################################################################
JOBID=${PBS_JOBID%%.*}
echo jobnumber is $JOBID
node_DIR="/localscratch/$JOBID"
mkdir $node_DIR
pwd
echo #############################################################################
#print the starting time of the job
echo start:
date

cd $PBS_O_WORKDIR

reference_FILE=${db_FILE##*/}
reference_DIR=${db_FILE%/*}

cp -v $reference_DIR/$reference_FILE $node_DIR/ 

FILE_NAME1=${FILE1##*/}
FILE_DIR=${FILE1%/*}
cp -v $FILE_DIR/$FILE_NAME1 $node_DIR/    ####<-------------set this

cd $node_DIR  

if [ $PAIRED = "Y" ]; then
	echo is paired, so get file2
	READ1_indicator="R1"
	READ2_indicator="R2"
	FILE_NAME2=${FILE_NAME1/"$READ1_indicator"/"$READ2_indicator"}
	echo cp $FILE_DIR/$FILE_NAME2 $node_DIR
	cp -v $FILE_DIR/$FILE_NAME2 $node_DIR
fi

echo /apps/packages/bio/trinity/2013-11-10/util/RSEM_util/run_RSEM_align_n_estimate.pl --transcripts $reference_FILE --seqType fq --left $FILE_NAME1 --right $FILE_NAME2  --thread_count 12
/apps/packages/bio/trinity/2013-11-10/util/RSEM_util/run_RSEM_align_n_estimate.pl --transcripts $reference_FILE --seqType fq --left $FILE_NAME1 --right $FILE_NAME2  --thread_count 12

rsem-plot-model $condition_1_NAME $condition_1_NAME.pdf

ls -lah
echo copying files
cp -v ./RSEM.isoforms.results $FILE_DIR/$FILE_NAME1.RSEM.isoforms.results
cp -v ./RSEM.genes.results $FILE_DIR/$FILE_NAME1.RSEM.gene.results
rm $node_DIR -r
echo files copied
cd $final_DIR
ls -lah

#print the end time of the job 
echo end:
date

echo ./$PBS_JOBID.out

exit;