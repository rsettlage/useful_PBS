#!/bin/bash

### qsub /groups/DAC/useful_PBS/PBS_fastqmcf_trim_R1R2_pairs.sh   # searches for all matched _R1_ and _R2_ fastq files (or fastq.gz)
###note, FASTQ can be gzipped or not
###########################################################################
## environment & variable setup
####### job customization
## name our process
#PBS -N FastqMcf_R1R2
## merge stdout and stderr
#PBS -j oe
## e-mail us when the job: '-M b' = begins, '-M a' = aborts, '-M e' = ends
###PBS -m a -M rsettlage@vbi.vt.edu
# Access group and queue
#PBS -W group_list=sfx
#PBS -q sandybridge_q
# Set the walltime, which is the maximum time your job can run in HH:MM:SS
#PBS -lwalltime=5:00:00
# Set the number of nodes, and the number of processors per node (up to 12)
#PBS -lnodes=1:ppn=1
####### end of job customization
# end of environment & variable setup
###########################################################################
#### add modules:
module load bio/ea-utils
module list
#end of add modules

###########################################################################
###print PBS script
PBS_script=$0
echo '#############################################################################'
#cat $PBS_script
echo '#############################################################################'
###########################################################################
echo start:
date
echo jobid is $PBS_JOBID
JOBID=${PBS_JOBID%%.*}
echo jobnumber is $JOBID
###########################################################################
################set up the directories on node and for results#########
# just stay in current directory, fast enough not to need local drive
pwd

##########################################################################
## begin execution stage # Below here enter the commands to start your job
#print the starting time of the job
echo start:
date
ADAPTER_FILE="/groups/DAC/vector_contaminant_fasta/Illumina_adapters.txt"

###for f in $PBS_O_WORKDIR/*_R1_*fastq $PBS_O_WORKDIR/*_R1_*fastq.gz
for f in $PBS_O_WORKDIR/GRL4553_9-20-13-BR2_CTTGTA_L008_R1.fastq $PBS_O_WORKDIR/*R1_*fastq.gz

do
	echo file = $f
	#ls -l $f
	o1="${f%.fastq.gz}_fastqmcf.fastq"
	f2=`perl -e '$f=shift; print $f if $f=~s/_R1_/_R2_/' $f`
	echo file2 = $f2
	if [ -f $f2 ]; then
		#ls -l $f2
		o2="${f2%.fastq.gz}_fastqmcf.fastq"
		fastq-mcf $ADAPTER_FILE  $f $f2 -o $o1 -o $o2
		echo "Trimmed output to $o1 and $o2"
	else
		echo "Could not find mate file $f2"
	fi
done
#print the end time of the job
date
echo finished!
exit

