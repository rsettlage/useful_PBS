#!/bin/bash

### for i in 71 77 81 87 91 ; do qsub /groups/DAC/useful_PBS/PBS_Velvethg_noCNY.sh -v KMER=$i ; done
####qsub /groups/DAC/useful_PBS/PBS_Velvethg_noCNY.sh -v KMER=91
####### job customization
## name our process
#PBS -N Vlt_968
## merge stdout and stderr
#PBS -j oe
## e-mail us when the job: '-M b' = begins, '-M a' = aborts, '-M e' = ends
#PBS -m a -M rsettlage@vbi.vt.edu
# Access group and queue
#PBS -W group_list=sfx
#PBS -q sandybridge_q
# Set the walltime, which is the maximum time your job can run in HH:MM:SS
#PBS -lwalltime=48:00:00
# Set the number of nodes, and the number of processors per node (up to 12)
#PBS -lnodes=1:ppn=12
####### end of job customization
echo /groups/DAC/useful_PBS/PBS_Velvethg_noCNY.sh
more /groups/DAC/useful_PBS/PBS_Velvethg_noCNY.sh
export MODULEPATH=/apps/packages/bio/modulefiles:$MODULEPATH
export MODULEPATH=/apps/modulefiles/stats:$MODULEPATH
module load bio/velvet/1.2.08-bin
module load amos
#print the starting time of the job
echo start:
date
####run Velvetg

kmer=$KMER

cd $PBS_O_WORKDIR
assembly_DIR=968_$KMER

shortPaired1="/groups/DAC/Morrill_Apr2015/working_data/968_S14_L001_R1_001_AT_QT.paired_matched.fastq.gz /groups/DAC/Morrill_Apr2015/working_data/968_S14_L001_R2_001_AT_QT.paired_matched.fastq.gz"

velveth $assembly_DIR $kmer -fastq.gz -separate -shortPaired $shortPaired1
#velveth ./ 81 -reuse_binary
cd $assembly_DIR
velvetg ./ -cov_cutoff auto -exp_cov auto -ins_length 200 -min_contig_lgth 200 -read_trkg yes -unused_reads yes -scaffolding yes -amos_file no
#velvetg ./ -cov_cutoff auto -exp_cov auto -min_contig_lgth 200 -read_trkg yes -unused_reads yes -amos_file yes

perl /groups/DAC/useful_perl/getContigSummary.pl contigs.fa -o contigs.fa.stats
##amos2ace velvet_asm.afg

#print the end time of the job
echo end:
date

echo ./$PBS_JOBID.out


exit
