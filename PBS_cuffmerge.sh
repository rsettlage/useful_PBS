#!/bin/bash

###usage qsub /groups/DAC/useful_PBS/PBS_cuffmerge.sh -v db_file=,gtf_file=,cufflinks_gtf_manifest_file=,out_dir=
####### job customization
## name our process
#PBS -N zz_AT_QT_CT
## merge stdout and stderr
#PBS -j oe
## e-mail us when the job: '-M b' = begins, '-M a' = aborts, '-M e' = ends
#PBS -m a -M rsettlage@vbi.vt.edu
# Access group and queue
#PBS -W group_list=sfx
#PBS -q sandybridge_q
#PBS -lwalltime=20:00:00
#PBS -lnodes=1:ppn=12
####PBS -lnodes=1:ppn=6:smp for sfxsmp nodes
####### end of job customization, adding modules
export MODULEPATH=/apps/packages/bio/modulefiles:$MODULEPATH
module load bio/cufflinks/2.2.1
###########
qsubfile=/groups/DAC/useful_PBS/PBS_cuffmerge.sh
more $qsubfile
###########################################################################
echo start:
date
echo jobid is $PBS_JOBID
JOBID=${PBS_JOBID%%.*}
echo jobnumber is $JOBID
###########################################################################
#print the starting time of the job
echo start:
date

####run cufflinks

cd $PBS_O_WORKDIR
mkdir $out_dir
cd $out_dir

##note that manifest file is simply a text listing, out_dir does not have the slash at the end
cuffmerge -p 12 -s $db_file -g $gtf_file $cufflinks_gtf_manifest_file -o $out_dir

#print the end time of the job
echo end:
date

echo ./$PBS_JOBID.out
exit;
