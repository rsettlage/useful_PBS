#!/bin/bash

###USEAGE
##qsub /groups/DAC/useful_PBS/blat_SFX_PBS.sh -v SUBJECT=,QUERY=

# Set the walltime, which is the maximum time your job can run in HH:MM:SS
#PBS -lwalltime=20:00:00
#PBS -lnodes=1:ppn=1
#PBS -m a -M rsettlage@vbi.vt.edu
# Access group and queue
#PBS -W group_list=grl
#PBS -q grl_q
#PBS -N lrna_CA_all
#PBS -j oe


# If modules are needed, source modules environment:
. /etc/profile.d/modules.sh

# Add any modules you might require:

module load bio/blat

cd $PBS_O_WORKDIR

# Below here enter the commands to start your job

# Simple single process examples:

export CMD="blat -noHead $SUBJECT $QUERY $QUERY.psl"
$CMD

echo $CMD

exit;
