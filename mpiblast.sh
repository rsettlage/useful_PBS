#!/bin/bash
## sample PBS/torque script for software which uses MPI

## all PBS/torque specific commands are prepended with '#PBS'
## please remember that those are interpreted only by PBS/torque!

###########################################################################
## environment & variable setup

####### job customization
## name our process
#PBS -N test_MPI_blastt
## merge stdout and stderr
#PBS -j oe
# write PBS output files to output directory
#PBS -o results/
## e-mail us when the job: '-M b' = begins, '-M a' = aborts, '-M e' = ends
# #PBS -m abe -M rsettlage@vbi.vt.edu
####### end of job customization


######## specify resource allocation

## ask for 1 hour of wall time
## if your job takes longer than that, it will be automatically killed
#PBS -l walltime=90:00:00

## ask for 2 nodes and 12 processor slots per node
## WARNING! if your job is not MPI aware, set it to: '-lnodes=1:ppn=1'
#PBS -lnodes=4:ppn=48

# request specific queue. in most cases you want to leave it as 'default'
#PBS -q default

######## additional sample resource reservations
## ask for 500MB of memory
## if your job takes more memory than that, it will be automatically killed
# #PBS -l mem=500mb

######## end of resource allocation

# end of environment & variable setup
###########################################################################

## actual command and its options to be executed by this job
#TARGET_DIR=data/
#TARGET_DB=month.aa
#QUERY=queries/queries.faa

QUERY=queries/zinc-finger.faa
TARGET_DIR=/common/databases/fasta
TARGET_DB=est_human
TYPE=blastn
RESULT_DIR=results/
EXTRA_OPTS="-e .001 -m 8 -v 1 -b 10 -M BLOSUM62 --removedb"

###########################################################################
## environment & variable setup
## likely there are no changes required in this section


## begin loading bright cluster manager modules, which are required on shadowfax
. /etc/profile.d/modules.sh

module add shared
module add hpl
module add openmpi
module add gotoblas/core/64
## end of loading modules


## check how many processor slots we have allocated
## and use that value later for 'mpiexec -np ${NUM_SLOTS}'
NUM_SLOTS=$(/bin/sed "s/  //g" < $PBS_NODEFILE | wc -l)

## end of environment & variable setup
###########################################################################


###########################################################################
## begin execution stage

## change to the working directory
cd ${PBS_O_WORKDIR}

# make sure we have all directories created
! [ -d ${RESULT_DIR} ] && mkdir ${RESULT_DIR}

## record what time we began this job
echo -e "\n--> begin ${PBS_JOBNAME} with id ${PBS_JOBID} at: " $(date) "\n"
echo -e "\n--> working in ${PBS_O_WORKDIR} with ${NUM_SLOTS} processor slots \n"

# staging of data for mpiblast

echo -e "\n--> peforming mpiformatdb \n"
DESTINATION=`egrep -i ^shared ~/.ncbirc | cut -f 2 -d '='`
mpiformatdb --nfrags=${NUM_SLOTS} -i ${TARGET_DIR}/${TARGET_DB} -n ${DESTINATION}

echo -e "\n--> timestamp: $(date) \n"

## execute our command
echo -e "\n--> peforming mpiblast \n"

time mpiexec -np ${NUM_SLOTS} -machinefile ${PBS_NODEFILE} mpiblast -d ${TARGET_DB} -i ${QUERY} -p ${TYPE} -o ${RESULT_DIR}/${PBS_JOBNAME}-${PBS_JOBID}-results.txt ${EXTRA_OPTS}

# clean up local databases
echo -e "\n--> peforming mpiblast_cleanup \n"
time mpiexec -np ${NUM_SLOTS} -machinefile ${PBS_NODEFILE} mpiblast_cleanup

## record what time we end this job
echo -e "\n--> end ${PBS_JOBID} at: " $(date) "\n"

## end execution stage
###########################################################################
