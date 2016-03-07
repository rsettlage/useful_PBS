#!/bin/bash

## all PBS/torque specific commands are prepended with '#PBS'
####
####find `pwd` -type f \( -iname '*.fa' \) -exec qsub /groups/DAC/useful_PBS/multi_PBS_qiime_FILE.sh -v db_file=/groups/DAC/16s/gg_12_10_otus/rep_set/97_otus.fasta,db_tax_file=/groups/DAC/16s/gg_12_10_otus/taxonomy/97_otu_taxonomy.txt,Query='{}' \;

###########################################################################
## PBS stuff
####### job customization
#PBS -N q_t
## merge stdout and stderr
#PBS -j oe
## e-mail us when the job: '-M b' = begins, '-M a' = aborts, '-M e' = ends
#PBS -m a -M rsettlage@vbi.vt.edu
## request queue sfx_q, sfxsmp with group list as sfx or sfxsmp, need to add a :smp to the nodes section
#PBS -W group_list=sfx
#PBS -q sandybridge_q
#PBS -lnodes=1:ppn=12
####PBS -l naccesspolicy=singlejob or #PBS -W x="NACCESPOLICY:SINGLEJOB"
######## specify resource allocation
## ask for 1 hour of wall time
#PBS -l walltime=200:00:00
######## additional sample resource reservations
## ask for 500MB of memory
###PBS -l mem=500mb
# end of PBS stuff
###########################################################################
## environment & variable setup
export MODULEPATH=/apps/packages/bio/modulefiles:$MODULEPATH
module load bio/qiime
module list
## end of loading modules
## end of environment & variable setup
###########################################################################
echo start:
date
echo jobid is $PBS_JOBID
JOBID=${PBS_JOBID%%.*}
echo jobnumber is $JOBID
SCRIPT=/groups/DAC/useful_PBS/multi_PBS_qiime_FILE.sh
more $SCRIPT
module list
###########################################################################
## begin execution stage

work_DIR=$PBS_O_WORKDIR
local_DIR=/localscratch/$JOBID
mkdir $local_DIR

echo originating directory is $work_DIR
echo node directory is $local_DIR
echo
echo

database_FILE=${db_file##*/}
database_DIR=${db_file%/*}
database_tax_FILE=${db_tax_file##*/}
database_tax_DIR=${db_tax_file%/*}
query_FILE=${Query##*/}
query_FILE_base=${query_FILE%.*}
query_DIR=${Query%/*}
result_DIR=$query_DIR/$query_FILE._qiime

result_suffix=_$Results  
result_FILE=$query_FILE$result_suffix
echo database is $database_FILE
echo query is $query_FILE
echo database directory is $database_DIR
echo result directory is $query_DIR
echo result will go in $result_DIR
echo

cd $work_DIR
mkdir $result_DIR

cp -v $database_DIR/$database_FILE $local_DIR/$database_FILE
cp -v $database_tax_DIR/$database_tax_FILE $local_DIR/$database_tax_FILE
cp -v $query_DIR/$query_FILE $local_DIR/$query_FILE
cd $local_DIR
echo directory we are working in is...
pwd
echo

query_EXT=${query_FILE##*.}
if [ ${query_EXT} == 'gz' ];then
	echo unzipping $query_FILE
	gzip -d $query_FILE
	new_query_FILE=${query_FILE%.*}
	echo new query file is $new_query_FILE
	query_FILE=$new_query_FILE
fi

ls -lah

echo
echo
pwd
echo local node directory listing:
ls -lahR
echo
echo ready to run qiime


echo "pick_otus.py -i $query_FILE -m usearch --word_length 64 --db_filepath $database_FILE -o usearch_results/ --minsize 2"
`pick_otus.py -i $query_FILE -m usearch --word_length 64 --db_filepath $database_FILE -o usearch_results/ --minsize 2`
##ls -lahR
echo
pick_otu_results_suffix="_otus.txt"
pick_otu_results=$query_FILE_base$pick_otu_results_suffix
pick_rep_set_out_suffix="_representative_set.fasta"
pick_rep_set_out=$query_FILE_base$pick_rep_set_out_suffix
echo "pick_rep_set.py -i usearch_results/$pick_otu_results -f $query_FILE -o $pick_rep_set_out"
`pick_rep_set.py -i usearch_results/$pick_otu_results -f $query_FILE -o $pick_rep_set_out`
##ls -lahR
echo
assign_taxonomy_out_suffix="_RDP_assigned_taxonomy"
assign_taxonomy_out=$query_FILE_base$assign_taxonomy_out_suffix
echo "assign_taxonomy.py -m rdp -i $pick_rep_set_out -r $database_FILE -t $database_tax_FILE --rdp_max_memory=3000 -o $assign_taxonomy_out"
`assign_taxonomy.py -m rdp -i $pick_rep_set_out -r $database_FILE -t $database_tax_FILE --rdp_max_memory=3000 -o $assign_taxonomy_out`
##ls -lahR
echo
assign_taxonomy_out_FILE_suffix="_representative_set_tax_assignments.txt"
assign_taxonomy_out_FILE=$assign_taxonomy_out/$query_FILE_base$assign_taxonomy_out_FILE_suffix
make_out_table_out_suffix="_otu_table.biom"
make_out_table_out=$assign_taxonomy_out/$query_FILE_base$make_out_table_out_suffix
echo "make_otu_table.py -i usearch_results/$pick_otu_results -t $assign_taxonomy_out_FILE -o $make_out_table_out"
`make_otu_table.py -i usearch_results/$pick_otu_results -t $assign_taxonomy_out_FILE -o $make_out_table_out`
##ls -lahR
echo
convert_biom_out_suffix="_otu_table.txt"
convert_biom_out=$assign_taxonomy_out/$query_FILE_base$convert_biom_out_suffix
echo "convert_biom.py -i $make_out_table_out -o $convert_biom_out -b"
`convert_biom.py -i $make_out_table_out -o $convert_biom_out -b`
##ls -lahR


echo 
echo files on node:
ls -lahR
echo copying result over using cp -v $local_DIR $result_DIR
##need to remove the initial files before copying over to save some space
rm $local_DIR/$database_FILE
rm $local_DIR/$query_FILE
rm $local_DIR/$database_tax_FILE
ls -lahR
cp -v -R ./ $result_DIR/


## end execution stage
###########################################################################
echo end:
date

rm -r $local_DIR/

echo $query_FILE
echo finished
exit
