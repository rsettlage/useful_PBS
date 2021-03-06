#!/bin/bash

## all PBS/torque specific commands are prepended with '#PBS'
####qsub /groups/DAC/useful_PBS/multi_PBS_blastplus_blastP_FILE.sh -v Results=blP.txt,db_file=/groups/DAC/vector_contaminants_file/adaptors_contaminants.fa,Query=test.fa
####find `pwd` -maxdepth 1 -iname '*fasta' -exec qsub /groups/DAC/useful_PBS/multi_PBS_blastplus_blastP_FILE.sh -v Results=blstP.txt,db_file=,Query='{}' \;

###########################################################################
## PBS stuff
####### job customization
#PBS -N CM_Blstp
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
#PBS -l walltime=50:00:00
######## additional sample resource reservations
## ask for 500MB of memory
###PBS -l mem=500mb
# end of PBS stuff
###########################################################################
## environment & variable setup
export MODULEPATH=/apps/packages/bio/modulefiles:$MODULEPATH
module load bio/blast/2.2.27
module load bio/blat
module list
## end of loading modules
## end of environment & variable setup
###########################################################################
echo start:
date
echo jobid is $PBS_JOBID
JOBID=${PBS_JOBID%%.*}
echo jobnumber is $JOBID
SCRIPT=/groups/DAC/useful_PBS/multi_PBS_blastplus_blastP_FILE.sh
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
query_FILE=${Query##*/}
query_DIR=${Query%/*}

result_suffix=.$Results  
result_FILE=$query_FILE$result_suffix
echo database is $database_FILE
echo query is $query_FILE
echo database directory is $database_DIR
echo result directory is $query_DIR
echo result will go in $result_FILE
echo

cd $work_DIR

cp -v $database_DIR/$database_FILE $local_DIR/$database_FILE
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
query_EXT=${query_FILE##*.}
if [ ${query_EXT} == 'fastq' ] || [ ${query_EXT} == 'fq' ];then
	fasta_suffix=.fasta
	new_query_FILE=$query_FILE$fasta_suffix
	cat $query_FILE | perl -e '$i=0;while(<>){if(/^\@/&&$i==0){s/^\@/\>/;print;}elsif($i==1){s/\./N/g;print;$i=-3}$i++;}' > $new_query_FILE
	query_FILE=$new_query_FILE
fi

ls -lah

if [ -e $database_FILE ];then
	####first look to see if it is zipped, if so, unzip it and set the extension approprately
	db_EXT=${database_FILE##*.}
	echo db is currently a $db_EXT file
	if [ ${db_EXT} == 'bz2' ];then
		echo unzipping $database_FILE
		bunzip2 -v $database_FILE
		is_FASTA=${database_FILE%.*}
		database_FILE=$is_FASTA
		echo new database is $database_FILE
	fi
	if [ ${db_EXT} == 'gz' ];then
		echo unzipping $database_FILE
		gzip -d $database_FILE
		is_FASTA=${database_FILE%.*}
		database_FILE=$is_FASTA
		echo new database is $database_FILE
	fi
	echo databasefile $database_FILE is ready, so index it   
	EXEC_CMD1="makeblastdb -in $database_FILE -dbtype nucl"  ###prot or nucl  <---------------set this----------####
	echo $EXEC_CMD1  
	${EXEC_CMD1}
else
	echo no file to index
	pwd
	ls -la
	exit
fi

echo
echo
pwd
echo directory listing:
ls -la
echo ready to run blast
echo tblastn for protein query using nuleotide database
echo tblastx for translated nucleotide database and query
echo blastp for protein database and query
echo blastx for protein database and nucleotide query
###EXEC_CMD2="blastx -outfmt 6 -num_threads 12 -evalue 1e-6 -num_descriptions 2 -db $database_FILE -query $query_FILE -out $result_FILE"
`tblastn -outfmt '6 qseqid sseqid pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore' -num_threads 12 -max_target_seqs 10 -word_size 3 -evalue 1e-4 -db $database_FILE -query $query_FILE -out $result_FILE`
###-word_size 19 for closer matches

##`blastp  -outfmt 6 -num_threads 12 -max_target_seqs 1 -query $query_FILE -db $database_FILE -out $result_FILE`

## for xml output
##`blastx -outfmt 5 -num_threads 4 -evalue 1e-4 -max_target_seqs 5 -db $database_FILE -query $query_FILE -out $result_FILE`

##for trinity, use this format
##`blastp -outfmt 6 -max_target_seqs 1 -num_threads 12 -db $database_FILE -query $query_FILE -out $result_FILE`

##for blast2go, use blastall and output xml table
#result_suffix=XML 
#result_FILE=$result_FILE.$result_suffix
#`/apps/packages/bio/blast/2.2.22/bin/blastall -a 12 -b 20 -v 20 -p blastx -e 1e-6 -m 0 -d $database_FILE -i $query_FILE -o $result_FILE`


##`awk '{ if(!x[$1]++) {print $0; bitscore=($14-6)} else { if($14>bitscore) print $0} }' $result_FILE >$result_FILE.trimmed`

echo 
echo files on node:
ls -lah
echo copying result over using cp -v $result_FILE $query_DIR/$result_FILE
cp -v $result_FILE $query_DIR/$result_FILE
##cp -v $result_FILE.trimmed $query_DIR/$result_FILE

## end execution stage
###########################################################################
echo end:
date

rm -r $local_DIR/

echo $query_FILE
echo finished
exit
