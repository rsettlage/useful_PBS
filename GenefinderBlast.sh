#!/bin/bash

## all PBS/torque specific commands are prepended with '#PBS'
####qsub /groups/DAC/useful_PBS/GenefinderBlast.sh -v Results=Ruhi_script.bln.txt,db_file=/groups/DAC/Morrill_Oct2013/CD46mg2inpBluescriptSK.fasta,Query=test.file
####find `pwd` -maxdepth 1 -iname '*R1*fastq.gz' -exec qsub /groups/DAC/useful_PBS/GenefinderBlast.sh -v Results=bln.txt,db_file=/groups/DAC/Morrill_Oct2013/CD46mg2inpBluescriptSK.fasta,Query='{}' \;

###########################################################################
## PBS stuff
####### job customization
#PBS -N Morrill_blst
## merge stdout and stderr
#PBS -j oe
## e-mail us when the job: '-M b' = begins, '-M a' = aborts, '-M e' = ends
#PBS -m a -M ruhi@vbi.vt.edu
## request queue sfx_q, sfxsmp with group list as sfx or sfxsmp, need to add a :smp to the nodes section
#PBS -W group_list=sfx
#PBS -q sfx_q
#PBS -lnodes=1:ppn=4
####PBS -l naccesspolicy=singlejob or #PBS -W x="NACCESPOLICY:SINGLEJOB"
######## specify resource allocation
## ask for 1 hour of wall time
#PBS -l walltime=20:00:00
######## additional sample resource reservations
## ask for 500MB of memory
###PBS -l mem=500mb
# end of PBS stuff
###########################################################################
## environment & variable setup
export MODULEPATH=/apps/packages/bio/modulefiles:$MODULEPATH
module load bio/blast/2.2.27
module load pigz/2.2.4    
module list
## end of loading modules
## end of environment & variable setup
###########################################################################
echo start:
date
echo jobid is $PBS_JOBID
JOBID=${PBS_JOBID%%.*}
echo jobnumber is $JOBID
SCRIPT=/groups/DAC/useful_PBS/GenefinderBlast.sh
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
query_FILE1=${Query##*/}
query_DIR=${Query%/*}

cd $work_DIR

cp -v $database_DIR/$database_FILE $local_DIR/$database_FILE
cp -v $query_DIR/$query_FILE1 $local_DIR/$query_FILE1

	echo is paired, so get file2
	READ1_indicator="_R1_"
	READ2_indicator="_R2_"
	query_FILE2=${query_FILE1/"$READ1_indicator"/"$READ2_indicator"}
	echo cp $query_DIR/$query_FILE2 $local_DIR/$query_FILE2
	cp -v $query_DIR/$query_FILE2 $local_DIR/$query_FILE2

result_suffix=.$Results  
result_FILE1=$query_FILE1$result_suffix
echo
echo query is $query_FILE1

result_FILE2=$query_FILE2$result_suffix
echo database is $database_FILE
echo query is $query_FILE2
echo database directory is $database_DIR
echo result directory is $query_DIR
echo result will go in $result_FILE
echo

cd $local_DIR
echo directory we are working in is...
pwd
echo

query_EXT=${query_FILE1##*.}      
if [ ${query_EXT} == 'gz' ];then
	echo unzipping $query_FILE1
	gzip -d $query_FILE1
	new_query_FILE1=${query_FILE1%.*}
	echo new query file is $new_query_FILE1
	query_FILE1=$new_query_FILE1
fi
query_EXT=${query_FILE1##*.}
if [ ${query_EXT} == 'fastq' ] || [ ${query_EXT} == 'fq' ];then
	fasta_suffix=.fasta
	reads_FILE1=$query_FILE1
	new_query_FILE1=$query_FILE1$fasta_suffix
	cat $query_FILE1 | perl -e '$i=0;while(<>){if(/^\@/&&$i==0){s/^\@/\>/;print;}elsif($i==1){s/\./N/g;print;$i=-3}$i++;}' > $new_query_FILE1
	query_FILE1=$new_query_FILE1
fi

query_EXT=${query_FILE2##*.}      
if [ ${query_EXT} == 'gz' ];then
	echo unzipping $query_FILE2
	gzip -d $query_FILE2
	new_query_FILE2=${query_FILE2%.*}
	echo new query file is $new_query_FILE2
	query_FILE2=$new_query_FILE2
fi
query_EXT=${query_FILE2##*.}
if [ ${query_EXT} == 'fastq' ] || [ ${query_EXT} == 'fq' ];then
	fasta_suffix=.fasta
	reads_FILE2=$query_FILE2
	new_query_FILE2=$query_FILE2$fasta_suffix
	cat $query_FILE2 | perl -e '$i=0;while(<>){if(/^\@/&&$i==0){s/^\@/\>/;print;}elsif($i==1){s/\./N/g;print;$i=-3}$i++;}' > $new_query_FILE2
	query_FILE2=$new_query_FILE2
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
	ls -lah
	exit
fi


echo
echo
pwd
echo directory listing:
ls -la
echo ready to run blast
### blastn is nuc vs nuc, blastx is nuc vs prot  ######## <----------------------------------------set this---------####
###EXEC_CMD2="blastx -outfmt 6 -num_threads 12 -evalue 1e-6 -num_descriptions 2 -db $database_FILE -query $query_FILE -out $result_FILE"
###EXEC_CMD2=`blastn -outfmt '6 qseqid sseqid pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore' -num_threads 12 -max_target_seqs 5 -evalue 1e-4 -db $database_FILE -query $query_FILE -out $result_FILE`
###EXEC_CMD2="tblastx -outfmt 6 -num_threads 12 -evalue 1e-6 -db $database_FILE -query $query_FILE -out $result_FILE"
###-word_size 19 for closer matches

##`blastn -outfmt 6 -num_threads 8 -evalue 1e-4 -db $database_FILE -query $query_FILE -out $result_FILE`
`blastn -outfmt '6 qseqid sseqid pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore' -num_threads 3 -evalue 1e-4 -db $database_FILE -query $query_FILE1 -out $result_FILE1`
`blastn -outfmt '6 qseqid sseqid pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore' -num_threads 3 -evalue 1e-4 -db $database_FILE -query $query_FILE2 -out $result_FILE2`


####vi zzruhi_bash.sh
####i

awk '!x[$1]++' $result_FILE1 > $result_FILE1.nodup
awk '!x[$1]++' $result_FILE2 > $result_FILE2.nodup
cat $result_FILE1.nodup $result_FILE2.nodup > R1R2.dups
awk '{print $1}' R1R2.dups | sort -k 1,1 | uniq -u > R1R2.uniq
awk '{print $1}' R1R2.dups | sort -k 1,1 | uniq -d > R1R2.notuniq
`sed -i 's/^/@/' R1R2.notuniq`
`sed -i 's/^/@/' R1R2.uniq`

ls -lah


#####hit esc
#####:wq

echo $reads_FILE1

awk -f /groups/DAC/useful_PBS/GenefinderSave_reads_file1_from_file2.awk R1R2.notuniq $reads_FILE1 > $result_FILE1.paired_plasmid_hits.fastq
awk -f /groups/DAC/useful_PBS/GenefinderSave_reads_file1_from_file2.awk R1R2.notuniq $reads_FILE2 > $result_FILE2.paired_plasmid_hits.fastq
awk -f /groups/DAC/useful_PBS/GenefinderSave_reads_file1_from_file2.awk R1R2.uniq $reads_FILE1 > $result_FILE1.one_end_plasmid_hits.fastq
awk -f /groups/DAC/useful_PBS/GenefinderSave_reads_file1_from_file2.awk R1R2.uniq $reads_FILE2 > $result_FILE2.one_end_plasmid_hits.fastq

 
echo files on node:
ls -lah


cp -v *bln.txt $query_DIR/
cp -v *plasmid_hits.fastq $query_DIR/
####cp -v R1R2* $query_DIR/


## end execution stage
###########################################################################
echo end:
date

rm -r $local_DIR/

echo $query_FILE
echo finished
exit
