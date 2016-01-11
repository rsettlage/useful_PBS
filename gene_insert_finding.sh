#!/bin/bash

## all PBS/torque specific commands are prepended with '#PBS'
####qsub /groups/DAC/useful_PBS/multi_PBS_blastplus_FILE.sh -v Results=bln.txt,db_file=/groups/DAC/vector_contaminants_file/adaptors_contaminants.fa,Query=test.fa
####find `pwd` -iname '*gz' -exec qsub /groups/DAC/useful_PBS/multi_PBS_blastplus_FILE.sh -v Results=bln.txt,db_file=,Query='{}' \;

###########################################################################
## PBS stuff
####### job customization
#PBS -N Peanut_blst
## merge stdout and stderr
#PBS -j oe
## e-mail us when the job: '-M b' = begins, '-M a' = aborts, '-M e' = ends
#PBS -m a -M rsettlage@vbi.vt.edu
## request queue sfx_q, sfxsmp with group list as sfx or sfxsmp, need to add a :smp to the nodes section
#PBS -W group_list=sfx
#PBS -q sandybridge_q
#PBS -lnodes=1:ppn=6
####PBS -l naccesspolicy=singlejob or #PBS -W x="NACCESPOLICY:SINGLEJOB"
######## specify resource allocation
## ask for 1 hour of wall time
#PBS -l walltime=80:00:00
######## additional sample resource reservations
## ask for 500MB of memory
###PBS -l mem=500mb
# end of PBS stuff

align_FILES=$FILE_NAME1
if [ $PAIRED = "Y" ]; then
	echo is paired, so get file2
	READ1_indicator="_R1_"
	READ2_indicator="_R2_"
	FILE_NAME2=${FILE_NAME1/"$READ1_indicator"/"$READ2_indicator"}
	echo cp $FILE_DIR1/$FILE_NAME2 $node_DIR
	cp -v $FILE_DIR1/$FILE_NAME2 $node_DIR
	align_FILES=$FILE_NAME1" "$FILE_NAME2
fi

####vi zzruhi_bash.sh
####i

pigz -d *paired_matched.fastq.gz
cat *R1*bln.txt > R1_hits.bln.txt
cat *R2*bln.txt > R2_hits.bln.txt
awk '!x[$1]++' R1_hits.bln.txt > R1_hits.bln.txt_nodup
awk '!x[$1]++' R2_hits.bln.txt > R2_hits.bln.txt_nodup
cat R1_hits.bln.txt_nodup R2_hits.bln.txt_nodup > R1R2_hits.bln.txt_dups
awk '{print $1}' R1R2_hits.bln.txt_dups |sort -k 1,1|uniq -u > R1R2_hits.bln.txt_uniq
awk '{print $1}' R1R2_hits.bln.txt_dups |sort -k 1,1|uniq -d > R1R2_hits.bln.txt_notuniq

#####hit esc
#####:wq

cat *R1*bln.txt > R1_hits.bln.txt
cat *R2*bln.txt > R2_hits.bln.txt
awk '!x[$1]++' R1_hits.bln.txt > R1_hits.bln.txt_nodup
awk '!x[$1]++' R2_hits.bln.txt > R2_hits.bln.txt_nodup
cat R1_hits.bln.txt_nodup R2_hits.bln.txt_nodup > R1R2_hits.bln.txt_dups
awk '{print $1}' R1R2_hits.bln.txt_dups |sort -k 1,1|uniq -u > R1R2_hits.bln.txt_uniq
awk '{print $1}' R1R2_hits.bln.txt_dups |sort -k 1,1|uniq -d > R1R2_hits.bln.txt_notuniq
sed -i 's/^/@/' ../R1R2_hits.bln.txt_notuniq
sed -i 's/^/@/' ../R1R2_hits.bln.txt_uniq

awk -f /groups/DAC/useful_programs/remove_reads_file1_from_file2.awk ../R1R2_hits.bln.txt_notuniq all_R1_lane2_AT_QT_paired_matched.fastq >R1_paired_plasmid_hits.fastq
awk -f /groups/DAC/useful_programs/remove_reads_file1_from_file2.awk ../R1R2_hits.bln.txt_notuniq all_R2_lane2_AT_QT_paired_matched.fastq >R2_paired_plasmid_hits.fastq
awk -f /groups/DAC/useful_programs/remove_reads_file1_from_file2.awk ../R1R2_hits.bln.txt_uniq all_R1_lane2_AT_QT_paired_matched.fastq >R1_one_end_plasmid_hits.fastq 
awk -f /groups/DAC/useful_programs/remove_reads_file1_from_file2.awk ../R1R2_hits.bln.txt_uniq all_R2_lane2_AT_QT_paired_matched.fastq >R2_one_end_plasmid_hits.fastq