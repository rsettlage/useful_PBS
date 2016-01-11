#!/bin/bash


####ssh sfx063 

##### get all the R1 and R2 sequences into separate files
zcat *R1*fastq.gz >R1.fastq
zcat *R2*fastq.gz >R2.fastq

export MODULEPATH=/apps/packages/bio/modulefiles:$MODULEPATH
bio/velvet/1.2.08-bin