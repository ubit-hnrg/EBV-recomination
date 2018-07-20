#!/bin/bash 
#this file is completely harcoded and provisory
inputfile='../reference_genomes/bwa/NC_007605.fasta'
queryfile='../msas/fullrandseqs_test/fullrandseqs_test.fa'
outputfile='./outputfile.fa'

# bwa does a local local aligment
bwa index $inputfile 
bwa mem $inputfile $queryfile > $outputfile