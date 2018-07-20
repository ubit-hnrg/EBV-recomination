#!/bin/bash

trim=0.25
for f in "$@";
    do
    fname=$(basename $f); # remove path
    fbname=${fname%.*};   # remove extension
    
    inputpath=/data/EBV/msas/$fbname/;
    msafile=$inputpath$fbname.clustal; 
    trimmedout=$inputpath$fbname'_trimmed'$trim'.clustal'
    mapptrim=$inputpath'colnumb.txt'

    echo $trimmedout 
    #echo $msafile;

    trimal -in $msafile -out $trimmedout -gt $trim -colnumbering > $mapptrim
    done
#trimal -in ./fullrandseqs_test/fullrandseqs_test.clustal -out ./fullrandseqs_test/trimed.clustal -gt 0.25 -colnumbering > ./fullrandseqs_test/colnumb.txt
