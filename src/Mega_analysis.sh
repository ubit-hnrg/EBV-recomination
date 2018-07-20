#!/bin/bash

for i in "$@"
do
case $i in
    -m=*|--msa=*)
    msafile="${i#*=}"

    ;;
    -o=*|--out_dir=*)
    outdir="${i#*=}"
    ;;
    
    -O=*|--out_file=*)
    outfile="${i#*=}"
    ;;

    *)
    
    ;;
esac
    
done

rateBySiteFile='/home/ariel/EBV/MEGA/site_by_site_rates_nucleotide.mao'

echo runing rate by site 
megacc -a $rateBySiteFile -d $msafile -o $outdir$outfile


