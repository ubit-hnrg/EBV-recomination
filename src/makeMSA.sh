#!/bin/bash

trim=0.05
ref_seq='/data/EBV/byACCIDs/NC_007605.fasta'    

for i in "$@"
do
case $i in
    -f=*|--file=*)
    file="${i#*=}"

    ;;
    -g=*|--gap=*)
    gap="${i#*=}"
    
    ;;
    
    *)
            # unknown option
    ;;
esac
    
done

echo $gap
echo $file

seqpath='/data/EBV/byACCIDs/'
fname=$(basename $file) # remove path
fbname=${fname%.*}   # remove extension

echo $fname
echo $fbname

# create output folder 
mkdir /data/EBV/msas/$fbname
outpath=/data/EBV/msas/$fbname/

## harcoded paths to outputfiles
echo merging fasta files for group $fbname;
input=$file
auxfa=$outpath'aux_'$fbname.fa;
#auxfa1=$outpath'aux1_'$fbname.fa;
outfa=$outpath$fbname.fa;
outmsa=$outpath$fbname'_msa_gap'$gap.fa;
echo $
echo $outpath
echo $auxfa
echo $outfa
echo $outmsa

sudo python /home/ariel/repo/TBI/EBV-db/src/cat-fasta-files.py --idfile=$input --path=$seqpath --outputfile=$auxfa;
cat $auxfa |tr 'N' '-' > $outfa
sed -i -e 's/\-C\_00/NC_00/g' $outfa

sudo rm $auxfa 
#sudo rm $auxfa1

echo runing msa-mafftfor group $fbname;
#sudo mafft --clustalout --thread 4 $outfa > $outclust;

# this option is accurate but its take a lot of memmory and time. For only 6 full EBV sequences it take 30Gb of memory at the very begining! 
#sudo mafft --localpair --lop -10 --lexp -0.1  $outfa > $outmsa

## here we are using the default seting --6merpair, where distance is computed based on 6mers. Nevertheless we are changing the open penalty score
sudo mafft --op $gap  $outfa > $outmsa

#echo runing rate4site
#sudo rate4site -s $outclust -o $outpath$fbname_rate4site.res

echo 'TERMINO EL ALINEAMIENTO EXITOSAMENTE'

ffname=$(basename $outmsa); # remove path
ffbname=${ffname%.*};   # remove extension
echo  $ffname
echo $ffbname
    
echo $inputpath
trimmedout=$outpath$ffbname'_trimmed'$trim'.fas'
echo $trimmedout
mapptrim=$outpath'colnumb.txt'
echo $mapptrim

## trim msa file
echo $trim
trimal -in $outmsa -out $trimmedout -gt $trim -colnumbering > $mapptrim

## run mappings 
python /home/ariel/repo/TBI/EBV-db/src/trimmed_2_ref.py -m $outmsa -r $ref_seq -M $mapptrim -o $outmsa'_whole_mapping.file'










