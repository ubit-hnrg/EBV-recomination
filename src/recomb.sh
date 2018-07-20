#!/bin/bash


gen='17494219'
mkdir $gen
group='N.America'
fastas='/data/EBV/by_gene/'$group/$gen'_msa.fa';
src_path='/home/ariel/repo/TBI/EBV-db/src/'
output_dir='./'

path_lk_file='/home/ariel/Projects/Gutierrez/EBV-recomb/recomb/lk_files/lk_n50_t0.1'

for msa in $fastas;
do
    python $src_path'get_mod_fasta_format.py' --msafile=$msa
    auxpath=${msa%.fa};
    modFasta=$auxpath'.modFasta'
    #base=$(basename $modFasta);
    #group=${base%_msa*};
    comm='-seq '$modFasta' -prefix ./'$gen'/'
    echo $comm
    convertLDHat $comm ;
done

nseqs=16;
lkpref='n'$nseqs'_';
#lkfile=$lkpref'new_lk.txt';
# notar que el numero de secuencias de esta tabla debe ser el doble (por ser un genoma diploide)
#lkgen -lk path_lk_file -nseq $nseqs -prefix $lkpref

lkfile='/home/ariel/Projects/Gutierrez/EBV-recomb/recomb/lk_files/n16_new_lk.txt'

pair='-seq ./'$gen'/sites.txt -loc ./'$gen'/locs.txt -lk '$lkfile' -prefix ./'$gen'/'
echo "runing '\n'"
echo $pair
pairwise $pair

# perfecto falta optimizar parametros que si no son inicializados desde el principio el prorgrama se tara y pregunta a mitad de camino




python ~/repo/TBI/EBV-db/src/get_mod_fasta_format.py --msafile=/data/EBV/msas/up150k/all150.fa_1_5000.fa
mv /data/EBV/msas/up150k/all150.modFasta /data/EBV/msas/up150k/all150_1_5000.modFasta
convertLDHat -seq /data/EBV/msas/up150k/all150_1_5000.modFasta -prefix 'allEBV_'


python ~/repo/TBI/EBV-db/src/get_mod_fasta_format.py --msafile=/data/EBV/msas/up150k/30_161_pos_0-5k.fa
convertLDHat -seq /data/EBV/msas/up150k/30_161_pos_0-5k.modFasta

lkgen -lk /home/ariel/Projects/Gutierrez/EBV-recomb/recomb/lk_files/lk_n192_t0.001 -nseq 262 -prefix n262

lkgen -lk /home/ariel/Projects/Gutierrez/EBV-recomb/recomb/lk_files/n -nseq 262 -prefix n262
interval -seq ./sites.txt -loc locs.txt -lk 
