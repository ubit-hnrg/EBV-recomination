import argparse
from Bio import Alphabet
from Bio import AlignIO

parser = argparse.ArgumentParser(description='convert msa file format')
parser.add_argument('-i','--inputfile',required=True)
parser.add_argument('-f','--inputformat',required=True)
parser.add_argument('-o','--outputfile',required=True)  
parser.add_argument('-F','--outputformat',required=True)


# parse args
args = parser.parse_args()
inputfile = args.inputfile
inputformat = args.inputformat
outputfile = args.outputfile
outputformat = args.outputformat



AlignIO.convert(inputfile,inputformat,outputfile,outputformat,alphabet=Alphabet.generic_dna)

# example how to run
#for msa in */*msa*.fa;do echo $msa;output=${msa%.*}.nex;echo $output;python ~/repo/TBI/EBV-db/src/align_format_convert.py -i $msa -f fasta -o $output -F nexus;done