import argparse
import pandas as pd
import os 
parser = argparse.ArgumentParser(description='cat fasta files given a ')

parser.add_argument('-i','--idfile',required=True)
parser.add_argument('-p','--path',required=True)
parser.add_argument('-o','--outputfile',required=True)


if(True):
    args = parser.parse_args()
    idfile = args.idfile
    path = args.path
    outputfile = args.outputfile
else:
    idfile = './data/EBV/reference_genomes/refIds.txt'
    path = './data/EBV/byACCIDs/'
    outputfile = './data/EBV/multiple_fastas/references.fa'
    
# use id list
ids = pd.read_csv(idfile,sep ='\t').iloc[:,0]

files = path + ids + '.fasta'

okfiles = [f for f in files if os.path.isfile(f)]

catfiles = ' '.join(okfiles)

os.system('cat ' + catfiles + '> '+ outputfile)