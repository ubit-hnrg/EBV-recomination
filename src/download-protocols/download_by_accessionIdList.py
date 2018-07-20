from __future__ import print_function
from Bio import Entrez
import pandas as pd
import time
from os import listdir
from os.path import isfile, join
import argparse

parser = argparse.ArgumentParser(description='Download NCBI sequences using an accesion id file list')

parser.add_argument('-i','--inputfile',required=True)
parser.add_argument('-o','--outputpath',required=True)
args = parser.parse_args()

#args = parser.parse_args()
inputfile = args.inputfile
outputpath = args.outputpath

accession_ids = pd.read_csv(inputfile)


def dump_accessionList(accList,outputpath,address = "arieljberenstein@gmail.com"):
    Entrez.email =address
    for genome_id in accList:
        record = Entrez.efetch(db="nucleotide", id=genome_id, rettype="gb", retmode="text")
        fasta = Entrez.efetch(db="nucleotide", id=genome_id, rettype="fasta", retmode="text")


        filename = '{}genBankRecord_{}.gb'.format(outputpath,genome_id)
        fastafile = '{}{}.fasta'.format(outputpath,genome_id)

        print('Writing:{}'.format(genome_id))
        
        with open(filename, 'w') as f:
            f.write(record.read())
        with open(fastafile, 'w') as h:
            h.write(fasta.read())
    return None



def download_accIds_from_list_oneFile(list_of_accession_ids,outputfile = './output.gb',address = "arieljberenstein@gmail.com"):
    Entrez.email =address
    genomeAccessions = list_of_accession_ids

    search           = " ".join(genomeAccessions)
    handle           = Entrez.read(Entrez.esearch(db="nucleotide", term=search, retmode="xml"))
    genomeIds        = handle['IdList']
    records          = Entrez.efetch(db="nucleotide", id=genomeIds, rettype="gb")

    file_out = open(outputfile, "w")    # store each genomes .gb in separate files
    file_out.write(records.read())
    file_out.close()
    return None
    

def main():
    IDs = accession_ids.iloc[:,0].tolist()
    #outputfile = outputpath + inputfile + '.gb'
    #download_accIds_from_list(IDs,outputfile = outputfile)
    dump_accessionList(IDs,outputpath = outputpath)

if __name__ == "__main__":
    main()



