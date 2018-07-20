from __future__ import print_function
from Bio import Entrez
import pandas as pd
import time
from os import listdir
from os.path import isfile, join
import argparse


parser = argparse.ArgumentParser(description='Download NCBI sequences using an explicit search_term')

parser.add_argument('-s','--search_term',required=True)
parser.add_argument('-o','--outputpath',required=True)
parser.add_argument('-r','--retmax',default = 2)

args = parser.parse_args()

#args = parser.parse_args()
search_term = args.search_term
outputpath = args.outputpath
RetMax = args.retmax

def get_accessions_of_idList(idList):
    handle = Entrez.efetch(db='nucleotide', id= ','.join(idList),retmode="xml")
    result = Entrez.read(handle)
    accessions = []
    for s in result:
        accessions.append(s['GBSeq_primary-accession'])
    return accessions


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


def main():
    handle = Entrez.esearch(db='nucleotide', term=search_term,RetMax = RetMax)
    genome_ids = Entrez.read(handle)['IdList']
    file_search_mapp =  outputpath + 'idList'+  search_term + '.txt'
    accList = get_accessions_of_idList(genome_ids)
    
    with open(file_search_mapp, 'w') as g:
        g.write('search_term : {} \n'.format(search_term))
        g.write("\n".join(accList))
        g.write("\n")

    dump_accessionList(accList,outputpath=outputpath)
        
if __name__ == "__main__":
    main()
