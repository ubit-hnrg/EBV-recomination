from Bio import Entrez
import pandas as pd
import time
from os import listdir
from os.path import isfile, join
import argparse


parser = argparse.ArgumentParser(description='Download NCBI sequences using an explicit search_term')

parser.add_argument('-s','--search_term',required=True)
parser.add_argument('-o','--outputfile',required=True)
parser.add_argument('-r','--retmax',default = 200000)

args = parser.parse_args()

#args = parser.parse_args()
search_term = args.search_term
outputfile = args.outputfile
RetMax = args.retmax

# example
#search_term = "gastric[All Fields] AND ((viruses[filter] OR 'Human gammaherpesvirus 4'[Organism]) AND ('150000'[SLEN] : '180000'[SLEN])"


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

        print 'Writing:{}'.format(genome_id)
        
        with open(filename, 'w') as f:
            f.write(record.read())
        with open(fastafile, 'w') as h:
            h.write(fasta.read())
    return None

def get_accessions_from_searchTerm(search_term,address = "arieljberenstein@gmail.com"):
    Entrez.email =address
    handle = Entrez.esearch(db='nucleotide', term=search_term,RetMax = RetMax)
    genome_ids = Entrez.read(handle)['IdList']
        
    #file_search_mapp =  outputpath + 'idList'+  search_term + '.txt'
    
    if len(genome_ids)!=0:
        accList = get_accessions_of_idList(genome_ids)
    else:
        accList =[]
    return(accList)

def write_accList(accList,search_term,fname):
    with open(fname, 'w') as g:
        g.write('search_term : {} \n'.format(search_term))
        g.write("\n".join(accList))
        g.write("\n")
    return(None)


def main():
    print search_term
    AccIds = get_accessions_from_searchTerm(search_term)
    if len(AccIds)!=0:
        write_accList(accList=AccIds,search_term=search_term,fname=outputfile)
    else: 
        print 'WARNING: none search term found. The output file was not written. Please, review your query'

    return(None)
    #dump_accessionList(accList,outputpath=outputpath)
    
    


if __name__ == "__main__":
    main()
