import pandas as pd
import glob

def parse_decimal(text):
    pieces = text.split(',')
    j1=['|'.join(pieces[i:i+2]) for i in xrange(0, len(pieces), 2)]
    return('.'.join(j1))

files = glob.glob('./*.csv')
for inputfile in files:
    outfile = inputfile.split('.csv')[0]+'_ok.txt'
    with open(inputfile, "r") as infile, open(outfile, "w") as ofile:
        for line in infile:
            ofile.write(parse_decimal(line))