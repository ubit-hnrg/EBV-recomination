#!/bin/bash

cat accession_ids_EBV_Correia_2017.txt ncbi_chiara2016_acc_ids.txt |cut -f1 -d . |sort |uniq |sort > correia-and-chaira-unique.txt
