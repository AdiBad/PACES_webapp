#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###################################################################
##           Generating files necessary to use STRING            ##
###################################################################
##  - fetch FASTA files for proteins                             ## 
###################################################################

""" 
Requires: filteredData.tsv (generated with filter.py)
Output: acetylation.tsv
"""

import pandas as pd
import requests

data = pd.read_csv('../Output/filteredData.tsv',sep='\t')
uniqueID = data['Protein'].unique()

"""
Function to create a multifasta file, containing the sequence information in fasta format accessed via UniProt
for each protein that is specified in an input array.

IMPORTANT: when the ID can not be found on UniProt, the protein will be missing from the multifasta file.
"""
def getFasta(array):
    try:
        all_fasta = open("../Output/filteredDataSeq.fasta", "w")
        for x in array:
            response = requests.get('https://www.uniprot.org/uniprot/' + x + '.fasta')
            all_fasta.write(response.text)
        all_fasta.close()
    except FileNotFoundError:
        all_fasta.close()
        getFasta()
        
getFasta(uniqueID)