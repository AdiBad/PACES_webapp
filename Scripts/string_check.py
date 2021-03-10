#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###################################################################
##      Checking if FASTA from interaction.py is complete        ##
###################################################################
##  - print all missing proteins                                 ## 
###################################################################

""" 
Requires: filteredData.tsv (generated with filter.py) & filteredDataSeq.fasta (generated with interaction.py)
Output: prints UniProt ID's of proteins of which the sequence could not be fetched
"""

import pandas as pd

data = pd.read_csv('../Output/filteredData.tsv',sep='\t')
uniqueID = data['Protein'].unique()

"""
Function to extract the gene name in a FASTA header line

Arguments:
    descripString: string in which gene should be located
"""
def findProtID(descripString):
    import re
    global result
    gene = re.search('\|(.+?)\|',descripString)
    if gene:
        result = gene.group(1)
    return result

"""
Function to extract all protein identifiers that are found in a multifasta file.
"""
def getAllFoundProteinIDs():
    result = []
    for line in open("../Output/filteredDataSeq.fasta"):
        global descrip
        if line.startswith(">"):
                descrip = line
                protID = findProtID(descrip)
                result.append(protID)
    return result

##Print all UniProt IDs that are in the filtered input data but are missing in the multifasta file
foundProt = getAllFoundProteinIDs()
print(set(foundProt) ^ set(uniqueID))