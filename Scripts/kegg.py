#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###################################################################
##             Annotating the proteins with KEGG                 ##
###################################################################
##  - fetch KEGG ID and KEGG pathways for all proteins in input  ## 
###################################################################

""" 
Requires: filteredData.tsv (generated with filter.py)
Output: pathways.tsv
"""

import pandas as pd
from bioservices.kegg import KEGG

data = pd.read_csv('../Output/filteredData.tsv',sep='\t')
uniqueID = data['Protein'].unique()
df = pd.DataFrame(data=uniqueID,columns=["uniprotID"])

k = KEGG()

"""
Function to get Kegg ID's when given a specific UniProt ID and a KEGG accessionpoint.
Makes use of KEGG.conv() function which converts foreign IDs to KEGG IDs.

Arguments:
    upId: UniProt ID (str)
    keggLink: accession point to KEGG database
"""
def getKeggId(upId,keggLink):
    result = keggLink.conv("pae","up:{}".format(upId))
    if result == "\n":
        keggId = "NA"
    else:
        for y in result.values():
            keggId = y.partition(":")[2]
    return keggId

"""
Function to iteratively apply getKeggId() to the values of a column containing uniprotID's
& then fill in the found KEGG IDs in a new column "keggID" of the same dataframe.

Arguments:
    df: DataFrame containing the UniProt ID's (pandas DataFrame)
    upIdCol: name of the column that conations the UniProt ID'd (string)
    keggLink: accession point to KEGG database
"""
def makeKeggCol(df,upIdCol,keggLink):
    i = 0
    for x in df[upIdCol]:
        df.at[i,"keggID"] = getKeggId(x,keggLink)
        i = i+1
"""
Function to retrieve the KEGG pathway, given a specific KEGG ID.
Makes use of KEGG.get_pathway_by_gene()

Arguments:
    kId: KEGG ID (string)
    keggLink: accession point to KEGG database
"""
def getKeggPath(kId,keggLink):
    if kId == "NA":
        result = "No KEGG ID given"
    else:
        path = keggLink.get_pathway_by_gene(kId,"pae")
        if path == None:
            result = "NA"
        else:
            result = path
    return result

"""
Function to parse a pathway dictionnary into a string which can be stored in a dataframe later.
When no pathway is retrieved returns: "No pathways"
Else: pathways in pathAccesionNumber:pathDescription format, seperating multiple pathways with  // 

Arguments: 
    diction: dictionary format pathway data (generated with getKeggPath)
"""
def parsePath(diction):
    global result
    if diction == "NA" or diction == "No KEGG ID given":
        result = "No pathways"
    else: 
        result = ""
        for (y,z) in zip(diction.values(),diction.keys()):
            result += z
            result += ":"
            result += y.lower()
            result += " // "
        result = result[:-4]
    return result

"""
Function to apply getKeggPath() and parsePath() iteratively to the values of a column containing KEGG IDs
& then fill in the parsed pathways in a new column "keggPathways" in the same dataframe

Arguments:
    df : dataframe with column containing KEGG identifiers
    kIdCol : name of column in df that contains KEGG identifiers
    keggLink: accession point to KEGG database
"""
def makeKeggPathCol(df,kIdCol,keggLink):
    i = 0
    for x in df[kIdCol]:
        df.at[i,"keggPathways"] = parsePath(getKeggPath(x,keggLink))
        i = i+1
        
        
makeKeggCol(df,"uniprotID",k)
makeKeggPathCol(df,"keggID",k)
df.to_csv('../Output/pathways.tsv', sep = '\t')