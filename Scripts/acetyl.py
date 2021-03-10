#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###################################################################
##         Extracting protein level acetylation info             ##
###################################################################
##  - reshape peptide level information to protein level         ## 
###################################################################

""" 
Requires: filteredData.tsv (generated with filter.py)
Output: acetylation.tsv
"""

import pandas as pd
import numpy as np

data = pd.read_csv('../Output/filteredData.tsv',sep='\t')
uniqueID = data['Protein'].unique()
dfP = pd.DataFrame(data=uniqueID,columns=["uniprotID"])


##Functions on the original (filtered) data: 

"""
Function to find out if this acetylation site was detected in both conditions, or one of both.
Conditions: + acetyltransferase gp13, control
Effect: returns this information in a new column

Arguments:
    df: DataFrame containing the filtered data (generated with filter.py)
"""

def acWhichConditions(df):
    i = 0
    for x in df['Intensity.L.']:
        if (x == 0) and (df.iloc[i]['Intensity.H.'] != 0):
            result = 'gp13'
        elif (x != 0) and (df.iloc[i]['Intensity.H.'] == 0):
            result = 'control'
        else:
            result = 'both'
        df.at[i,'Condition'] = result
        i = i + 1
        
"""
Function that applies a logarithm of base 2 of the Ratio.H.L column (providing a log fold change)
Effect: returns this information in a new column

IMPORTANT: we compare +gp13/control --> positive result: more acetylation in +gp13

Arguments:
    df: DataFrame containing the filtered data (generated with filter.py)
"""
def acLogFold(df):
    i = 0
    for x in df['Ratio.H.L.Normalized']:
        if np.isnan(x):
            df.at[i,'logFoldChange'] = np.NaN
        else:
            df.at[i,'logFoldChange'] = np.log2(x)
        i = i + 1
        
##Functions to generate acetylation.tsv: 

"""
Function that returns all rows of the filtered data that have as protein name the specified input string

Arguments:
    descripString: UniProt identifier for which one wants to fetch information
"""
def findRowsByProt(descripString):
    rowsWithProt = data[data['Protein']==descripString]
    return rowsWithProt
        
"""
Function to get the number of acetylated peptides for a specified input protein ID.
Reasoning: Row per acetylated peptide + only 1 acetylation site at each peptide detected --> #rows = #aectylation sites

Arguments:
    protId: UniProt identifier for protein of which one wants to calculate the number of acetylation sites
"""
def calcAcSites(protId):
    protDat = findRowsByProt(protId)
    return protDat.shape[0]    
        
"""
Function to extract the gene names based on a specified input protein ID, by making use of REGEX
Reasoning: if mulitple possibilities for gene name are given, the first one is assumed to be correct

Arguments:
    protId: UniProt identifier for protein of which one wants to calculate the number of acetylation sites
"""
def findGeneName(protId):
    import re
    global result
    
    protDat = findRowsByProt(protId)
    regionToSearch = protDat.iloc[0]['Protein.Descriptions']
    gene = re.search('GN=(.+?) PE',regionToSearch)
    
    if gene:
        result = gene.group(1)
    return result

"""
Function to check in which conditions the peptides of 1 protein are detected
Conditions: 'control', 'both' ('gp13' is never observed)
Effect: returns string stating either what the condition is for all peptides

Arguments:
    protId: UniProt identifier for protein of which one wants to calculate the number of acetylation sites
"""
def acCondProt(protId):
    protDat = findRowsByProt(protId)
    result = pd.unique(protDat['Condition'])
    if result.shape[0] == 1:
        return "all peptides: " + result[0]
    else:
        result = np.array(protDat['Condition'])
        peptString = ""
        for i in range(0,result.shape[0]):
            pept = "peptide " + str(i+1) + ": " + result[i] + " // "
            peptString = peptString + pept
        return peptString[:-4]

"""
Funtion to return the peptide sequences for a given protein ID
Effect: returns a string with all the peptides and their sequence, seperated by //

Arguments:
    protId: UniProt identifier for protein of which one wants to calculate the number of acetylation sites
"""
def getDetectedPeptides(protId):
    protDat = findRowsByProt(protId)
    result = pd.DataFrame(data=protDat['Modified.Sequence'])
    peptString = ""
    for i in range(0,result.shape[0]):
        pept = "peptide " + str(i+1) + ": " + result.iloc[i]['Modified.Sequence'] + " // "
        peptString = peptString + pept
    return peptString[:-4]
 
"""
Function to return the logFoldChanges for all peptides for a given protein ID.
Effect: returns a string with all the peptides and their log fold change, seperated by //

Arguments:
    protId: UniProt identifier for protein of which one wants to calculate the number of acetylation sites
"""
def getLogFoldPeptides(protId):
    protDat = findRowsByProt(protId)
    result = pd.DataFrame(data=protDat['logFoldChange'])
    peptString = ""
    for i in range(0,result.shape[0]):
        pept = "peptide " + str(i+1) + ": " + str(result.iloc[i]['logFoldChange']) + " // "
        peptString = peptString + pept
    return peptString[:-4]
    
"""
Function to see global picture in terms of increase/decrease of logFoldChange of peptides of a given protein.

Arguments:
    protId: UniProt identifier for protein of which one wants to calculate the number of acetylation sites
"""
def logFCProt(protId):
    protDat = findRowsByProt(protId)
    arr = np.array(protDat['logFoldChange'])
    arrNoNa = arr[np.logical_not(np.isnan(arr))]
    sign = np.sign(arrNoNa)
    
    appendage = ""
    if arr.shape[0] != arrNoNa.shape[0]:
        appendage = ' !NaN in peptide(s)'
    if np.unique(sign).shape[0] == 1:
        if "-" in str(sign[0]):
            result = 'negative'
        else:
            result = 'positive'
    else:
        if arrNoNa.shape[0] > 0:
            result = 'depends on peptide'
        if arrNoNa.shape[0] == 0:
            result =''
    return result+appendage
    
"""
Create extra columns in dataframe where you assign these values:
    gene name, number of acetylation sites, peptides, condition in which the peptide/protein was detected, peptide and protein level log fold change
    
Arguments:
    df: DataFrame containing the UniProt ID's (pandas DataFrame)
    descripCol: name of column in df that contains all UniProt ID's
"""
def makeExtraCol(df,descripCol):
    i = 0
    for x in df[descripCol]:
        df.at[i,"geneName"] = findGeneName(x)
        df.at[i,"numAcSites"] = calcAcSites(x)
        df.at[i,"peptides"] = getDetectedPeptides(x)
        df.at[i,"detectCondition"] = acCondProt(x)
        df.at[i,"peptLogFC"] = getLogFoldPeptides(x)
        df.at[i,"protLogFC"] = logFCProt(x)
        i = i+1   
    
##Applying functions to the filtered dataframe:
acWhichConditions(data)
acLogFold(data)

##Creating the new dataframe wchich contains acetylation information on protein level:
makeExtraCol(dfP,'uniprotID')
dfP.to_csv('../Output/acetylation.tsv',sep='\t')