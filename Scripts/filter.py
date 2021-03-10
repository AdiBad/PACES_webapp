#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###################################################################
##                    Filtering the input data                   ##
###################################################################
##  - remove peptides with no measured intensity, or PEP >= 0.05 ## 
###################################################################

""" 
Arguments: acetylation data (as tab seperated textfile)
Output: filteredData.tsv
"""

import sys
import pandas as pd

allData = pd.read_csv(str(sys.argv[1]),sep='\t',decimal=',')
filtData = allData[allData["PEP"] < 0.05]
filtData = filtData[filtData["Intensity."]!=0]
filtData.to_csv('../Output/filteredData.tsv', sep = '\t')

