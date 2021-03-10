#!/usr/bin/env python3
# -*- coding: utf-8 -*-

##############################################################
##  Making nodeDf dataframe                                 ##
##############################################################
##  - contains all attribute values for nodes and edges     ##
##  - contains inforation for interaction table             ##
##############################################################

""" 
Requires:
    string_interactions.tsv (manually generated)
    287.protein.actions.v11.0.txt.gz (given)
    string_mapping.tsv (manually generated)
    pathways.tsv (generated with kegg.py)
    acetylation.tsv (generated with acetyl.py)
Output: 
    nodeDf.tsv
    ackegg.tsv
"""

import pandas as pd
import numpy as np

# nodeDf = node dataframe
# This will be used to make nodes and edges for the graph-table, and to make the interaction table
# Data in this dataframe will be set as attributes of nodes and edges, to color the nodes, give labels, etc.
nodeDf = pd.read_csv('../String_man/string_interactions.tsv', sep='\t')

# Read in file with reactions information
interactions = pd.read_csv('../287.protein.actions.v11.0.txt.gz', sep='\t')

# Remove unnecessary columns from nodeDf
drop = ['neighborhood_on_chromosome', 'gene_fusion', 'phylogenetic_cooccurrence', 'homology', 'coexpression', 
                        'experimentally_determined_interaction', 'database_annotated', 'automated_textmining']
for col in drop:
    del nodeDf[col]

# merge nodeDf and interactions dataframes by common string id
merged_df = nodeDf.merge(interactions, how='left', left_on=['node1_string_id', 'node2_string_id'], right_on=['item_id_a', 'item_id_b'])

# Remove unnecessary columns from merged_df
drop = ['action', 'is_directional', 'a_is_acting', 'score', 'item_id_a', 'item_id_b']
for col in drop:
    del merged_df[col]

# rename column from 'mode' to 'interaction'
merged_df.rename(columns={'mode': 'interaction'}, inplace=True)

# remove all identical rows (across all columns)
merged_df.drop_duplicates(inplace=True)             # 13971 --> 11077

# replace 'NaN' with 'unknown'
merged_df.interaction.replace(np.NaN, 'unknown', inplace=True)

"""
The merged_df has almost identical rows where the only difference is the interaction, since two proteins can have multiple ways of interacting. In the following step, these different rows are merged together by grouping using all columns except for the interaction type, making this interaction into a set (= list with unique values).
This is then concatenated into a string, comma seperated.

Illustration:
The first dataframe will be transformed into the second one.
-----------------------------------------------------------------------------------------------------
| node 1    | node 2    | node1_string_id   | node2_string_id   | combined_score    | interaction   |
|----------------------------------------------------------------------------------------------------
| Acad10      acsA1	      287.DR97_5620	      287.DR97_1056	      0.671	              binding       |
| Acad10      DR97_149    287.DR97_5620       287.DR97_149        0.811               binding       |
| Acad10      DR97_149    287.DR97_5620       287.DR97_149        0.811               reaction      |
| ycgB	      ygaU	      287.DR97_3555	      287.DR97_2546	      0.590	              unknown       |
-----------------------------------------------------------------------------------------------------

---------------------------------------------------------------------------------------------------------
| node 1    | node 2    | node1_string_id   | node2_string_id   | combined_score    | interaction       |
|--------------------------------------------------------------------------------------------------------
| Acad10      acsA1	      287.DR97_5620	      287.DR97_1056	      0.671	              binding           |
| Acad10      DR97_149    287.DR97_5620       287.DR97_149        0.811               binding, reaction |
| ycgB	      ygaU	      287.DR97_3555	      287.DR97_2546	      0.590	              unknown           |
---------------------------------------------------------------------------------------------------------
"""
merged_df_ag = merged_df.groupby(['node1', 'node2', 'node1_string_id', 'node2_string_id', 'combined_score'])['interaction'].apply(set).reset_index()
merged_df_ag['interaction'] = [', '.join(map(str, s)) for s in merged_df_ag['interaction']]

## Making dictionaries
# Making uniprotID dictionary --> 0(1)
uniprot_info = pd.read_csv('../String_man/string_mapping.tsv', sep='\t')
uniprot_dict = dict(zip(uniprot_info['stringId'], uniprot_info['queryItem'].str.extract(r'(?<=\|)(.*)(?=\|)', expand=False)))
# Making kegg_id dictionary --> 0(1)
kegg_info = pd.read_csv('../Output/pathways.tsv',sep='\t')
kegg_dict = dict(zip(kegg_info['uniprotID'], kegg_info['keggID']))


# Using dictionaries to add identifiers for uniprot and kegg to the nodeDf
# This was not added before merging the interactions into one string, as some KEGG ID's have no value
# The groupby method from pandas would remove these rows, which is undesirable. Thus it is added after this step
merged_df_ag['node1_uniprot'] = merged_df_ag['node1_string_id'].apply(lambda x: uniprot_dict.get(x))
merged_df_ag['node2_uniprot'] = merged_df_ag['node2_string_id'].apply(lambda x: uniprot_dict.get(x))
merged_df_ag['node1_kegg'] = merged_df_ag['node1_uniprot'].apply(lambda x: kegg_dict.get(x))
merged_df_ag['node2_kegg'] = merged_df_ag['node2_uniprot'].apply(lambda x: kegg_dict.get(x))

# write nodeDf (merged_df_ag) dataframe to file
merged_df_ag.to_csv('../Output/nodeDf.tsv', sep='\t', index = False)


######################################################################################
##  Adding KEGG pathways info to acetylation dataframe                              ##
######################################################################################
##  - Data from acetylation.tsv will be used for acetylation table                  ##
##  - KEGG pathways info necessary to be able to filter individual proteins on this ##
######################################################################################

# Read in acetylation data
acetylation_pre = pd.read_csv('../Output/acetylation.tsv', sep='\t')

# Add KEGG pathway data by merging
acetylation = acetylation_pre.merge(kegg_info, how='left', left_on='uniprotID', right_on='uniprotID')

# removing rows with indexes
del acetylation['Unnamed: 0_x']
del acetylation['Unnamed: 0_y']

# write to file
acetylation.to_csv('../Output/ackegg.tsv', sep='\t', index=False)