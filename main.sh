#!/bin/bash

#####################################################################################################
## Title:   Visualization of the acetylome changes in Pseudomonas aeruginosa upon phage infection  ##
## Authors: Stefaan Verwimp, Aditya Badola, Hannelore Longin, Ben De Maesschalck                   ##
## Part: Preprocessing pipeline - automatic                                                        ##
#####################################################################################################

cd Scripts
python filter.py ../input.txt
wait
python kegg.py 
python acetyl.py 
python interaction.py

