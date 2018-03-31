""" 
	Remove additional zipped files from data set
"""

import os
import sys

dataFolder = 'Scharrhud_Data/medline'
listfiles = os.listdir(dataFolder)

for medlineFile in listfiles:
	if medlineFile.endswith('.gz'):
		path = dataFolder + "/" + medlineFile
		os.remove(path);
	else:
		print("Not removing file: " + medlineFile);