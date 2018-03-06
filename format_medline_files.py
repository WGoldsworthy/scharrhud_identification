"""
	Format the medline files from xml to individual text files
	Author: William Goldsworthy
	Date: 5.3.18
"""

import os
from Bio import Entrez

def fetch_medline_repo_files(filename):
    print("Fetching Medline Baseline Repo");
    filename = "Scharrhud_Data/medline/" + filename;
    medlineBaseRepo = open(filename, "r");
    testData = Entrez.read(medlineBaseRepo);
    testDataFolder = 'Scharrhud_Data/TestData/';

    for record in testData['PubmedArticle']:
    	pmid = record['MedlineCitation']['PMID'];

    	try: 
    		recordFileName = testDataFolder + str(pmid) + '.txt';
    		recordFile = open(recordFileName, "w");

    		if type(record['MedlineCitation']['Article']['Abstract']['AbstractText']) is list:
    			recordFile.write(record['MedlineCitation']['Article']['Abstract']['AbstractText'][0]);
    		else:
    			recordFile.write(record['MedlineCitation']['Article']['Abstract']['AbstractText']);
    	
    	except KeyError:
    		# Only want to include files that have an abstract
    		os.remove(recordFileName);
    print("Finished Fetching Medline files for " + filename);

medline_folder = "Scharrhud_Data/medline"
medline_files = os.listdir(medline_folder);

for medlinefile in medline_files:
	if medlinefile.endswith(".xml"):
		fetch_medline_repo_files(medlinefile);