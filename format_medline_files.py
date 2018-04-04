"""
	Format the medline files from xml to individual text files
	Author: William Goldsworthy
	Date: 5.3.18
"""

import sys, getopt
import os
from Bio import Entrez

class CommandLine:
    def __init__(self):

        opts, args = getopt.getopt(sys.argv[1:], 'p:t:');
        opts = dict(opts);

        if '-p' in opts:
            self.path = opts['-p'];
        else:
            print("Use -p flag to set medline folder");
            sys.exit();

        if '-t' in opts:
            self.testDataPath = opts['-t'];
        else:
            print("Use the -t flag to set the test data folder path");
            sys.exit();

config = CommandLine();

def fetch_medline_repo_files(filename):
    print("Fetching Medline Baseline Repo" + filename);
    # filename = "Scharrhud_Data/medline/" + filename;
    filename = config.path + filename;
    medlineBaseRepo = open(filename, "r");
    testData = Entrez.read(medlineBaseRepo);

    # testDataFolder = 'Scharrhud_Data/TestData/';

    # testDataFolder = '/Volumes/EXTERNALDRV/TestData/';
    testDataFolder = config.testDataPath;

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

# medline_folder = "Scharrhud_Data/medline"
medline_folder = config.path;
medline_files = os.listdir(medline_folder);

for medlinefile in medline_files:
	if medlinefile.endswith(".xml"):
		fetch_medline_repo_files(medlinefile);