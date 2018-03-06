# Document Identification using text mining for Scharrhud database
# Author: William Goldsworthy
# Date: 6.11.17

import sys, getopt
import csv
from Bio import Entrez
import os
from Bio import Medline
from itertools import repeat
import nltk, string
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import linear_kernel
import operator
import numpy as np
import timeit

class CommandLine:
    def __init__(self):
        print("Scharrhud Document Identification.")

        opts, args = getopt.getopt(sys.argv[1:], 'hqmo:r:')
        opts = dict(opts)
            
        if '-h' in opts:
            self.print_help();

        if '-q' in opts:
            self.fetch_query_abstracts = True;
        else:
            self.fetch_query_abstracts = False;

        if '-m' in opts:
            self.fetch_medline = True;
        else: 
            self.fetch_medline = False;

        if '-o' in opts:
            self.write_to_output_file = True;
            self.output_file = opts['-o']
        else:
            self.output_file = None

        if '-r' in opts:
            self.num_results = opts['-r']
        else:
            self.num_results = 50

    def print_help(self):
        help = """\
            ------------------------------------------------------------
            USE: python <PROGNAME> (options)
            OPTIONS:
                -h : print this help message
                -q : fetch the query abstracts and write to a query file. 
                     If already have the abstracts and file, this can be turned off to reduce runtime.
                -m : fetch the Medline Baseline Repo. These are all the documents
                     that will be tested against. If you need to redownload and process all of the files
                     into individual text files then I recommend using the format_medline_files.py script.
                -o FILE : write results to an output FILE (Default: stdout)
                -r NUMBER : Number of top results you want to return (Defualt: 50)
            ------------------------------------------------------------\
            """
        print(help);
        sys.exit()

    def fetch_query(self):
        print("Fetch query Abstract");
        fetch_query_abstracts();

start = timeit.default_timer();

# Used to control what functions run based on the commandline arguments
config = CommandLine();

# Sometime nltk needs the resource 'punkt' to use word_normalize function
# Can use the line below to download and enable
# nltk.download('punkt')

# The text file contains all of the Pids of the documents that are part of the
# ScharrHud database. This extracts them and puts them into a list.

Pids = []
i=0
idsfile = open("Scharrhud_Data/PM_ids.txt", "r");
for pmid in idsfile:
	# Strip 'PM:' from the strings
	pmid = pmid[3:];
	Pids.append(pmid);

stemmer = nltk.stem.lancaster.LancasterStemmer();
remove_punctuation_map = dict((ord(char), None) for char in string.punctuation);

def stem_tokens(tokens):
    return [stemmer.stem(item) for item in tokens]

def normalize(text):
    return stem_tokens(nltk.word_tokenize(text.lower().translate(remove_punctuation_map)))

vectorizer = TfidfVectorizer(tokenizer=normalize, stop_words='english');

Entrez.email = 'wjgoldsworthy1@sheffield.ac.uk';


# Fetches the abstracts fro each of the documents in the scharrhud database (query).
# Writes all abstracts into a single file which will be used as the query for the cosine similarity
def fetch_query_abstracts():
    print("Fetching Query Abstracts");
    handle = Entrez.efetch(db="pubmed", id=Pids, rettype="medline", retmode="xml")

    # Use this for XML parsing
    records = Entrez.read(handle);

    query_file = open('Scharrhud_Data/query/query.txt', "w");

    for record in records['PubmedArticle']:
    	try:
    		abstractText = record['MedlineCitation']['Article']['Abstract']['AbstractText'];

    		for string in abstractText:
    			query_file.write(string);
    			
    	except KeyError:
    		print('Error - No Abstract for record');


    handle.close();
    print("Finished Fetching Query Abstracts");


if config.fetch_query_abstracts:
    fetch_query_abstracts();

# Section for handling the Test data i.e. The documents that the query will be checked against.

# TODO: move this into a different file. require_relative
# The function should parse all medline xml files not just a single.
def fetch_medline_repo_files():
    print("Fetching Medline Baseline Repo");
    medlineBaseRepo = open('Scharrhud_Data/medline17n0001.xml', "r");
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
    print("Finished Fetching Medline Baseline Repo");

if config.fetch_medline:
    fetch_medline_repo_files();

path = 'Scharrhud_Data/TestData/';

doc_ids = os.listdir(path);

num_docs = len(doc_ids);
print("There are " + str(num_docs) + " documents in the test data.");

documents = [open(path+f, 'r').read() for f in doc_ids];

print("Calculating tfidf scores");
tfidf = vectorizer.fit_transform(documents);
print("Finished calculating tfidf scores");

query_file = open('Scharrhud_Data/query/query.txt', "r");

query = vectorizer.transform(query_file);

feature_names = vectorizer.get_feature_names();

print("Calculating cosine similarity values");
cosine_sim = linear_kernel(query, tfidf).flatten();
print("Finished calculating cosine similarity scores");

doc_sim_values = dict(zip(doc_ids, cosine_sim));

ranked = sorted(doc_sim_values.items(), key=operator.itemgetter(1), reverse=True);

# Function to write results to an output file specified in the commandline arguments
# As default this is off and results are printed to console.
def write_to_output_file(ranked):
    output_file = open(config.output_file, "w");
    for x in range(0, int(config.num_results)):
        output_file.write(str(x + 1) + ": ");
        output_file.write(' '.join(map(str , ranked[x])));
        output_file.write('\n');

if config.write_to_output_file:
    write_to_output_file(ranked);
else:
    print(ranked[0:config.num_results]);

# Function to display tfidf scores of the features. 
# Not sure whether to add this as an option, results can be very large and mostly meaningless
def display_scores(vectorizer, tfidf_result):
    
    scores = zip(vectorizer.get_feature_names(),
                 np.asarray(tfidf_result.sum(axis=0)).ravel())
    sorted_scores = sorted(scores, key=lambda x: x[1])
    for item in sorted_scores:
        print("{0:50} Score: {1}".format(item[0], item[1]));

end = timeit.default_timer();

print("Elapsed time: " + str( (end - start)/60) )

# display_scores(vectorizer, tfidf);