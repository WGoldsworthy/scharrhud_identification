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
import io

class CommandLine:
    def __init__(self):
        print("Scharrhud Document Identification.")

        opts, args = getopt.getopt(sys.argv[1:], 'hqo:r:p:')
        opts = dict(opts)
            
        if '-h' in opts:
            self.print_help();

        if '-q' in opts:
            self.fetch_query_abstracts = True;
        else:
            self.fetch_query_abstracts = False;

        if '-o' in opts:
            self.write_to_output_file = True;
            self.output_file = opts['-o']
        else:
            self.output_file = None

        if '-r' in opts:
            self.num_results = opts['-r']
        else:
            self.num_results = 50

        if '-p' in opts:
            self.test_path = opts['-p']
        else:
            self.test_path = 'Scharrhud_Data/test_classifier/';

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

# Used to control what functions run based on the commandline arguments
config = CommandLine();

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


def fetch_query_abstracts():

	handle = Entrez.efetch(db="pubmed", id=Pids, rettype="medline", retmode="xml")

	# Use this for XML parsing
	records = Entrez.read(handle);
	query_file = open('Scharrhud_Data/test_classifier/query/query.txt', "w");

	testDataFolder = 'Scharrhud_Data/test_classifier/TestData/'

	for record in records['PubmedArticle']:
		# print(record['MedlineCitation']['Article']['ArticleDate'][0]['Year']);
		try:
			abstractText = record['MedlineCitation']['Article']['Abstract']['AbstractText'];

			try:
				yearOfPub = record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year'];
			except:
				yearOfPub = record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['MedlineDate'];
				yearOfPub = yearOfPub[:4];


			if int(yearOfPub) < 2012:
				for string in abstractText:
					query_file.write(string);
			else:
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
				
		except KeyError:
			print('Error - No Abstract for record');

	handle.close();

if config.fetch_query_abstracts:
	fetch_query_abstracts();

# medlineBaseRepo = open('Scharrhud_Data/medline17n0001.xml', "r");
# testData = Entrez.read(medlineBaseRepo);

# for record in testData['PubmedArticle']:
# 	pmid = record['MedlineCitation']['PMID'];

# 	try: 
# 		recordFileName = testDataFolder + str(pmid) + '.txt';
# 		recordFile = open(recordFileName, "w");

# 		if type(record['MedlineCitation']['Article']['Abstract']['AbstractText']) is list:
# 			recordFile.write(record['MedlineCitation']['Article']['Abstract']['AbstractText'][0]);
# 		else:
# 			recordFile.write(record['MedlineCitation']['Article']['Abstract']['AbstractText']);
	
# 	except KeyError:
# 		# Only want to include files that have an abstract
# 		os.remove(recordFileName);


path = 'Scharrhud_Data/test_classifier/TestData/';

doc_ids = os.listdir(path);

documents = [io.open(path+f, 'rb').read() for f in doc_ids];

tfidf = vectorizer.fit_transform(documents);

query_file = open('Scharrhud_Data/test_classifier/query/query.txt', "r");

query = vectorizer.transform(query_file);

# feature_names = vectorizer.get_feature_names();

cosine_sim = linear_kernel(query, tfidf).flatten();

doc_sim_values = dict(zip(doc_ids, cosine_sim));

ranked = sorted(doc_sim_values.items(), key=operator.itemgetter(1), reverse=True);

# print(ranked[0:10]);

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