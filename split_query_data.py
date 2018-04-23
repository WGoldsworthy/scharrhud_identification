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
import scipy.sparse
from scipy.sparse import csr_matrix
import pickle

class CommandLine:
	def __init__(self):
		print("Scharrhud Document Identification.")

		opts, args = getopt.getopt(sys.argv[1:], 'hqo:r:p:m:s:n')
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
			self.write_to_output_file = False;
			self.output_file = None

		if '-r' in opts:
			self.num_results = opts['-r']
		else:
			self.num_results = 50

		if '-p' in opts:
			self.test_path = opts['-p']
		else:
			self.test_path = 'Scharrhud_Data/test_classifier/';

		if '-n' in opts:
			self.save = True
		else:
			self.save = False;

		if '-m' in opts:
			self.preload = True;
			self.tfidf = opts['-m'];
		else:
			self.preload = False;

		if '-s' in opts:
			self.bySimilarity = True;
			self.simValue = opts['-s'];
		else:
			self.bySimilarity = False;

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

s_file = open('Pubmed_stopwords.txt','r')
stop_words = [];
for word in s_file:
	stop_words.append(word.strip());

stemmer = nltk.stem.lancaster.LancasterStemmer();
remove_punctuation_map = dict((ord(char), None) for char in string.punctuation);

def stem_tokens(tokens):
	# remove stop words
	filtered = [w for w in tokens if not w in stop_words]
	#remove numbers
	new_items = [item for item in filtered if not item.isdigit()]
	return [stemmer.stem(item) for item in tokens]

def normalize(text):
	return stem_tokens(nltk.word_tokenize(text.lower().translate(remove_punctuation_map)))

if not config.preload:
	vectorizer = TfidfVectorizer(tokenizer=normalize, stop_words=stop_words);

Entrez.email = 'wjgoldsworthy1@sheffield.ac.uk';

def save_sparse_csr(filename, array):
	np.savez(filename, data=array.data, indices=array.indices, indptr=array.indptr, shape=array.shape);

def load_sparse_csr(filename):
	loader = np.load(filename + '.npz');
	return csr_matrix((loader['data'], loader['indices'], loader['indptr']), shape=loader['shape']);


def fetch_query_abstracts():

	handle = Entrez.efetch(db="pubmed", id=Pids, rettype="medline", retmode="xml")

	# Use this for XML parsing
	records = Entrez.read(handle);
	query_file = open('Scharrhud_Data/test_classifier/query/query.txt', "w");
	pmids_2013_file = open('Scharrhud_Data/2013_docs.txt', 'w');

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


			if int(yearOfPub) < 2013:
				for string in abstractText:
					query_file.write(string);
			else:
				pmid = record['MedlineCitation']['PMID'];
				pmids_2013_file.write(pmid + '\n');
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

if config.preload:
	tfidf = load_sparse_csr(config.tfidf);
	vect_file = open('vectorizer.pk', 'rb');
	vectorizer = pickle.load(vect_file);
else:
	documents = [io.open(path+f, 'rb').read() for f in doc_ids];
	tfidf = vectorizer.fit_transform(documents);

	if config.save:
		with open('vectorizer.pk', 'wb') as fin:
			pickle.dump(vectorizer, fin);
		print("Finished calculating tfidf scores");
		save_sparse_csr('tfidf', tfidf);

query_file = open('Scharrhud_Data/test_classifier/query/query.txt', "r");

query = vectorizer.transform(query_file);

cosine_sim = linear_kernel(query, tfidf).flatten();

doc_sim_values = dict(zip(doc_ids, cosine_sim));

ranked = sorted(doc_sim_values.items(), key=operator.itemgetter(1), reverse=True);

# Function to write results to an output file specified in the commandline arguments
# As default this is off and results are printed to console.
def write_to_output_file(ranked):
	output_file = open(config.output_file, "w");
	if config.bySimilarity:
		x = 0;
		for item in ranked:
			if item[1] >= float(config.simValue):
				output_file.write(str(x + 1) + ": ");
				output_file.write(' '.join(map(str , ranked[x])));
				output_file.write('\n');
				x = x + 1;
			else:
				break;
	else:
		for x in range(0, int(config.num_results)):
			output_file.write(str(x + 1) + ": ");
			output_file.write(' '.join(map(str , ranked[x])));
			output_file.write('\n');

if config.write_to_output_file:
	write_to_output_file(ranked);
else:
	print(ranked[0:config.num_results]);