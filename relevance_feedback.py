# Relevance Feedback Implementation
# Author: William Goldsworthy
# Date: 2.4.18

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
import shutil
import scipy.sparse
from scipy.sparse import csr_matrix
import pickle

class CommandLine:
	def __init__(self):
		print("Scharrhud Document Identification.")

		opts, args = getopt.getopt(sys.argv[1:], 'hq:o:r:p:i:m:s:nt:')
		opts = dict(opts)
			
		if '-h' in opts:
			self.print_help();

		# CHANGED TO SPECIFY QUERY FILE
		if '-q' in opts:
			self.query = opts['-q']
			#self.fetch_query_abstracts = True;
		else:
			#self.fetch_query_abstracts = False;
			self.query = 'Scharrhud_Data/relevance_feedback/query.txt';

		if '-o' in opts:
			self.write_to_output_file = True;
			self.output_file = opts['-o']
		else:
			self.write_to_output_file = False;
			self.output_file = None

		if '-r' in opts:
			self.num_results = int(opts['-r'])
		else:
			self.num_results = 50

		if '-p' in opts:
			self.test_path = opts['-p']
		else:
			self.test_path = 'Scharrhud_Data/relevance_feedback/TestData/';

		if '-i' in opts:
			self.iterations = opts['-i']
		else:
			self.iterations = 1;

		if '-n' in opts:
			self.save = True;
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

		if '-t' in opts:
			self.termsToAdd = True;
			self.terms = opts['-t'];
		else:
			self.terms = 20;


	def print_help(self):
		help = """\
			------------------------------------------------------------
			USE: python <PROGNAME> (options)
			OPTIONS:
				-h : print this help message
				-q : fetch the query abstracts and write to a query file. 
					 If already have the abstracts and file, this can be turned off to reduce runtime.
				-o FILE : write results to an output FILE (Default: stdout)
				-r NUMBER : Number of top results you want to return (Defualt: 50)
				-p PATH : The path to your test data (Default: Scharrhud_Data/TestData/)
			------------------------------------------------------------\
			"""
		print(help);
		sys.exit()

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

s_file = open('PubMed_stopwords.txt','r')
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


def improve_query(results, iteration):
	doc_names = [x[0] for x in results]
	docs = [open(path+f, 'r').read() for f in doc_names]
	
	vectorizer2 = TfidfVectorizer(tokenizer=normalize, stop_words=stop_words);
	terms = vectorizer2.fit_transform(docs);
	
	feature_names = vectorizer2.get_feature_names();
	scores = zip(feature_names, np.asarray(terms.sum(axis=0)).ravel())
	sorted_scores = sorted(scores, key=lambda x: x[1])
	
	with open('Scharrhud_Data/relevance_feedback/results/query' + str(iteration) + '.txt', 'r') as query_file:
		query = query_file.read();

	new_query = open('Scharrhud_Data/relevance_feedback/results/query' + str(iteration_num + 1) + '.txt', 'w');
	new_query.write(query);

	numTerms = int(config.terms) * (-1);
	for item in sorted_scores[numTerms:]:
		new_query.write(item[0]);
		new_query.write(' ');

	new_query.close();
	query_file.close();
	
	# for item in sorted_scores:
	#  print("{0:50} Score: {1}".format(item[0], item[1]));
	# terms_file = open('Scharrhud_Data/relevance_feedback/results/terms_to_add.txt', "w");
	# for item in sorted_scores[-300:]:
	# 	terms_file.write(item[0]);
	# 	terms_file.write('\n');

# Fetches the abstracts fro each of the documents in the scharrhud database (query).
# Writes all abstracts into a single file which will be used as the query for the cosine similarity
def fetch_query_abstracts():
	print("Fetching Query Abstracts");
	handle = Entrez.efetch(db="pubmed", id=Pids, rettype="medline", retmode="xml")

	# Use this for XML parsing
	records = Entrez.read(handle);

	query_file = open('Scharrhud_Data/relevance_feedback/query.txt', "w");

	for record in records['PubmedArticle']:
		try:
			abstractText = record['MedlineCitation']['Article']['Abstract']['AbstractText'];

			for string in abstractText:
				query_file.write(string);
				
		except KeyError:
			print('Error - No Abstract for record');


	handle.close();
	print("Finished Fetching Query Abstracts");


# if config.fetch_query_abstracts:
# 	fetch_query_abstracts();

if config.test_path:
	path = config.test_path;

iteration_num = 0;
# shutil.copy2('Scharrhud_Data/relevance_feedback/query.txt', 'Scharrhud_Data/relevance_feedback/results/query0.txt');
# shutil.copy2('Scharrhud_Data/relevance_feedback/results/query200.txt', 'Scharrhud_Data/relevance_feedback/results/query0.txt');

shutil.copy2('Scharrhud_Data/test_classifier/query/query.txt', 'Scharrhud_Data/relevance_feedback/results/query0.txt');

doc_ids = os.listdir(path);

num_docs = len(doc_ids);
print("There are " + str(num_docs) + " documents in the test data.");

while int(config.iterations) > 0:

	print('[ITERATION NUMBER: ' + str(iteration_num) +']');

	if config.preload:
		print("Beginning Load of saved Matrix")
		tfidf = load_sparse_csr(config.tfidf);
		vect_file = open('vectorizer.pk', 'rb');
		vectorizer = pickle.load(vect_file)
		print("Loaded Saved Matrix");
	else:
		print("Calculating tfidf scores");
		documents = [open(path+f, 'r').read() for f in doc_ids];
		tfidf = vectorizer.fit_transform(documents);


		if config.save:
			print("Saving tfidf objects")
			with open('vectorizer.pk', 'wb') as fin:
				pickle.dump(vectorizer, fin);
			print("Finished calculating tfidf scores");
			save_sparse_csr('tfidf', tfidf);

	
	query_file = open('Scharrhud_Data/relevance_feedback/results/query' + str(iteration_num) +'.txt', "r");

	query = vectorizer.transform(query_file);

	feature_names = vectorizer.get_feature_names();

	print("Calculating cosine similarity values");
	cosine_sim = linear_kernel(query, tfidf).flatten();
	print("Finished calculating cosine similarity scores");

	doc_sim_values = dict(zip(doc_ids, cosine_sim));

	ranked = sorted(doc_sim_values.items(), key=operator.itemgetter(1), reverse=True);

	improve_query(ranked[0:config.num_results], iteration_num);

	config.iterations = int(config.iterations) - 1;
	iteration_num += 1;



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