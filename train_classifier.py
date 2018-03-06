import sys, getopt
import csv
from ui_helper import UI
from Bio import Entrez
import os
from Bio import Medline
from itertools import repeat
import nltk, string
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import linear_kernel
import operator
import numpy as np

ui = UI();

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

# handle = Entrez.efetch(db="pubmed", id=Pids, rettype="medline", retmode="xml")

# # Use this for XML parsing
# records = Entrez.read(handle);

# # abstracts_folder = 'Scharrhud_Data/query';


# query_file = open('Scharrhud_Data/test_classifier/query/query.txt', "w");

# vectorizer = TfidfVectorizer(tokenizer=normalize, stop_words='english');

testDataFolder = 'Scharrhud_Data/test_classifier/TestData/'

# for record in records['PubmedArticle']:
# 	# print(record['MedlineCitation']['Article']['ArticleDate'][0]['Year']);
# 	try:
# 		abstractText = record['MedlineCitation']['Article']['Abstract']['AbstractText'];

# 		try:
# 			yearOfPub = record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year'];
# 		except:
# 			yearOfPub = record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['MedlineDate'];
# 			yearOfPub = yearOfPub[:4];


# 		if int(yearOfPub) < 2012:
# 			for string in abstractText:
# 				query_file.write(string);
# 		else:
# 			pmid = record['MedlineCitation']['PMID'];
# 			try: 
# 				recordFileName = testDataFolder + str(pmid) + '.txt';
# 				recordFile = open(recordFileName, "w");

# 				if type(record['MedlineCitation']['Article']['Abstract']['AbstractText']) is list:
# 					recordFile.write(record['MedlineCitation']['Article']['Abstract']['AbstractText'][0]);
# 				else:
# 					recordFile.write(record['MedlineCitation']['Article']['Abstract']['AbstractText']);
			
# 			except KeyError:
# 				# Only want to include files that have an abstract
# 				os.remove(recordFileName);
			
# 	except KeyError:
# 		print('Error - No Abstract for record');


# handle.close();

medlineBaseRepo = open('Scharrhud_Data/medline17n0001.xml', "r");
testData = Entrez.read(medlineBaseRepo);

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


path = 'Scharrhud_Data/test_classifier/TestData/';

doc_ids = os.listdir(path);

documents = [open(path+f, 'r').read() for f in doc_ids];

tfidf = vectorizer.fit_transform(documents);

query_file = open('Scharrhud_Data/query/query.txt', "r");

query = vectorizer.transform(query_file);

feature_names = vectorizer.get_feature_names();

cosine_sim = linear_kernel(query, tfidf).flatten();

doc_sim_values = dict(zip(doc_ids, cosine_sim));

ranked = sorted(doc_sim_values.items(), key=operator.itemgetter(1), reverse=True);

print(ranked[0:10]);

def display_scores(vectorizer, tfidf_result):
    
    scores = zip(vectorizer.get_feature_names(),
                 np.asarray(tfidf_result.sum(axis=0)).ravel())
    sorted_scores = sorted(scores, key=lambda x: x[1])
    for item in sorted_scores:
        print("{0:50} Score: {1}".format(item[0], item[1]));



# display_scores(vectorizer, tfidf);