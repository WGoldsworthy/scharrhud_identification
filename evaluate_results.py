# Evaluation scores for precision, recall and F-measure
# Author: William Goldsworthy
# Date: 29.3.17

import sys, getopt
import os
from Bio import Entrez
from Bio import Medline

class CommandLine:
	def __init__(self):

		opts, args = getopt.getopt(sys.argv[1:], 'g:r:');
		opts = dict(opts);

		if '-g' in opts:
			self.gold = opts['-g'];
		else:
			print('Use -g flag to set gold standard file');
			sys.exit();

		if '-r' in opts:
			self.results = opts['-r'];
		else:
			print("Use -r flag to set results file");
			sys.exit();

config = CommandLine();

def create_relevant_list():
	Pids = []
	i=0
	idsfile = open("Scharrhud_Data/PM_ids.txt", "r");
	for pmid in idsfile:
		# Strip 'PM:' from the strings
		pmid = pmid[3:];
		Pids.append(pmid);


	Entrez.email = 'wjgoldsworthy1@sheffield.ac.uk';

	handle = Entrez.efetch(db="pubmed", id=Pids, rettype="medline", retmode="xml")

	# Use this for XML parsing
	records = Entrez.read(handle);

	gold_file = open("Scharrhud_Data/2013_docs.txt", "w");

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
				print("Less than 2012");
			else:
				pmid = record['MedlineCitation']['PMID'];
				gold_file.write(pmid);
				gold_file.write("\n");
				
		except KeyError:
			print('Error - No Abstract for record');


	handle.close();

def check_against_gold():
	gold = open(config.gold, "r");
	results = open(config.results, "r");

	recall = 0
	precision = 0

	A = 0;
	B = 0;
	C = 0;

	relevant_docs = [];
	for pid in gold:
		pid = pid.strip('\n');
		relevant_docs.append(pid);

	retrieved_docs = [];
	for line in results:
		words = line.split();
		retrieved_docs.append(words[1][:-4]);

	for relevant in relevant_docs:
		if relevant in retrieved_docs:
			A = A + 1;
		else:
			C = C + 1;

	for retrieved in retrieved_docs:
		if retrieved not in relevant_docs:
			B = B + 1;

	print("A: " + str(A))
	print("B: " + str(B))
	print("C: " + str(C))

	recall = A / (A + C);
	precision = A / (A + B);
	f_measure = (2 * precision * recall) / (precision + recall);

	print("Recall: " + str(recall) )
	print("Precision: " + str(precision) )
	print("F-Measure: " + str(f_measure) )

check_against_gold();
# create_relevant_list();