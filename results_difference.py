# Percentage difference in results
# Author: William Goldsworthy
# Date: 3.4.18

import sys, getopt
import os
import difflib

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

results1 = open(config.gold, 'r');
results2 = open(config.results, 'r');

def formatter(results):
	result_list = []
	for line in results:
		words = line.split();
		result_list.append(words[1]);
	return result_list

results1 = formatter(results1);
results2 = formatter(results2);

sm = difflib.SequenceMatcher(None, results1, results2);
print(sm.ratio());