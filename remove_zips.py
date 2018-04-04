""" 
	Remove additional zipped files from data set
"""

import os
import sys, getopt

class CommandLine:
	def __init__(self):

		opts, args = getopt.getopt(sys.argv[1:], 'p:');
		opts = dict(opts);

		if '-p' in opts:
			self.path = True
			self.path = opts['-p'];
		else:
			print('Set a path with -p flag');
			sys.exit();

config = CommandLine();


dataFolder = config.path
listfiles = os.listdir(dataFolder)

for medlineFile in listfiles:
	if medlineFile.endswith('.gz'):
		path = dataFolder + "/" + medlineFile
		os.remove(path);
	else:
		print("Not removing file: " + medlineFile);