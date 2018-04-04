"""
	Download all of the medline files
	Author: William Goldsworthy
	Date: 5.3.18
"""

import os
import urllib.request
import sys, getopt
import gzip
import shutil

class CommandLine:
	def __init__(self):

		opts, args = getopt.getopt(sys.argv[1:], 'p:')
		opts = dict(opts);

		if '-p' in opts:
			self.path = True;
			self.path = opts['-p'];
		else:
			self.path = False;

config = CommandLine();

def download_medline_repo(config):
	print("Downloading Medline Repo");
	# At the time of writing there are 892 medline files available on the 2017 repo.
	
	if config.path:
		path = config.path;
	else:
		path = 'Scharrhud_Data/medline/';

	url = "https://mbr.nlm.nih.gov/Download/Baselines/2017/medline17n0"
	# Range was changed here to only download years 2008-2015 (580-800)
	for i in range(580, 800):
		num_len = len(str(i));
		if num_len == 1:
			filename = "00" + str(i) + ".xml.gz";
		elif num_len == 2:
			filename = "0" + str(i) + ".xml.gz";
		else:
			filename = str(i) + ".xml.gz";

		url_request = url + filename;
		print(url_request);
		filename = path + 'medline17n0' + filename;
		f = open(filename, 'wb');
		u = urllib.request.urlopen(url_request);
		file_size = int(u.getheader('Content-Length'));
		print("Downloading: " + filename + " Bytes: " + str(file_size));


		file_size_dl = 0
		block_sz = 8192
		while True:
		    buffer = u.read()
		    if not buffer:
		        break

		    file_size_dl += len(buffer)
		    f.write(buffer)
		    status = r"%10d  [%3.2f%%]" % (file_size_dl, file_size_dl * 100. / file_size)
		    status = status + chr(8)*(len(status)+1)
		    print(status),

		f.close()

		print("Unzipping: " + filename);
		with gzip.open(filename, 'rb') as zipfile:
			unzipFilename = filename[:-3];
			with open(unzipFilename, 'wb') as unzipfile:
				shutil.copyfileobj(zipfile, unzipfile);


warning = """
This script will download all of the medline repo files for 2017.
This is approximately 30GB of data:
""";

print(warning);

if config.path:
	print("Your selected path is: " + config.path);
else:
	print("You have not specified a path with the -p flag. The data will be downloaded" +
		  " to Scharrhud_Data/medline/ as default");

x = input("Are you sure you want to continue? Y/N and return \n");

if x == "Y" or x == "y":
	download_medline_repo(config);
else:
	exit();