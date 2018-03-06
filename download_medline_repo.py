"""
	Download all of the medline files
	Author: William Goldsworthy
	Date: 5.3.18
"""

import os
import urllib.request

def download_medline_repo():
	print("Downloading Medline Repo");
	# At the time of writing there are 892 medline files available on the 2017 repo.
	url = "https://mbr.nlm.nih.gov/Download/Baselines/2017/medline17n0"
	for i in range(1, 892):
		num_len = len(str(i));
		if num_len == 1:
			filename = "00" + str(i) + ".xml.gz";
		elif num_len == 2:
			filename = "0" + str(i) + ".xml.gz";
		else:
			filename = str(i) + ".xml.gz";

		url_request = url + filename;
		print(url_request);
		u = urllib.request.urlopen(url_request);
		file_size = int(u.getheader('Content-Length'));
		print("Downloading: " + filename + " Bytes: " + str(file_size));


warning = """
This script will download all of the medline repo files for 2017.
This is approximately 30GB of data:
""";

print(warning);

x = input("Are you sure you want to continue? Y/N and return \n");

if x == "Y" or x == "y":
	download_medline_repo();
else:
	exit();