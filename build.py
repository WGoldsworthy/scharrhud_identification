# Build Script for the system.
# Set up file structure and do any intial set up required.
# Author: William Goldsworthy
# Date: 4.4.18

import os
import sys
from shutil import copyfile

def build():
	if not os.path.exists('Scharrhud_Data/medline/'):
		os.makedirs('Scharrhud_Data/medline/');
		print('[INFO] Medline Folder Created.')
	else:
		print('[INFO] Medline Folder already exists.');

	if not os.path.exists('Scharrhud_Data/relevance_feedback/'):
		os.makedirs('Scharrhud_Data/relevance_feedback/');
		copyfile('Scharrhud_Data/query/query.txt', 'Scharrhud_Data/relevance_feedback/query.txt');
		print('[INFO] Relevance Feedback Folder Created.')
	else:
		print('[INFO] Relevance Feedback Folder already exists.')

	if not os.path.exists('Scharrhud_Data/relevance_feedback/medline/'):
		os.makedirs('Scharrhud_Data/relevance_feedback/medline/');
		print('[INFO] Relevance Feedback Medline Folder Created.')
	else:
		print('[INFO] Relevance Feedback Medline Folder already exists.')

	if not os.path.exists('Scharrhud_Data/relevance_feedback/results/'):
		os.makedirs('Scharrhud_Data/relevance_feedback/results/');
		print('[INFO] Relevance Feedback Results Folder Created.')
	else:
		print('[INFO] Relevance Feedback Results Folder already exists.')

	if not os.path.exists('Scharrhud_Data/relevance_feedback/TestData/'):
		os.makedirs('Scharrhud_Data/relevance_feedback/TestData/');
		print('[INFO] Relevance Feedback Test Data Folder Created.')
	else:
		print('[INFO] Relevance Feedback Test Data Folder already exists.')

	if not os.path.exists('Scharrhud_Data/test_classifier/'):
		os.makedirs('Scharrhud_Data/test_classifier/');
		print('[INFO] Test Classifier Folder Created.')
	else:
		print('[INFO] Test Classifier Folder already exists.')

	if not os.path.exists('Scharrhud_Data/test_classifier/query/'):
		os.makedirs('Scharrhud_Data/test_classifier/query/');
		copyfile('Scharrhud_Data/query/query.txt', 'Scharrhud_Data/test_classifier/query/query.txt');
		print('[INFO] Test Classifier Query Folder Created.')
	else:
		print('[INFO] Test Classifier Query Folder already exists.')

	if not os.path.exists('Scharrhud_Data/test_classifier/TestData/'):
		os.makedirs('Scharrhud_Data/test_classifier/TestData/');
		print('[INFO] Test Classifier Test Data Folder Created.')
	else:
		print('[INFO] Test Classifier Test Data Folder already exists.');

	if not os.path.exists('Scharrhud_Data/TestData/'):
		os.makedirs('Scharrhud_Data/TestData/');
		print('[INFO] Test Data Folder Created.')
	else:
		print('[INFO] Test Data Folder already exists.')

build();