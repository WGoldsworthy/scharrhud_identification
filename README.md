# Medical Document Identification using text mining for the ScHARRHUD Database

This project aims to use text mining techniques to identify medical studies for the ScHARRHUD database at the University of Sheffield.

### Getting Started

Firstly clone or download the Repository.

```
git clone https://github.com/WGoldsworthy/scharrhud_identification.git
```

Once you have the repository locally, you need to run the build script which will set up the folders needed to house the data.

```
python3 build.py
```

With the folder structure set up, the next step is to acquire the test data. For this project, we were using the Medline Base Repository for only the years between 2005 and 2013. Using the download_medline_repo.py script, this will retrieve all the files and unzip them into youe chosen folder. By default these files will be downloaded to the Scharrhud_data/medline/ folder, however you can specifiy a different path using the -p flag if you need to. 

```
python3 download_medline_repo.py -p Scharrhud_Data/medline/
```

Once the files have downloaded, they need to be formatted into individual text files so that they can be passed into the system. This is a time consuming process and so this step was put into an individual script.

Use the -p flag to define the path to the medline files.
Use the -t flag to define the path to which you want the individual text files to be stored.

```
python3 format_medline_files.py -p Scharrhud_Data/medline/ -t Scharrhud_Data/TestData/
```

Now that all the data is available, the system should be ready to run.


### Prerequisites

