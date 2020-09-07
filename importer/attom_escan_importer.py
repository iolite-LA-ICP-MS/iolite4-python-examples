# A python-based importer for iolite 4 starts with some metadata
#/ Type: Importer
#/ Name: Attom Peak EScan Importer
#/ Authors: Joe Petrus and Bence Paul
#/ Description: An importer for loading AttoM Peak EScan .csv files
#/ References: None
#/ Version: 1.0
#/ Contact: support@iolite-software.com

import time
import numpy as np
import pandas as pd
import re
from datetime import datetime


"""

NOTE: a valid importer needs to have the following functions to work:

setFileName(fileName)
correct_format()
import_data()

These three functions will be called by iolite at import time.
You can also add the following function to make selecting your files easier:

accepted_files()

This function should return a string with file extensions accepted by this
importer, e.g. ".csv". Can also return a semi-colon separated list,
e.g. ".csv;.dat").

Other than that, you can add whatever you like to your importer.

Qt imports can be done through 'iolite', e.g.
from iolite.QtCore import QRegularExpression

Before the functions are called a few additional objects are
added to the module

data        an interface to iolite's C++ data. E.g. you can get
            existing time series data or selection groups with it
            as well as make new ones.

IoLog       an interface to iolite's logging facility. You can add
            messages with, e.g., IoLog.debug('My message')

importer    an interface to the corresponding C++ class in iolite
            from which some properties can be accessed, e.g. fileName
"""

m_fileName = ""

def setFileName(fileName):
    m_fileName = fileName

def accepted_files():
    return ".csv"


def correct_format():
    """
    This method will be called by iolite when the user selects a file
    to import. Typically, it uses the provided name (stored in
    plugin.fileName) and parses as much of it as necessary to determine
    if this importer is appropriate to import the data. For example,
    although X Series II and Agilent data are both comma separated
    value files, they can be distinguished by the characteristic
    formatting in each. In our implementation, distinguishing the
    two is done with 'regular expressions'parsing of the first 
    several lines of the file.

    Keep in mind that there is nothing stopping you from just
    returning True (thus supporting all files!) or simply checking
    the file extension, but such generic checks can yield unexpected
    results. You cannot be sure which order the various importer
    plugins will be checked for compatibility.

    This method must return either True or False.
    """
    IoLog.debug("correct_format for Attom Peak EScan .csv file importer called on file = %s"%(importer.fileName))

    if importer.fileName.endswith('csv'):
        nrf_regex = r"Analysis Type: \"PeakEScan\""
        with open(importer.fileName) as f:
            match = re.search(nrf_regex, f.read(), re.MULTILINE)
            if match:
                return True

    return False

def import_data():
    """
    This method uses the provided file name (stored in plugin.fileName),
    parses its contents, and registers time series data with iolite

    Importer progress can be updated via the 'message' and 'progress'
    signals. These will be displayed in the iolite interface.

    When finished, the 'finished' signal should be emitted.
    """
    IoLog.debug("import_data called on file = %s"%(importer.fileName))

    #find the date and time using a regular expression:
    date_regex = r"Date: \"(?P<hours>\d{2}):(?P<minutes>\d{2}):(?P<seconds>\d{2})\s+(?P<day>\d{2})\/(?P<month>\d{2})\/(?P<year>\d{4})\""

    #find the sample name with a regex too:
    name_regex = r'Run:\s+"([\w._-]+)"$'

    with open(importer.fileName) as f:
        file_contents = f.read()
        match_date = re.search(date_regex, file_contents, re.MULTILINE)
        match_name = re.search(name_regex, file_contents, re.MULTILINE)

    if not match_date:
        IoLog.error("Couldn't find a match for the date time stamp")
        importer.message('Error during import')
        importer.progress(100)
        importer.finished()

    timestring = ",".join(match_date.groups())
    start_time = datetime.strptime(timestring, '%H,%M,%S,%d,%m,%Y')
    start_time_in_s = start_time.timestamp()
    IoLog.debug(f'Start time is: {start_time}')

    sample_name = ""
    if not match_name:
        IoLog.warning("Couldn't find a match for the sample name for file " + importer.fileName)
        sample_name = "Not found"
    else:
        sample_name = match_name.group(1).lstrip()
        IoLog.debug("Sample name: " + sample_name)

    #look through file to find the first line of data, which in these files
    # is preceded by a line with "========" written in it.
    first_line_of_data = 0
    names_list = []
    with open(importer.fileName, 'r') as file:
        for line in file:
            if "========" in line:
                next(file)                   #skip blank line
                names_list = next(file).split(",")
                break
            
            first_line_of_data += 1

    first_line_of_data += 4             #the first line of data is actually 4 lines below the ==== line

    IoLog.debug(f"Found first line of data on line: {first_line_of_data}")
    IoLog.debug(f"Found these headers: {names_list}")

    #Read csv file into pandas DataFrame 
    df = pd.read_csv(importer.fileName, skiprows=first_line_of_data, header=None, names=names_list)

    #add start time in seconds to time column:
    df['time'] = df['Timestamp'].add(start_time_in_s)
    # Add blank row at end of dataframe here to create a blank between each file
    df = df.append(pd.Series(), ignore_index=True)
    # And add 1 millisecond to last time for the timing of this blank value
    df['time'].iloc[-1] = df['time'].iloc[-2] + 0.001

    #Set the Step at this last blank row to 3 so that it will be included when we subset the table to just Step == 3 below
    df['Step'].iloc[-1] = 3

    #Drop "Mode.." and "Mass Step..." columns
    modecols_to_drop = [col for col in df.columns if "Mode" in col]
    masscols_to_drop = [col for col in df.columns if "Mass Step" in col]
    cols_to_drop =  modecols_to_drop + masscols_to_drop
    df = df.drop(cols_to_drop, axis=1)

    #Now, take only Step 3:
    df = df[df.Step == 3]

    #df.to_excel("___xx__FindMe.xlsx")

    #can now drop the Step and Timestamp columns:
    df = df.drop(['Step','Timestamp'],axis=1)

    #create mass list to identify columns:
    mass_list = [77,78,82,83,84,85,86,87,88,89]

    #And clean up names (convention is to label channels [Symbol][Mass]:
    names_list = df.columns.tolist()
    rename_dict = {}
    for channel in names_list:
        rename_dict[channel] = channel.replace('(','').replace(')','')

    df.rename(columns=rename_dict, inplace=True)

    metadata = {'machineName': 'Nu Plasma AttoM - PeakEScan Mode'}

    for column_name, mass in zip(df.columns.tolist(), mass_list):
        if column_name == 'time':
            continue

        data.addDataToInput(column_name, df['time'].values, df[column_name].values, metadata)
        # and set the mass for this channel:
        data.timeSeries(column_name).setProperty("Mass", mass)

    # Now calculate Total Beam:
    data.calculateTotalBeam()
    
    #create File and Sample objects for Files and Samples Browsers
    data.createFileSampleMetadata(sample_name, start_time, datetime.fromtimestamp(df['time'].iloc[-1]), importer.fileName)
    data.createImportedFileMetadata(start_time, datetime.fromtimestamp(df['time'].iloc[-1]), importer.fileName, datetime.now(), len(df.index), df.columns.tolist()[0:-1])

    importer.message('Finished')
    importer.progress(100)
    importer.finished()
