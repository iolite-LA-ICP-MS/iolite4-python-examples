# A python-based importer for iolite 4 starts with some metadata
#/ Name: NuPlasma Sr isotope importer
#/ Type: Importer
#/ Authors: Joe Petrus and Bence Paul
#/ Description: An importer for loading Nu Plasma U-Th data files (not .run files)
#/ References: None
#/ Version: 1.0
#/ Contact: support@iolite-software.com

"""
see intro.py for explanation of functions
"""

import time
import numpy as np
import pandas as pd
import re
from datetime import datetime

m_fileName = ""

def setFileName(fileName):
    m_fileName = fileName

def accepted_files():
    return ".run"

def correct_format():

    IoLog.debug("correct_format for Sr .run importer called on file = %s"%(importer.fileName))

    if importer.fileName.endswith('run'):
        IoLog.debug("Checking within .run file...")
        nrf_regex = r"Memory mapped\\Laser_Sr\.nrf"
        with open(importer.fileName) as f:
            match = re.search(nrf_regex, f.read(), re.MULTILINE)
            if match:
                return True
            else:
                IoLog.debug(".run file did not match the Sr nrf regex.")

    return False


def import_data():
    """
    see intro.py for explanation of functions
    """
    IoLog.debug("import_data in Sr importer called on file = %s"%(importer.fileName))

    #find the date and time using a regular expression:
    time_regex = r"\"(\d{2})\/(\d{2})\/(\d{4})\",\"(\d{2}):(\d{2}):(\d{2})\s([aApP][mM])\"[\n\r]"
    name_regex = r"Laser_Sr\.nrf\"\n\"(.*)\"\n"
    with open(importer.fileName) as f:
        file_contents = f.read()
        match_time = re.search(time_regex, file_contents, re.MULTILINE)
        match_name = re.search(name_regex, file_contents, re.MULTILINE)

    if not match_time:
        IoLog.error("Couldn't find a match for the date time stamp")
        importer.message('Error during import')
        importer.progress(100)
        importer.finished()
        return

    timestring = ",".join(match_time.groups())
    start_time = datetime.strptime(timestring, '%d,%m,%Y,%H,%M,%S,%p')
    start_time_in_s = start_time.timestamp()
    IoLog.debug("Start time is: " + start_time.strftime('%Y-%m-%d %H:%M:%S'))

    sample_name = ""
    if not match_name:
        IoLog.warning("Couldn't find a match for the sample name for file " + importer.fileName)
        sample_name = "Not found"
    else:
        sample_name = match_name.group(1).lstrip()
        IoLog.debug("Sample name: " + sample_name)

    '''
        Nu Plasma data .txt files contain no channel header information, so the columns that appear depends on the nrf file used to collect the data.
        You need to know what collector each column represents and give it a name below.
    '''

    names_list = ['Sr89','Sr88','Sr87','Sr86.5','Sr86','Sr85.5','Sr85','Sr84.5','Sr84','Sr83','Sr82','Sr81','Point_number','time','Zeroes_Column']
    mass_list = ['89','88','87','86.5','86','85.5','85','84.5','84','83','82','81']

    #For the time being, assume data always starts on line 58 (zero-index).
    df = pd.read_csv(importer.fileName, skiprows=58, skipfooter=3, header=None, names=names_list, engine='python')

    df = df.drop(['Point_number','Zeroes_Column'], axis=1)

    #add start time in seconds to time column:
    df['time'] = df['time'].add(start_time_in_s)

    for column, mass in zip(names_list, mass_list):
        if column is 'time':
            continue
        data.addDataToInput(column, df['time'], df[column], {"machineName": "Nu Plasma"})

    # Now calculate Total Beam:
    data.calculateTotalBeam()


    #createFileSampleMetadata(QString sampleName, PyObject fileStartTime, PyObject fileEndTime, QString filePath)
    data.createFileSampleMetadata(sample_name, start_time, datetime.fromtimestamp(df['time'].iloc[-1]), importer.fileName)
    data.createImportedFileMetadata(start_time, datetime.fromtimestamp(df['time'].iloc[-1]), importer.fileName, datetime.now(), len(df.index), names_list[0:-3])

    #importer.addImportedMassSpecFile(start_time, datetime.fromtimestamp(df['time'].iloc[-1]), importer.fileName, datetime.now(), len(df.index), names_list, "Nu Plasma")

    importer.message('Finished')
    importer.progress(100)
    importer.finished()
