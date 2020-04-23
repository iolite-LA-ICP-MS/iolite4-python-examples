# A python-based importer for iolite 4 starts with some metadata
#/ Type: Importer
#/ Name: IsotopX TIMSDP file importer
#/ Authors: Bence Paul and Joe Petrus
#/ Description: An importer for loading IsotopX TIMSDP files
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
    return ".TIMSDP"

def correct_format():

    if importer.fileName.endswith('TIMSDP'):
        return True

    return False

def import_data():
    """
    see intro.py for explanation of functions
    """
    print("\n\n")

    IoLog.debug("import_data in IsotopX TIMSDP importer called on file = %s"%(importer.fileName))

    #find the date and time using a regular expression:
    time_regex = r"AnalysisStart,(\d{1,2}\/\d{1,2}\/\d{4}\s+\d{2}:\d{2}:\d{2})"
    name_regex = "SampleID, (.*)"
    with open(importer.fileName) as f:
        file_contents = f.read()
        match_time = re.search(time_regex, file_contents, re.MULTILINE)
        match_name = re.search(name_regex, file_contents, re.MULTILINE)

    file_contents = file_contents.split('\n')

    if not match_time:
        IoLog.error("Couldn't find a match for the date time stamp")
        importer.message('Error during import')
        importer.progress(100)
        importer.finished()
        return

    timestring = match_time.group(1)
    sample_name = match_name.group(1)
    if len(sample_name) < 1:
        sample_name = "[blank]"

    start_time = datetime.strptime(timestring, '%d/%m/%Y %H:%M:%S')
    start_time_in_s = start_time.timestamp()
    IoLog.debug("Start time is: " + start_time.strftime('%Y-%m-%d %H:%M:%S'))

    # Get all metadata as dictionary:
    metadata = {}

    cycles_line = 0
    line_counter = 0
    blocks_line = 0
    read_metadata = True        # will be set to false once we get to the #USERTABLES line

    for line in file_contents:
        #TODO: add collector metadata
        
        if re.match(r"#CYCLES", line):
            cycles_line = line_counter
        
        if re.match(r"#BLOCKS", line):
            blocks_line = line_counter
            break       
        
        if re.match(r"#USERTABLES", line):
            read_metadata = False
        
        # This will read metadata until we get to the USERTABLES line:
        if len(line.split(',')) > 1 and read_metadata:
            metadata[line.split(',')[0]] = line.split(',')[1]
        
        line_counter += 1

    names_list = file_contents[cycles_line + 1].split(',')
    footer_lines = len(file_contents) - blocks_line
    
    df = pd.read_csv(importer.fileName, skiprows=cycles_line + 2, skipfooter=footer_lines, header=None, names=names_list, engine='python')

    #add start time in seconds to time column:
    df['Time'] = df['Time'].add(start_time_in_s)

    # Add blank row at end of dataframe here to create a blank between each file
    df = df.append(pd.Series(), ignore_index=True)

    # And add 1 millisecond to last time for the timing of this blank value
    df['Time'].iloc[-1] = df['Time'].iloc[-2] + 0.001

    # Would drop unnecessary columns here (e.g. those that are all zeros etc)
    #e.g. df = df.drop(['Cycle'], axis=1)

    # clean up channel list to give proper names etc:
    channels_list = []
    for name in names_list:
        if name in ['Cycle', 'Time']:
            continue
        
        # add 'm' to start of channel with just mass in thier name:
        if re.match(r'^\d+$', name):
            channels_list.append('m'+name)
        else:
            channels_list.append(name)

    # now add data as channels:
    for channel in names_list[2:]:
        if re.match(r'^\d+$', channel):
            data.addDataToInput('m'+channel, df['Time'], df[channel], {"machineName": "IsotopX TIMS", "mass": channel})
        else:
            data.addDataToInput(channel, df['Time'], df[channel], {"machineName": "IsotopX TIMS"})

    # Now calculate Total Beam:
    data.calculateTotalBeam()

    data.createFileSampleMetadata(sample_name, start_time, datetime.fromtimestamp(df['Time'].iloc[-1]), importer.fileName)
    f = data.createImportedFileMetadata(start_time, datetime.fromtimestamp(df['Time'].iloc[-1]), importer.fileName, datetime.now(), len(df.index), channels_list)
    f.setProperty("metadata", metadata)

    importer.message('Finished')
    importer.progress(100)
    importer.finished()
