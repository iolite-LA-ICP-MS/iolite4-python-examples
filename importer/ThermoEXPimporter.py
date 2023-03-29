# A python-based importer for iolite 4 starts with some metadata
#/ Type: Importer
#/ Name: Thermo .exp file importer
#/ Authors: Bence Paul and Joe Petrus
#/ Description: An importer for loading .exp files from Neptunes and Tritons
#/ References: None
#/ Version: 1.0
#/ Contact: support@iolite-software.com

import numpy as np
import pandas as pd
import re
from datetime import datetime


m_fileName = ""

def setFileName(fileName):
    m_fileName = fileName

def accepted_files():
    return ".exp"

def correct_format():

    if importer.fileName.endswith('exp'):
        return True

    return False

def import_data():
    """
    see intro.py for explanation of functions
    """
    print("\n\n")

    IoLog.debug("import_data in Neptune/Triton .exp importer called on file = %s"%(importer.fileName))

    # Not going to convert the date/timestamp here like normal.
    # Going to do it for each line below
    name_regex = r"Sample ID: (.*)"
    with open(importer.fileName) as f:
        file_contents = f.read()
        match_name = re.search(name_regex, file_contents, re.MULTILINE)

    sample_name = match_name.group(1).strip()
    if len(sample_name) < 1:
        sample_name = "[blank]"

    file_contents = file_contents.split('\n')

    # ToDo: collect more metadata in file?
    metadata = {}

    line_counter = 0
    headers_line = 0
    read_metadata = True        # will be set to false once we get to the #USERTABLES line

    for line in file_contents:
        if line.startswith('Triton'):
            DATE_FORMAT = '%Y/%m/%d'
            machine_name = 'Triton'

        if line.startswith('Neptune'):
            DATE_FORMAT = '%d/%m/%Y'
            machine_name = 'Neptune'

        if line.startswith('Analysis date:'):
            datestring = line.split(':')[1]
            datestring = datestring.strip('\t ')

        if line.startswith('Cycle\tTime'):
            header_line = line_counter

        #if we find the '*** Cup' line, we know we've reached the end of the data...
        if line.startswith('***\t'):
            footer_lines = len(file_contents) - line_counter - 1
            break

        line_counter += 1

    names_list = file_contents[header_line].split('\t')

    df = pd.read_csv(importer.fileName, sep='\t', skiprows=header_line + 1, skipfooter=footer_lines, header=None, names=names_list, engine='python')

    def convert_time(time_str, date_str, date_format):
        secs = datetime.strptime(date_str + time_str + '000', date_format + '%H:%M:%S:%f').timestamp()
        return secs

    # Convert the time in the Time column to seconds since epoch for iolite
    df['Time'] = df['Time'].apply(convert_time, args=(datestring, DATE_FORMAT))

    # Add blank row at end of dataframe here to create a blank between each file
    df = df.append(pd.Series(), ignore_index=True)

    # And add 1 millisecond to last time for the timing of this blank value
    df['Time'].iloc[-1] = df['Time'].iloc[-2] + 0.001
    # Define start and end time here
    start_time = datetime.fromtimestamp(df['Time'].iloc[0])
    end_time = datetime.fromtimestamp(df['Time'].iloc[-1])

    # clean up channel list to give proper names etc:
    channels_list = []
    for name in names_list:
        if name in ['Cycle', 'Time']:
            continue

        # Don't include blank column, which is last column in df
        if len(name) < 1:
            continue

        # add 'm' to start of channel with just mass in their name:
        if re.match(r'^[\d.]+$', name):
            channels_list.append('m'+name)
        else:
            channels_list.append(name)

    # now add data as channels:
    for channel in names_list[2:]:
        if len(channel) < 1:
            continue # There is a blank column as the last column

        elif re.match(r'^[\d.]+\w+[+]{2}$', channel): # Handle ++ channels, such as 173Yb++
            mass = float(re.match(r'^([\d.]+)', channel).group(1)) / 2.
            data.addDataToInput(channel, df['Time'], df[channel], {"machineName": machine_name})
            ch = data.timeSeries(channel)
            ch.setProperty('Mass', mass)


        elif re.match(r'^[\d.]+$', channel):  # Handle channels such as '83.5'
            data.addDataToInput('m'+channel, df['Time'], df[channel], {"machineName": machine_name, "mass": channel})

        else:
            data.addDataToInput(channel, df['Time'], df[channel], {"machineName": machine_name})

    # Now calculate Total Beam:
    data.calculateTotalBeam()

    data.createFileSampleMetadata(sample_name, start_time, end_time, importer.fileName)
    f = data.createImportedFileMetadata(start_time, end_time, importer.fileName, datetime.now(), len(df.index), channels_list)
    f.setProperty("metadata", metadata)

    importer.message('Finished')
    importer.progress(100)
