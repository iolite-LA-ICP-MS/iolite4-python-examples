# A python-based importer for iolite 4 starts with some metadata
#/ Name: NuPlasma U-Th importer
#/ Authors: Joe Petrus and Bence Paul
#/ Description: An importer for loading Nu Plasma U-Th data files (not .run files)
#/ References: None
#/ Version: 1.0
#/ Contact: support@iolite-software.com

import time
import numpy as np
import pandas as pd
import re
from datetime import datetime

constants = {
 '229Th': {'Decay constant': 9.44342e-05, 'Isotope Atomic mass': 229.0317553},
 '230Th': {'Decay constant': 9.15771e-06, 'Isotope Atomic mass': 230.0331266},
 '231Pa': {'Decay constant': 2.11583e-05, 'Isotope Atomic mass': 231.0358789},
 '232Th': {'Decay constant': 4.93343e-11, 'Isotope Atomic mass': 232.0380504},
 '233U': {'Decay constant': 4.35394e-06, 'Isotope Atomic mass': 233.0396282},
 '234U': {'Decay constant': 2.82629e-06, 'Isotope Atomic mass': 234.0409456},
 '235U': {'Decay constant': 9.84864e-10, 'Isotope Atomic mass': 235.0439231},
 '236U': {'Decay constant': 2.95964e-08, 'Isotope Atomic mass': 236.0455619},
 '238U': {'Decay constant': 1.55125e-10, 'Isotope Atomic mass': 238.0507826},
 '238/235 nat': 137.88,
 '229/233 spk': 0.3607895,
 '230/229 spk': 0.00005,
 '234/233 spk': 0.0003524,
 '236/233 spk': 14.463645,
 'assumed IC2 gain': 0.92,
 'assumed Raw 238/235': 140.5,
 'CPS/Volts': 62500000.0,
 '(ThH/Th)/(UH/U)': 0.0,
 'NP080 IC0 DT': 19.9,
 'NP080 IC1 DT': 20.0
 }



m_fileName = ""

def setFileName(fileName):
    m_fileName = fileName

def accepted_files():
    return ".txt"

def correct_format():
    """
    This method will be called by iolite when the user selects a file
    to import. Typically, it uses the provided name (stored in
    plugin.fileName) and parses as much of it as necessary to determine
    if this importer is appropriate to import the data. For example,
    although X Series II and Agilent data are both comma separated
    value files, they can be distinguished by the characteristic
    formatting in each. In our implementation, distinguishing the
    two is done with 'regular expressions' (QRegularExpression)
    parsing of the first several lines of the file.

    Keep in mind that there is nothing stopping you from just
    returning True (thus supporting all files!) or simply checking
    the file extension, but such generic checks can yield unexpected
    results. You cannot be sure which order the various importer
    plugins will be checked for compatibility.

    This method must return either True or False.
    """
    IoLog.debug("correct_format for U-Th importer called on file = %s"%(importer.fileName))

    if importer.fileName.endswith('txt'):
        IoLog.debug("Checking within .run file...")
        nrf_regex = r"UThAge_II\.nrf"
        with open(importer.fileName) as f:
            match = re.search(nrf_regex, f.read(), re.MULTILINE)
            if match:
                return True
            else:
                IoLog.debug(".run file did not match the U-Th nrf regex.")

    return False

def import_data():
    """
    This method uses the provided file name (stored in plugin.fileName),
    parses its contents, and registers time series data with iolite by
    emitting the timeSeriesData signal.

    Note that emitting signals here simply means calling the corresponding
    function, e.g. importer.message('My message')

    Importer progress can be updated via the 'message' and 'progress'
    signals. These will be displayed in the iolite interface.

    When finished, the 'finished' signal should be emitted.
    """
    IoLog.debug("import_data in U-Th importer called on file = %s"%(importer.fileName))

    #find the date and time using a regular expression:
    time_regex = r"\"Started analysis at (?P<hour>\d+):(?P<minutes>\d+) [\w ]+, (?P<month>\w+) (?P<day>\d+), (?P<year>\d+)\""
    name_regex = r"\"Sample Name is (.*)\""
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
    start_time = datetime.strptime(timestring, '%H,%M,%B,%d,%Y')
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

    names_list = ['H5 base 1','H4 base 1','H3 base 1','H2 base 1','H1 base 1','Ax base 1','L1 base 1','L2 base 1','L3 base 1','IC0 base 1','L4 base 1','L5 base 1','L6 base 1',
                    'IC1 base 1','L7 base 1','H5 base 2','H4 base 2','H3 base 2','H2 base 2','H1 base 2','Ax base 2','L1 base 2','L2 base 2','L3 base 2','IC0 base 2','L4 base 2',
                    'L5 base 2','L6 base 2','IC1 base 2','L7 base 2','H5 Cycle 1','H4 Cycle 1','H3 Cycle 1','H2 Cycle 1','H1 Cycle 1','Ax Cycle 1','L1 Cycle 1','L2 Cycle 1',
                    'L3 Cycle 1','IC0 Cycle 1','L4 Cycle 1','L5 Cycle 1','L6 Cycle 1','IC1 Cycle 1','L7 Cycle 1','H5 Cycle 2','H4 Cycle 2','H3 Cycle 2','H2 Cycle 2','H1 Cycle 2',
                    'Ax Cycle 2','L1 Cycle 2','L2 Cycle 2','L3 Cycle 2','IC0 Cycle 2','L4 Cycle 2','L5 Cycle 2','L6 Cycle 2','IC1 Cycle 2','L7 Cycle 2', 'measurement']

    #For the time being, assume data always starts on line 21 (zero-index). NOTE: had to use skiprows=40... not sure why, perhaps new line characters in file? Anyway, this seems to work for test file...
    df = pd.read_csv(importer.fileName, skiprows=41, skipfooter=3, header=None, names=names_list)

    #Testing accessing constants:
    #IoLog.debug("Here's the constants: " + str(constants['234U']['Decay constant']) + "\n" + str(constants['NP080 IC0 DT']))

    #df.to_excel('/Users/bence/iolite4_xyz/iolite4/tests/UTh/test_import.xlsx')

    return

    #add start time in seconds to time column:
    df['time'] = df['time'].add(start_time_in_s)

    for column, mass in zip(names_list, mass_list):
        if column in ['measurement','time','type_col']:
            continue
        channel = data.createTimeSeries(column, data.Input, df['time'].values, df[column].values)
        channel.setProperty("Mass", mass)
        channel.setProperty("Units", "volts")


    # Now calculate Total Beam:
    data.calculateTotalBeam()

    # IoLog.debug("Start Time: " + start_time.strftime("%d,%m,%Y,%I,%M,%S,%p"))
    # IoLog.debug("End Time: " + datetime.fromtimestamp(df['time'].iloc[-1]).strftime("%d,%m,%Y,%I,%M,%S,%p"))
    # IoLog.debug("Import Time: " + datetime.now().strftime("%d,%m,%Y,%I,%M,%S,%p"))
    # IoLog.debug("No Of Points: {}".format(len(df.index)))
    # IoLog.debug("Mass list: " + ', '.join(names_list))
    #
    # importer.addImportedMassSpecFile(start_time, datetime.fromtimestamp(df['time'].iloc[-1]), importer.fileName, datetime.now(), len(df.index), ', '.join(names_list), "Nu Plasma")

    importer.message('Finished')
    importer.progress(100)
    importer.finished()
