# A python-based importer for iolite 4 starts with some metadata
#/ Name: NuPlasma .run importer
#/ Authors: Joe Petrus and Bence Paul
#/ Description: An importer for loading Nu Plasma .run files
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
    IoLog.debug("correct_format called on file = %s"%(importer.fileName))

    if importer.fileName.endswith('run'):
        return True

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
    IoLog.debug("import_data called on file = %s"%(importer.fileName))

    #find the date and time using a regular expression:
    regex = r"^\"(?P<day>\d{2})\/(?P<month>\d{2})\/(?P<year>\d{4})\",\"(?P<hour>\d{2}):(?P<minute>\d{2}):(?P<seconds>\d{2}) (?P<ampm>[AaPpMm]{2})\"$"

    datafile = ""
    with open(importer.fileName) as f:
        match = re.search(regex, f.read(), re.MULTILINE)

    if not match:
        IoLog.error("Couldn't find a match for the date time stamp")
        importer.message('Error during import')
        importer.progress(100)
        importer.finished()

    timestring = ",".join(match.groups())
    start_time = datetime.strptime(timestring, '%d,%m,%Y,%I,%M,%S,%p')

    #IoLog.information("Datetime string AFTER CONVERSION is " + start_time.strftime('%Y-%m-%d %H:%M:%S'))

    IoLog.information("TIMESTAMP is " + str(start_time.timestamp()))

    #look through file to find the first line of data, which in .run files
    # is preceded by a line with "Spare Text" written in it.
    first_line_of_data = 0
    with open(importer.fileName, 'r') as file:
        for line in file:
            if "Spare text" in line:
                break
            first_line_of_data += 1

    first_line_of_data += 1
    IoLog.debug("Found spare text on line: {}".format(first_line_of_data))


    '''
        Nu Plasma .run files contain no header data, so the columns that appear depends on the nrf file used to collect the data.
        You need to know what collector each column represents and give it a name below.
        This particular example is for loading Sr data (in particular, using our Laser_Sr.nrf format). We could write a check in our
        correct_format() function above to check that the nrf used is Laser_Sr.nrf
    '''

    names_list = ['Sr88','Sr87','Sr86','Rb85','Sr84','CaAr83','CaAr82','m81','m80','m79','m78','m77','measurement','time','type_col']

    df = pd.read_csv(importer.fileName, skiprows=first_line_of_data, header=None, names=names_list)

    IoLog.information(str(df))



    # Normally you would parse the file specified by plugin.fileName
    # into some numpy arrays and then pass those to iolite. Rather
    # than parse some data, for this example, we'll just generate
    # some random data.

    # num_points = 1000
    # t = np.linspace(time.time(), time.time()+60*60, num=num_points)
    #
    # for i in range(10):
    #     # Update our task information
    #     importer.message('Now doing %i'%i)
    #     importer.progress(i*10)
    #
    #     # Sleep a bit so it doesn't happen too fast!
    #     # For demonstration purposes only...
    #     time.sleep(1)
    #
    #     # Make some random data
    #     d = np.random.randn(num_points)
    #
    #     # Note that the channel type is specified through an
    #     # enum, here data.Input
    #     data.createTimeSeries('Channel%i'%i, data.Input, t, d)

    importer.message('Finished')
    importer.progress(100)
    importer.finished()
