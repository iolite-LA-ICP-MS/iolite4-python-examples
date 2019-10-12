# A python-based importer for iolite 4 starts with some metadata
#/ Type: Importer
#/ Name: Importer Introduction
#/ Authors: Joe Petrus and Bence Paul
#/ Description: This is an example python-based plugin for importing data into iolite 4
#/ References: None
#/ Version: 1.0
#/ Contact: support@iolite-software.com

import time
import numpy as np

"""
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

    if importer.fileName.endswith('ioe'):
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


    # Normally you would parse the file specified by plugin.fileName
    # into some numpy arrays and then pass those to iolite. Rather
    # than parse some data, for this example, we'll just generate
    # some random data.

    num_points = 1000
    t = np.linspace(time.time(), time.time()+60*60, num=num_points)

    for i in range(10):
        # Update our task information
        importer.message('Now doing %i'%i) 
        importer.progress(i*10)

        # Sleep a bit so it doesn't happen too fast!
        # For demonstration purposes only...
        time.sleep(1)

        # Make some random data
        d = np.random.randn(num_points)

        # Note that the channel type is specified through an
        # enum, here data.Input
        data.createTimeSeries('Channel%i'%i, data.Input, t, d)

    importer.message('Finished')
    importer.progress(100)
    importer.finished()
