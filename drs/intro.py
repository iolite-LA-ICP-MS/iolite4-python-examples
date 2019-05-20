# A python-based data reduction scheme for iolite 4 starts with some metadata
#/ Name: Data Reduction Scheme Introduction
#/ Authors: Joe Petrus and Bence Paul
#/ Description: This is an example python-based plugin for processing data in iolite 4
#/ References: None
#/ Version: 1.0
#/ Contact: support@iolite-software.com

"""
Qt imports can be done through 'iolite', e.g.
from iolite.QtGui import QLabel

Before functions are called a few additional objects are
added to the module:

data	an interface to iolite's C++ data. E.g. you can get
        existing time series data or selection groups with it
        as well as make new ones.

IoLog	an interface to iolite's logging facility. You can add
        messages with, e.g., IoLog.debug('My message')

drs		an interface to the PythonDRS C++ class in iolite from 
		which some built-in features can be accessed, e.g.,
		baselineSubtract(group, channels, mask)
"""

from iolite import QtGui
from iolite.TimeSeriesData import TimeSeriesData
from time import sleep
import numpy as np


def runDRS():
	"""
	This method will be called by iolite when the user clicks 
	Crunch Data in the DRS window or as part of a processing
	template. It should transform 'input' data into 'output'
	data using the provided settings.

	DRS progress can be updated via the 'message' and 'progress'
	signals. These will be displayed in the iolite interface.

	When finished, the 'finished' signal should be emitted.

	As an example, we will do baseline subtraction of all
	input channels using a DRS helper function.

	"""

	drs.message("Starting baseline subtract DRS...")
	drs.progress(0)

	# Get settings
	settings = drs.settings()
	print(settings)   

	indexChannel = data.timeSeries(settings["IndexChannel"])
	maskChannel = data.timeSeries(settings["MaskChannel"])
	cutoff = settings["MaskCutoff"]
	trim = settings["MaskTrim"]

	# Create debug messages for the settings being used
	IoLog.debug("indexChannelName = %s" % indexChannel.name())
	IoLog.debug("maskChannelName = %s" % maskChannel.name())
	IoLog.debug("maskCutoff = %f" % cutoff)
	IoLog.debug("maskTrim = %f" % trim)

	# Setup index time
	drs.message("Setting up index time...")
	drs.progress(5)
	drs.setIndexChannel(indexChannel)

	# Setup the mask
	drs.message("Making mask...")
	drs.progress(10)
	mask = drs.createMaskFromCutoff(maskChannel, cutoff, trim)

	# Interp onto index time and baseline subtract
	drs.message("Interpolating onto index time and baseline subtracting...")
	drs.progress(25)

	allInputChannels = data.timeSeriesList(TimeSeriesData.tsInput)

	for counter, channel in enumerate(allInputChannels):
		drs.message("Baseline subtracting %s" % channel.name())
		drs.progress(25 + 75*counter/len(allInputChannels))
		sleep(5) # Sleeping only so that the progress can be observed

		drs.baselineSubtract(data.selectionGroup("Baseline"), [allInputChannels[counter]], mask, 25, 100)

	drs.message("Finished!")
	drs.progress(100)
	drs.finished()
	

def settingsWidget():
	"""
	This function puts together a user interface to configure the DRS.
	
	It is important to have the last line of this function call:
	drs.setSettingsWidget(widget)
	"""

	widget = QtGui.QWidget()
	formLayout = QtGui.QFormLayout()
	widget.setLayout(formLayout)

	timeSeriesNames = data.timeSeriesNames(TimeSeriesData.tsInput)
	defaultChannelName = ""
	if timeSeriesNames:
		defaultChannelName = timeSeriesNames[0]

	drs.setDefaultSetting("IndexChannel", defaultChannelName)
	drs.setDefaultSetting("MaskChannel", defaultChannelName)
	drs.setDefaultSetting("MaskCutoff", 100000.0)
	drs.setDefaultSetting("MaskTrim", 0.0)

	settings = drs.settings()

	indexComboBox = QtGui.QComboBox(widget)
	indexComboBox.addItems(data.timeSeriesNames(TimeSeriesData.tsInput))
	indexComboBox.setCurrentText(settings["IndexChannel"])
	indexComboBox.currentTextChanged.connect(lambda t: drs.setSetting("IndexChannel", t))

	maskComboBox = QtGui.QComboBox(widget)
	maskComboBox.addItems(data.timeSeriesNames(TimeSeriesData.tsInput))
	maskComboBox.setCurrentText(settings["MaskChannel"])
	maskComboBox.currentTextChanged.connect(lambda t: drs.setSetting("MaskChannel", t))

	maskLineEdit = QtGui.QLineEdit(widget)
	maskLineEdit.setText(settings["MaskCutoff"])
	maskLineEdit.textChanged.connect(lambda t: drs.setSetting("MaskCutoff", float(t)))

	maskTrimLineEdit = QtGui.QLineEdit(widget)
	maskTrimLineEdit.setText(settings["MaskTrim"])
	maskTrimLineEdit.textChanged.connect(lambda t: drs.setSetting("MaskTrim", float(t)))

	formLayout.addRow("Index channel", indexComboBox)
	formLayout.addRow("Mask channel", maskComboBox)
	formLayout.addRow("Mask cutoff", maskLineEdit)
	formLayout.addRow("Mask trim", maskTrimLineEdit)

	drs.setSettingsWidget(widget)
