# A python-based data reduction scheme for iolite 4 starts with some metadata:
#/ Name: Hf isotopes example
#/ Authors: Joe Petrus and Bence Paul
#/ Description: This is an example python-based plugin for processing Hf data in iolite 4
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
from iolite.SelectionGroup import SelectionGroup
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

	drs.message("Starting Hf isotopes DRS...")
	drs.progress(0)

	# Get settings
	settings = drs.settings()
	print(settings)   

	indexChannel = data.timeSeries(settings["IndexChannel"])
	maskChannel = data.timeSeries(settings["MaskChannel"])
	rmName = settings["ReferenceMaterial"]
	cutoff = settings["MaskCutoff"]
	trim = settings["MaskTrim"]
	HfTrue = settings["HfTrue"]
	Yb31 = settings["Yb31"]
	Yb63 = settings["Yb63"]
	age = settings["Age"]
	propErrors = settings["PropagateError"]

	# Create debug messages for the settings being used
	IoLog.debug("indexChannelName = %s" % indexChannel.name)
	IoLog.debug("maskChannelName = %s" % maskChannel.name)
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
		drs.message("Baseline subtracting %s" % channel.name)
		drs.progress(25 + 50*counter/len(allInputChannels))

		drs.baselineSubtract(data.selectionGroup("Baseline_1"), [allInputChannels[counter]], mask, 25, 75)


	HfLuYb176 = data.timeSeries("Hf176_CPS").data()
	Hf177 = data.timeSeries("Hf177_CPS").data()
	Hf178 = data.timeSeries("Hf178_CPS").data()
	Hf179 = data.timeSeries("Hf179_CPS").data()
	Yb171 = data.timeSeries("Yb171_CPS").data()
	Yb173 = data.timeSeries("Yb173_CPS").data()
	Lu175 = data.timeSeries("Lu175_CPS").data()


	HfFract = (np.log(HfTrue / (Hf179/Hf177))) / (np.log(178.946 / 176.943))*mask
	YbFract = (np.log(Yb31 /(Yb173/Yb171))) / (np.log(172.938222 / 170.936338))*mask
	Yb176 = Yb173 * (Yb63 / (np.power((175.942576 / 172.938222) , YbFract)))
	Lu176 = Lu175 * (0.02656 / (np.power((175.942694 / 174.9408) , YbFract)))
	Hf176c = HfLuYb176 - Yb176 - Lu176
	LuHf176 = HfLuYb176 - Yb176
	YbHf176 = HfLuYb176 - Lu176
	Yb_PPM_on_176 = (Yb176 / HfLuYb176) * 1000000*mask
	Lu_PPM_on_176 = (Lu176 / HfLuYb176) * 1000000*mask
	Hf176_177_Raw = (HfLuYb176 / Hf177) * np.power((175.941 / 176.943) , HfFract)*mask
	Hf176_177_Corr = (Hf176c/ Hf177) * np.power((175.941 / 176.943) , HfFract)*mask
	Hf176_177_LuCorr = (YbHf176 / Hf177) * np.power((175.941 / 176.943) , HfFract)*mask
	Hf176_177_YbCorr = (LuHf176 / Hf177) * np.power((175.941 / 176.943) , HfFract)*mask
	Hf178_177 = (Hf178 / Hf177) * np.power((177.944 / 176.943) , HfFract)*mask
	Lu176_Hf177_Raw = (Lu176 / Hf177)*mask
	Lu176_Hf177_Corr = (Lu176 / Hf177) * np.power((175.942694 / 176.943), (0.5*(HfFract + YbFract)))*mask
	Yb176_Hf177_Raw = (Yb176 / Hf177)*mask
	Yb176_Hf177_Corr = (Yb176 / Hf177) * np.power((175.942576 / 176.943), (0.5*(HfFract + YbFract)))*mask
	TotalHfBeam =  Hf178 / 0.27297

	data.createTimeSeries("Hf176_177_Corr", data.Output, indexChannel.time(), Hf176_177_Corr)
	data.createTimeSeries("Hf178_177", data.Output, indexChannel.time(), Hf178_177)
	
	StdSpline_Hf176_177 = data.spline(rmName, "Hf176_177_Corr").data()
	StdValue_Hf176_177 = data.referenceMaterialData(rmName)["176Hf/177Hf"].value()

	print("StdSpline_Hf176_177 mean = %f"%StdSpline_Hf176_177.mean())
	print("StdValue_Hf176_177 = %f"%StdValue_Hf176_177)

	StdCorr_Hf176_177 = (Hf176_177_Corr)* StdValue_Hf176_177 / StdSpline_Hf176_177

	data.createTimeSeries("StdCorr_Hf176_177", data.Output, indexChannel.time(), StdCorr_Hf176_177)

	#StdSpline_Hf178_177 = data.spline(rmName, "Hf178_177").data()
	#StdValue_Hf178_177 = data.referenceMaterialData(rmName)["178Hf/177Hf"].value()
	#StdCorr_Hf178_177= (Hf178_177)* StdValue_Hf178_177 / StdSpline_Hf178_177

	if propErrors:
		groups = [s for s in data.selectionGroupList() if s.type != SelectionGroup.sgBaseline]
		data.propagateErrors(groups, [data.timeSeries("StdCorr_Hf176_177")], data.timeSeries("Hf176_177_Corr"), rmName)


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

	rmNames = data.selectionGroupNames(SelectionGroup.sgReferenceMaterial)

	drs.setDefaultSetting("IndexChannel", defaultChannelName)
	drs.setDefaultSetting("ReferenceMaterial", "Z_Plesovice")
	drs.setDefaultSetting("MaskChannel", defaultChannelName)
	drs.setDefaultSetting("MaskCutoff", 0.1)
	drs.setDefaultSetting("MaskTrim", 0.0)
	drs.setDefaultSetting("HfTrue", 0.7325)
	drs.setDefaultSetting("Yb31", 1.132685)
	drs.setDefaultSetting("Yb63", 0.796218)
	drs.setDefaultSetting("Age", 0)
	drs.setDefaultSetting("PropagateError", False)

	settings = drs.settings()

	indexComboBox = QtGui.QComboBox(widget)
	indexComboBox.addItems(timeSeriesNames)
	indexComboBox.setCurrentText(settings["IndexChannel"])
	indexComboBox.currentTextChanged.connect(lambda t: drs.setSetting("IndexChannel", t))
	formLayout.addRow("Index channel", indexComboBox)

	rmComboBox = QtGui.QComboBox(widget)
	rmComboBox.addItems(rmNames)
	rmComboBox.setCurrentText(settings["ReferenceMaterial"])
	rmComboBox.currentTextChanged.connect(lambda t: drs.setSetting("ReferenceMaterial", t))
	formLayout.addRow("Reference material", rmComboBox)    

	maskComboBox = QtGui.QComboBox(widget)
	maskComboBox.addItems(data.timeSeriesNames(TimeSeriesData.tsInput))
	maskComboBox.setCurrentText(settings["MaskChannel"])
	maskComboBox.currentTextChanged.connect(lambda t: drs.setSetting("MaskChannel", t))
	formLayout.addRow("Mask channel", maskComboBox)

	maskLineEdit = QtGui.QLineEdit(widget)
	maskLineEdit.setText(settings["MaskCutoff"])
	maskLineEdit.textChanged.connect(lambda t: drs.setSetting("MaskCutoff", float(t)))
	formLayout.addRow("Mask cutoff", maskLineEdit)

	maskTrimLineEdit = QtGui.QLineEdit(widget)
	maskTrimLineEdit.setText(settings["MaskTrim"])
	maskTrimLineEdit.textChanged.connect(lambda t: drs.setSetting("MaskTrim", float(t)))
	formLayout.addRow("Mask trim", maskTrimLineEdit)

	hfTrueLineEdit = QtGui.QLineEdit(widget)
	hfTrueLineEdit.setText(settings["HfTrue"])
	hfTrueLineEdit.textChanged.connect(lambda t: drs.setSetting("HfTrue", float(t)))
	formLayout.addRow("Hf true", hfTrueLineEdit)

	Yb31LineEdit = QtGui.QLineEdit(widget)
	Yb31LineEdit.setText(settings["Yb31"])
	Yb31LineEdit.textChanged.connect(lambda t: drs.setSetting("Yb31", float(t)))
	formLayout.addRow("173Yb/171Yb", Yb31LineEdit)

	Yb63LineEdit = QtGui.QLineEdit(widget)
	Yb63LineEdit.setText(settings["Yb63"])
	Yb63LineEdit.textChanged.connect(lambda t: drs.setSetting("Yb63", float(t)))
	formLayout.addRow("176Yb/173Yb", Yb63LineEdit)

	ageLineEdit = QtGui.QLineEdit(widget)
	ageLineEdit.setText(settings["Age"])
	ageLineEdit.textChanged.connect(lambda t: drs.setSetting("Age", float(t)))
	formLayout.addRow("Age", ageLineEdit)

	propCheckBox = QtGui.QCheckBox(widget)
	propCheckBox.setChecked(settings["PropagateError"])
	propCheckBox.toggled.connect(lambda t: drs.setSetting("PropagateError", bool(t)))
	formLayout.addRow("PropagateError", propCheckBox)	

	drs.setSettingsWidget(widget)
