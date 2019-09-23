#/ Name: Trace Elements Sum-Norm
#/ Authors: Joe Petrus and Bence Paul
#/ Description: Sum-normalizing trace elements data reduction scheme
#/ References: None
#/ Version: 1.0
#/ Contact: support@iolite-software.com

from iolite import QtGui
from iolite.Qt import Qt
import numpy as np


def runDRS():
	drs.message("Starting baseline subtract DRS...")
	drs.progress(0)

	# Get settings
	settings = drs.settings()
	print(settings)   

	indexChannel = data.timeSeries(settings["Elements"][0])

	# Create debug messages for the settings being used
	IoLog.debug("indexChannelName = %s" % indexChannel.name)

	# Setup index time
	drs.message("Setting up index time...")
	drs.progress(5)
	drs.setIndexChannel(indexChannel)

	# Interp onto index time and baseline subtract
	drs.message("Interpolating onto index time and baseline subtracting...")
	drs.progress(25)

	allInputChannels = data.timeSeriesList(data.Input)
	
	commonProps = {'DRS': drs.name()}

	# Baseline subtract and calculate SQ concentrations for each channel
	for counter, channel in enumerate(allInputChannels):
		drs.message("Baseline subtracting %s" % channel.name)
		drs.progress(25 + 25*counter/len(allInputChannels))

		drs.baselineSubtract(data.selectionGroup("Baseline"), [allInputChannels[counter]], None, 25, 50)

		try:
			rm = data.referenceMaterialData(settings["External"])
			channel_ppm_in_rm = rm[channel.property('Element')].valueInPPM()
			channel_spline = data.spline(settings['External'], channel.name).data()
			channel_cps = data.timeSeries(channel.name + "_CPS").data()
			channel_ppm = channel_cps * channel_ppm_in_rm / channel_spline
			data.createTimeSeries(channel.name + " ppm", data.Output, None, channel_ppm, commonProps)
		except KeyError:
			print('Could not calculate SQ channel for %s'%(channel.name))

	# Work out normalizing factor
	channels_for_norm = [c + ' ppm' for c in settings['Elements']]
	channels_sum = np.zeros(len(indexChannel.data()))
	for channel in channels_for_norm:
		# Need to add in oxide ability!!
		channels_sum += data.timeSeries(channel).data()
	
	factor = (channels_sum/1e4)/settings['Value']
	data.createTimeSeries("NormalizationFactor", data.Intermediate, None, factor, commonProps)
	
	channels_to_adjust = [c for c in data.timeSeriesList(data.Output) if 'ppm' in c.name]
	print(channels_to_adjust)
	for channel in channels_to_adjust:
		print('Adjusting %s'%(channel.name))
		d = channel.data()/factor
		channel.setData(d)

			

	drs.message("Finished!")
	drs.progress(100)
	drs.finished()
	

def settingsWidget():
	widget = QtGui.QWidget()
	formLayout = QtGui.QFormLayout()
	formLayout.setFieldGrowthPolicy(QtGui.QFormLayout.FieldsStayAtSizeHint)
	formLayout.setFormAlignment(Qt.AlignHCenter | Qt.AlignTop)
	widget.setLayout(formLayout)

	rmComboBox = QtGui.QComboBox(widget)	
	rmComboBox.addItems(data.referenceMaterialNames())
	rmComboBox.currentTextChanged.connect(lambda s: drs.setSetting("External", str(s)))
	formLayout.addRow("External reference material", rmComboBox)

	elementsList = QtGui.QListWidget(widget)
	elementsList.addItems([c for c in data.timeSeriesNames(data.Input) if 'TotalBeam' not in c])
	elementsList.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
	def updateElements():
		elements = [wi.text() for wi in elementsList.selectedItems()]
		drs.setSetting("Elements", elements)
	elementsList.itemSelectionChanged.connect(updateElements)
	formLayout.addRow("Elements to normalize", elementsList)

	oxideCheckBox = QtGui.QCheckBox(widget)
	oxideCheckBox.toggled.connect(lambda b: drs.setSetting("Oxides", bool(b)))
	formLayout.addRow("Oxides?", oxideCheckBox)

	valueLineEdit = QtGui.QLineEdit(widget)
	valueLineEdit.textChanged.connect(lambda v: drs.setSetting("Value", float(v)))
	formLayout.addRow("Value (wt. %)", valueLineEdit)
	
	# Restore settings
	try:
		settings = drs.settings()
		rmComboBox.setCurrentText(settings["External"])
		elementsList.selectionModel().clear()
		for el in settings["Elements"]:
			matches = elementsList.findItems(el, Qt.MatchFixedString)
			if matches:
				matches[0].setSelected(True)
		oxideCheckBox.setChecked(settings["Oxides"])
		valueLineEdit.setText(str(settings["Value"]))
	except KeyError:
		pass

	drs.setSettingsWidget(widget)
