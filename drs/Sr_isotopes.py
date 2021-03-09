#/ Type: DRS
#/ Name: Sr Isotopes Example
#/ Authors: Joe Petrus and Bence Paul
#/ Description: A Sr isotopes example
#/ References: None
#/ Version: 1.0
#/ Contact: support@iolite-software.com

from iolite import QtGui
import numpy as np

def runDRS():

	drs.message("Starting Sr isotopes DRS...")
	drs.progress(0)

	# Get settings
	settings = drs.settings()
	print(settings)

	indexChannel = data.timeSeries(settings["IndexChannel"])
	rmName = settings["ReferenceMaterial"]
	maskOption = settings["Mask"]
	maskChannel = data.timeSeries(settings["MaskChannel"])
	cutoff = settings["MaskCutoff"]
	trim = settings["MaskTrim"]
	# Lambda87 = settings["Lambda87"]
	# Age = settings["Age"]
	RbBias = settings["RbBias"]
	CaArBias = settings["CaArBias"]
	propErrors = settings["PropagateError"]

	# Create debug messages for the settings being used
	IoLog.debug("indexChannelName = %s" % indexChannel.name)
	IoLog.debug("Masking data  = True" if maskOption else "Masking data  = False")
	IoLog.debug("maskChannelName = %s" % maskChannel.name)
	IoLog.debug("maskCutoff = %f" % cutoff)
	IoLog.debug("maskTrim = %f" % trim)
	#IoLog.debug("Lambda87 = %f" % Lambda87)
	#IoLog.debug("Age = %f" % Age)
	IoLog.debug("RbBias = %f" % RbBias)
	IoLog.debug("CaArBias = %f" % CaArBias)
	IoLog.debug("PropagateErrors = True" if propErrors else "PropagateErrors = False")

	# Setup index time
	drs.message("Setting up index time...")
	drs.progress(5)
	drs.setIndexChannel(indexChannel)

	# Setup the mask
	if maskOption:
		drs.message("Making mask...")
		drs.progress(10)
		mask = drs.createMaskFromCutoff(maskChannel, cutoff, trim)
		data.createTimeSeries('mask', data.Intermediate, indexChannel.time(), mask)
	else:
		mask = np.ones_like(indexChannel.data())
		data.createTimeSeries('mask', data.Intermediate, indexChannel.time(), mask)


	# Interp onto index time and baseline subtract
	drs.message("Interpolating onto index time and baseline subtracting...")
	drs.progress(25)

	allInputChannels = data.timeSeriesList(data.Input)
	blGrp = None

	if len(data.selectionGroupList(data.Baseline)) > 1:
		print("There are more than one baseline groups. Sr DRS cannot procede...")
		drs.message("Finished!")
		drs.progress(100)
		drs.finished()
		return
	else:
		blGrp = data.selectionGroupList(data.Baseline)[0]

	for counter, channel in enumerate(allInputChannels):
		drs.message("Baseline subtracting %s" % channel.name)
		drs.progress(25 + 50*counter/len(allInputChannels))

		drs.baselineSubtract(blGrp, [allInputChannels[counter]], mask, 25, 75)

	drs.message("Subtracting interferences...")
	drs.progress(40)

	SrCaAr88 = data.timeSeriesList(data.Intermediate, {'Mass': '88'})[0].data()
	SrCaAr86 = data.timeSeriesList(data.Intermediate, {'Mass': '86'})[0].data()

	SrCaAr86 = data.timeSeriesList(data.Input, {'Mass': '86'})[0].data()

	Sr86channel = data.timeSeriesList(data.Input, {'Mass': '86'})[0]

	SrRb87 = data.timeSeriesList(data.Intermediate, {'Mass': '87'})[0].data()
	Rb85 = data.timeSeriesList(data.Intermediate, {'Mass': '85'})[0].data()
	SrCaAr84 = data.timeSeriesList(data.Intermediate, {'Mass': '84'})[0].data()
	CaAr83 = data.timeSeriesList(data.Intermediate, {'Mass': '83'})[0].data()
	CaAr82 = data.timeSeriesList(data.Intermediate, {'Mass': '82'})[0].data()

	PFract = (np.log(8.37520938 / (SrCaAr88/SrCaAr86))) / (np.log(87.9056/85.9093))*mask
	PSrCaAr86 = SrCaAr86 - (CaAr82 * .004 / .647) / np.power((85.9160721 / 81.9210049), PFract)
	PSrCaAr88 = SrCaAr88 - (CaAr82 * .187 / .647) / np.power((87.9149151 / 81.9210049), PFract)
	Fract = (np.log(8.37520938/(PSrCaAr88/PSrCaAr86))) / (np.log(87.9056/85.9093))*mask
	RbFract = Fract * RbBias
	Rb87 = (Rb85 * 27.8346 / 72.1654) / (np.power((86.90918 / 84.9118), RbFract))
	Sr87 = SrRb87 - Rb87
	CaArFract = Fract * CaArBias
	Ca84 = CaAr82 * 3.22411 / np.power((83.917989 / 81.921122), CaArFract)
	Sr84 = SrCaAr84 - Ca84
	Sr88 = SrCaAr88 - ((CaAr82 * .187 / .647) / np.power((87.9149151 / 81.9210049), Fract))
	Sr86 = SrCaAr86 - ((CaAr82 * .004 / .647) / np.power((85.9160721 / 81.9210049), Fract))

	drs.message("Calculating ratios...")
	drs.progress(60)

	Sr8786_Uncorr = (SrRb87 / SrCaAr86) * np.power((86.9089 / 85.9093), Fract) * mask
	Sr8786_Corr = (Sr87 / Sr86) * np.power((86.9089 / 85.9093), Fract) * mask
	Rb87Sr86ratio = (Rb87 / Sr86) * np.power((86.9089 / 85.9093), Fract) * mask
	Sr8486_Uncorr = (SrCaAr84 / SrCaAr86) * np.power((83.9134 / 85.9093), Fract) * mask
	Sr8486_Corr = (Sr84 / Sr86) * np.power((83.9134 / 85.9093), Fract) * mask
	Rb87asPPM = (Rb87 / SrRb87) * 1000000 * mask
	CaAr84asPPM = (Ca84 / SrCaAr84) * 100000 * mask
	TotalSrBeam = Sr88 + Sr84 + Sr86 + Sr87 * mask
	inRun8283_Ratio = CaAr82 / CaAr83	* mask	#This is a measure of REE vs. CaAr contributions. If pure CaAr = 4.9. If pure REE it's not entirely predictable, based on solution analysis of Indian perovsite = 3.6, Indian WR = 2.4.

	#Gather up intermediate channels and add them as time series:
	int_channel_names = ['PFract', 'PSrCaAr86', 'PSrCaAr88', 'Fract', 'RbFract', 'Rb87', 'Sr87', 'CaArFract', 'Ca84', 'Sr84', 'Sr88', 'Sr86', 'Rb87asPPM', 'CaAr84asPPM', 'inRun8283_Ratio']
	int_channels = [PFract, PSrCaAr86, PSrCaAr88, Fract, RbFract, Rb87, Sr87, CaArFract, Ca84, Sr84, Sr88, Sr86, Rb87asPPM, CaAr84asPPM, inRun8283_Ratio]
	for name, channel in zip(int_channel_names, int_channels):
		data.createTimeSeries(name, data.Intermediate, indexChannel.time(), channel)

	output_channels_names = ['TotalSrBeam', 'Sr8786_Uncorr', 'Sr8786_Corr', 'Rb87Sr86ratio', 'Sr8486_Uncorr', 'Sr8486_Corr']
	output_channels = [TotalSrBeam, Sr8786_Uncorr, Sr8786_Corr, Rb87Sr86ratio, Sr8486_Uncorr, Sr8486_Corr]
	for name, channel in zip(output_channels_names, output_channels):
		data.createTimeSeries(name, data.Output, indexChannel.time(), channel)

	drs.message("Correcting ratios...")
	drs.progress(80)

	StdSpline_Sr87_86 = data.spline(rmName, "Sr8786_Corr").data()
	StdValue_Sr87_86 = data.referenceMaterialData(rmName)["87Sr_86Sr"].value()

	print("StdSpline_Sr87_86 mean = %f"%StdSpline_Sr87_86.mean())
	print("StdValue_Sr87_86 = %f"%StdValue_Sr87_86)

	StdCorr_Sr87_86 = (Sr8786_Corr) * StdValue_Sr87_86 / StdSpline_Sr87_86
	data.createTimeSeries('StdCorr_Sr87_86', data.Output, indexChannel.time(), StdCorr_Sr87_86)


	if propErrors:
		drs.message("Propagating errors...")
		drs.progress(90)

		groups = [s for s in data.selectionGroupList() if s.type != data.Baseline]
		data.propagateErrors(groups, [data.timeSeries("StdCorr_Sr87_86")], data.timeSeries("Sr8786_Corr"), rmName)


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

	timeSeriesNames = data.timeSeriesNames(data.Input)
	defaultChannelName = ""
	if timeSeriesNames:
		defaultChannelName = timeSeriesNames[0]

	rmNames = data.selectionGroupNames(data.ReferenceMaterial)

	drs.setSetting("IndexChannel", defaultChannelName)
	drs.setSetting("ReferenceMaterial", "C_MMC")
	drs.setSetting("Mask", False)
	drs.setSetting("MaskChannel", defaultChannelName)
	drs.setSetting("MaskCutoff", 0.1)
	drs.setSetting("MaskTrim", 0.0)
	#drs.setSetting("Lambda87", 1.42e-11)
	#drs.setSetting("Age", 0.)
	drs.setSetting("RbBias", 1.)
	drs.setSetting("CaArBias", 1.)
	drs.setSetting("PropagateError", False)

	settings = drs.settings()

	indexComboBox = QtGui.QComboBox(widget)
	indexComboBox.addItems(timeSeriesNames)
	indexComboBox.setCurrentText(settings["IndexChannel"])
	indexComboBox.currentTextChanged.connect(lambda t: drs.setSetting("IndexChannel", t))
	formLayout.addRow("Index channel", indexComboBox)

	rmComboBox = QtGui.QComboBox(widget)
	rmComboBox.addItems(rmNames)
	if settings["ReferenceMaterial"] in rmNames:
		rmComboBox.setCurrentText(settings["ReferenceMaterial"])
	else:
		rmComboBox.setCurrentText(rmNames[0])
		drs.setSetting("ReferenceMaterial", rmNames[0])
	rmComboBox.currentTextChanged.connect(lambda t: drs.setSetting("ReferenceMaterial", t))
	formLayout.addRow("Reference material", rmComboBox)

	verticalSpacer = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
	formLayout.addItem(verticalSpacer)

	maskCheckBox = QtGui.QCheckBox(widget)
	maskCheckBox.setChecked(settings["Mask"])
	maskCheckBox.toggled.connect(lambda t: drs.setSetting("Mask", bool(t)))
	formLayout.addRow("Mask", maskCheckBox)

	maskComboBox = QtGui.QComboBox(widget)
	maskComboBox.addItems(data.timeSeriesNames(data.Input))
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

	verticalSpacer2 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
	formLayout.addItem(verticalSpacer2)

	# lambda87LineEdit = QtGui.QLineEdit(widget)
	# lambda87LineEdit.setText(settings["Lambda87"])
	# lambda87LineEdit.textChanged.connect(lambda t: drs.setSetting("Lambda87", float(t)))
	# formLayout.addRow("Lambda87", lambda87LineEdit)

	# ageLineEdit = QtGui.QLineEdit(widget)
	# ageLineEdit.setText(settings["Age"])
	# ageLineEdit.textChanged.connect(lambda t: drs.setSetting("Age", float(t)))
	# formLayout.addRow("Age", ageLineEdit)

	rbBiasLineEdit = QtGui.QLineEdit(widget)
	rbBiasLineEdit.setText(settings["RbBias"])
	rbBiasLineEdit.textChanged.connect(lambda t: drs.setSetting("RbBias", float(t)))
	formLayout.addRow("Rb Bias", rbBiasLineEdit)

	caArBiasLineEdit = QtGui.QLineEdit(widget)
	caArBiasLineEdit.setText(settings["CaArBias"])
	caArBiasLineEdit.textChanged.connect(lambda t: drs.setSetting("CaArBias", float(t)))
	formLayout.addRow("CaAr Bias", caArBiasLineEdit)

	propCheckBox = QtGui.QCheckBox(widget)
	propCheckBox.setChecked(settings["PropagateError"])
	propCheckBox.toggled.connect(lambda t: drs.setSetting("PropagateError", bool(t)))
	formLayout.addRow("Propagate Errors?", propCheckBox)

	drs.setSettingsWidget(widget)
