#/ Type: DRS
#/ Name: Sm-Nd Isotopes with Downhole fractionation correction 
#/ Authors: Vitor Barrote, Chris Fisher and Joe Petrus
#/ Description: Sm-Nd I4 version of the I3 DRS created by Chris Fisher and Chad Paton
#/ DOI: 10.5281/zenodo.5512126
#/ References: Fisher et al. 2020
#/ Version: 1.0
#/ Contact: vitorbarrote@hotmail.com

from iolite import QtGui
from iolite.Qt import Qt, QColor
from iolite.ui import IolitePlotPyInterface as Plot
import numpy as np
from scipy.optimize import curve_fit, leastsq

def downholeFunc(t, a, b, c, d):
	return a + b*t + c*np.exp(-d * t)

def runDRS():

	drs.message("Starting Sm-Nd Downhole fract DRS...")
	drs.progress(0)

# Get settings
	settings = drs.settings()
	print(settings)   

	indexChannel = data.timeSeries(settings["IndexChannel"])
	maskChannel = data.timeSeries(settings["MaskChannel"])
	rmName = settings["ReferenceMaterial"]
	cutoff = settings["MaskCutoff"]
	trim = settings["MaskTrim"]
	NdTrue = settings["NdTrue"]
	Sm147_149 = settings["Sm147_149"]
	Sm144_149 = settings["Sm144_149"]
	Sm148_149 = settings["Sm148_149"]
	Sm150_149 = settings["Sm150_149"]
	#age = settings["Age"]
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

	#error prop?
	commonProps = {'DRS': drs.name()}

	# Setup the mask
	drs.message("Making mask...")
	drs.progress(10)
	mask = drs.createMaskFromCutoff(maskChannel, cutoff, trim)

	
# Interp onto index time and baseline subtract
	drs.message("Interpolating onto index time and baseline subtracting...")
	drs.progress(25)

	allInputChannels = data.timeSeriesList(data.Input)
	
	for counter, channel in enumerate(allInputChannels):
		drs.message("Baseline subtracting %s" % channel.name)
		drs.progress(25 + 50*counter/len(allInputChannels))

		drs.baselineSubtract(data.selectionGroup("Baseline"), [allInputChannels[counter]], mask, 25, 75)



#Change the name of the channel depending of what the output from the instrument is. I turned off NdCe142 because we can't correct for Ce.

	#NdCe142 = 	data.timeSeries("142_CPS").data()
	Nd143 = data.timeSeries("Nd143_CPS").data()
	NdSm144 = data.timeSeries("Nd144_CPS").data()
	Nd145 = data.timeSeries("Nd145_CPS").data()
	Nd146 = data.timeSeries("Nd146_CPS").data()
	Sm147 = data.timeSeries("Sm147_CPS").data()
	NdSm148 = data.timeSeries("Nd148_CPS").data()
	Sm149 = data.timeSeries("Sm149_CPS").data()
	NdSm150 = data.timeSeries("Nd150_CPS").data()
	
	SmFract = (np.log(Sm147_149/ (Sm147 / Sm149))) / (np.log(146.9149 / 148.9172 ))	*mask
	Sm144 = Sm149 * Sm144_149 *(( 148.9172 / 143.912) ** SmFract) 	*mask
	Sm148 = Sm149 * Sm148_149 *(( 148.9172 / 147.9148 ) ** SmFract) *mask
	Sm150 = Sm149 * Sm150_149 *(( 148.9172 / 149.91727 ) ** SmFract)*mask
	Nd144c = (NdSm144 - Sm144)*mask
	Nd148c = (NdSm148 - Sm148)*mask
	Nd150c = (NdSm150 - Sm150)
	NdFract = (np.log(NdTrue / (Nd146 / Nd144c))) / (np.log(145.9131 / 143.9098))*mask
	NdFractJNdi= (np.log(NdTrue / (Nd146 / NdSm144))) / (np.log(145.9131 / 143.9098))*mask	
	Nd143_144_Raw = (Nd143 / Nd144c) *mask
	Nd143_144_Corr = (Nd143 / Nd144c) * (142.9098 / 143.9098) ** NdFract *mask
	#Nd143_144_AgeCorr=Nd143_144_Corr-Sm147_Nd144_Raw*(exp(1000000*lambdaSm *age)-1)
	Nd143_144_Corr_JNdi = (Nd143 / NdSm144) * (142.9098 / 143.9098) ** NdFractJNdi *mask
	Nd145_144_Raw = (Nd145 / Nd144c) *mask
	Nd145_144_Corr = (Nd145 / Nd144c) * (144.9126 / 143.9098) ** NdFract *mask
	Nd145_144_Corr_JNdi = (Nd145 / NdSm144) * (144.9126 / 143.9098) ** NdFractJNdi *mask
	Nd148_144_Raw = (Nd148c / Nd144c) *mask
	Nd148_144_Corr = (Nd148c / Nd144c) * ( 147.916889 / 143.9098) ** NdFract *mask
	Nd148_144_Corr_JNdi = (NdSm148 / NdSm144) * (147.916889 / 143.9098) ** NdFractJNdi *mask
	Nd150_144_Raw = (Nd150c / Nd144c) *mask
	Nd150_144_Corr = (Nd150c / Nd144c) * ( 149.920887 / 143.9098) ** NdFract *mask
	Nd150_144_Corr_JNdi = (NdSm150 / NdSm144) * (149.920887 / 143.9098) ** NdFractJNdi *mask
	Sm147_Nd144_Corr = (Sm147 / Nd144c) *mask
	Nd146_144_Raw = (Nd146 / Nd144c)*mask
	Nd146_144_Raw_JNdi = (Nd146 / NdSm144) *mask
	#Nd142 = Nd146Beam * (0.272 / 0.172)	*  (145.9131 / 141.907719) ** NdFract *mask
	#Ce142 = (NdCe142Beam - Nd142) *mask
	#Ce_AbN = Ce142 * (1 / .1108) *mask
	NdFract_SmFract =  NdFract/ SmFract *mask
	TotalNd = Nd146/0.172 *mask

#DF correction

	method = drs.setting('BeamSecondsMethod')
	print(f'Using beam seconds method {method}')

	try:
		channel = data.timeSeries(drs.setting('BeamSecondsChannel'))
		value = float(drs.setting('BeamSecondsValue'))
	except Exception as e:
		if 'threshold' in method or 'change' in method:
			raise Exception(e)
	else:
			pass

	if 'log' in method:
		drs.createBeamSecondsFromLaserLog()
	elif 'samples' in method:
		drs.createBeamSecondsFromSamples()
	elif 'threshold' in method:
		drs.createBeamSecondsFromCutoff(channel, value)
	elif 'change' in method:
		drs.createBeamSecondsFromJump(channel, value)

	beamSeconds = data.timeSeries('BeamSeconds').data()

	timeStep = indexChannel.time()[1] - indexChannel.time()[0]
	startTrimSec = settings["StartTrim"]
	startTrimIndex = int(round(startTrimSec/timeStep))
	endTrimSec = settings["EndTrim"]
	endTrimIndex = int(round(endTrimSec/timeStep))

# Clear previous plots
	settings['FitWidget'].clearGraphs()

	drs.message('Working on ' + "Sm147/Nd144")
	drs.progress(50)

	data.createTimeSeries("Sm147_Nd144_Corr", data.Output, indexChannel.time(), Sm147_Nd144_Corr)
	rawRatio = data.timeSeries('Sm147_Nd144_Corr').data()
	rawRatio[np.isinf(rawRatio)] = np.nan
	ts = data.createTimeSeries("Sm147/Nd144", data.Intermediate, indexChannel.time(), rawRatio, commonProps)

	DHF = data.compileDownhole(data.selectionGroup(settings['ReferenceMaterial']), ts)
	DHFnans = np.isnan(DHF[1])
	DHFt = DHF[0][~DHFnans]
	DHFr = DHF[1][~DHFnans]

	if startTrimIndex != 0:
			DHFt = DHFt[startTrimIndex:]
			DHFr = DHFr[startTrimIndex:]

	if endTrimIndex != 0:
			DHFt = DHFt[:-endTrimIndex]
			DHFr = DHFr[:-endTrimIndex]
        
	params, cov = curve_fit(downholeFunc, DHFt, DHFr, ftol=1e-5)
	dc = rawRatio/(1 + (params[1]/params[0])*beamSeconds + (params[2]/params[0])*np.exp(-params[3]*beamSeconds))
	data.createTimeSeries('DC '+"Sm147/Nd144", data.Intermediate, indexChannel.time(), dc, commonProps)

	plot = settings['FitWidget']
	g = plot.addGraph()
	g.setData(DHFt, DHFr)
	g2 = plot.addGraph()
	g2.setColor(QColor(255, 0, 0))
	g2.setData(DHFt, downholeFunc(DHFt, params[0], params[1], params[2], params[3]))
	plot.left().label = "Sm147/Nd144"
	plot.bottom().label = 'Time (s)'
	plot.setToolsVisible(False)
	plot.rescaleAxes()
	plot.replot()

#Create output channels and correct to the RM ******************************************************
	rm = data.referenceMaterialData(settings["ReferenceMaterial"])
	rmValue = rm["147Sm/144Nd"].value()
	rmSpline = data.spline(settings['ReferenceMaterial'], 'Sm147_Nd144_Corr').data()
	finalRatio = (rmValue/rmSpline)*dc
	data.createTimeSeries('Final ' + "Sm147/Nd144", data.Output, indexChannel.time(), finalRatio, commonProps)

	data.createTimeSeries("Nd143_144_Corr", data.Output, indexChannel.time(), Nd143_144_Corr)
	
	StdSpline_Nd143_144 = data.spline(rmName, "Nd143_144_Corr").data()
	StdValue_Nd143_144 = data.referenceMaterialData(rmName)["143Nd/144Nd"].value()
	StdCorr_Nd143_144 = (Nd143_144_Corr)* StdValue_Nd143_144 / StdSpline_Nd143_144
	data.createTimeSeries("StdCorr_Nd143_144", data.Output, indexChannel.time(), StdCorr_Nd143_144)

	data.createTimeSeries("Nd143_144_Corr_JNdi", data.Output, indexChannel.time(), Nd143_144_Corr_JNdi)

	data.createTimeSeries("Sm147_Nd144_Corr", data.Output, indexChannel.time(), Sm147_Nd144_Corr)
	StdSpline_Sm147_Nd144 = data.spline(rmName, "Sm147_Nd144_Corr").data()
	StdValue_Sm147_Nd144 = data.referenceMaterialData(rmName)["147Sm/144Nd"].value()

	data.createTimeSeries("Nd145_144_Corr", data.Output, indexChannel.time(), Nd145_144_Corr)
	StdSpline_Nd145_144 = data.spline(rmName, "Nd145_144_Corr").data()
	StdValue_Nd145_144 = data.referenceMaterialData(rmName)["145Nd/144Nd"].value()

	print("StdSpline_Nd143_144 mean = %f"%StdSpline_Nd143_144.mean())
	print("StdValue_Nd143_144 = %f"%StdValue_Nd143_144)


#error propagation from the RM correction
	if propErrors:
		groups = [s for s in data.selectionGroupList() if s.type != data.Baseline]
		data.propagateErrors(groups, [data.timeSeries("StdCorr_Nd143_144")], data.timeSeries("Nd143_144_Corr"), rmName)
		data.propagateErrors(groups, [data.timeSeries("Final Sm147/Nd144")], data.timeSeries("DC Sm147/Nd144"), rmName)

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
	formLayout.setFieldGrowthPolicy(QtGui.QFormLayout.AllNonFixedFieldsGrow )
	formLayout.setFormAlignment(Qt.AlignHCenter | Qt.AlignTop)
	widget.setLayout(formLayout)


	timeSeriesNames = data.timeSeriesNames(data.Input)
	defaultChannelName = ""
	if timeSeriesNames:
		defaultChannelName = timeSeriesNames[0]

	rmNames = data.selectionGroupNames(data.ReferenceMaterial)
	
	drs.setSetting("IndexChannel", "Nd146")
	drs.setSetting("ReferenceMaterial", "M_STK")
	drs.setSetting("MaskChannel", "Nd146")
	drs.setSetting("MaskCutoff", 0.1)
	drs.setSetting("MaskTrim", 0.0)
	drs.setSetting("StartTrim", 0.1)
	drs.setSetting("EndTrim", 0.1)
	drs.setSetting("NdTrue", 0.7219)
	drs.setSetting("Sm147_149", 1.08680)
	drs.setSetting("Sm144_149", 0.22332)
	drs.setSetting("Sm148_149", 0.81419)
	drs.setSetting("Sm150_149", 0.53366)
	#drs.setSetting("Age", 0)
	drs.setSetting("BeamSecondsMethod", "Laser log")
	drs.setSetting("BeamSecondsChannel", "Nd146")
	drs.setSetting("BeamSecondsValue", 50)
	drs.setSetting("PropagateError", True)

	settings = drs.settings()
	
	indexComboBox = QtGui.QComboBox(widget)
	indexComboBox.setFixedWidth(150)
	indexComboBox.addItems(timeSeriesNames)
	indexComboBox.setCurrentText(settings["IndexChannel"])
	indexComboBox.currentTextChanged.connect(lambda t: drs.setSetting("IndexChannel", t))
	formLayout.addRow("Index channel", indexComboBox)

	rmComboBox = QtGui.QComboBox(widget)
	rmComboBox.setFixedWidth(150)
	rmComboBox.addItems(data.referenceMaterialNames())
	rmComboBox.currentTextChanged.connect(lambda s: drs.setSetting("ReferenceMaterial", str(s)))
	formLayout.addRow("Reference material", rmComboBox)

	maskComboBox = QtGui.QComboBox(widget)
	maskComboBox.setFixedWidth(150)
	maskComboBox.addItems(data.timeSeriesNames(data.Input))
	maskComboBox.setCurrentText(settings["MaskChannel"])
	maskComboBox.currentTextChanged.connect(lambda t: drs.setSetting("MaskChannel", t))
	formLayout.addRow("Mask channel", maskComboBox)

	maskLineEdit = QtGui.QLineEdit(widget)
	maskLineEdit.setFixedWidth(150)
	maskLineEdit.setText(settings["MaskCutoff"])
	maskLineEdit.textChanged.connect(lambda t: drs.setSetting("MaskCutoff", float(t)))
	formLayout.addRow("Mask cutoff", maskLineEdit)

	#maskTrimLineEdit = QtGui.QLineEdit(widget)
	#maskTrimLineEdit.setText(settings["MaskTrim"])
	#maskTrimLineEdit.textChanged.connect(lambda t: drs.setSetting("MaskTrim", float(t)))
	#formLayout.addRow("Mask trim", maskTrimLineEdit)
 
	startLineEdit = QtGui.QLineEdit(widget)
	startLineEdit.setFixedWidth(150)
	formLayout.addRow("Start trim (s)", startLineEdit)
	startLineEdit.textEdited.connect(lambda s: drs.setSetting("StartTrim", float(s)))

	endLineEdit = QtGui.QLineEdit(widget)
	endLineEdit.setFixedWidth(150)
	formLayout.addRow("End trim (s)", endLineEdit)
	endLineEdit.textEdited.connect(lambda s: drs.setSetting("EndTrim", float(s)))

	ndTrueLineEdit = QtGui.QLineEdit(widget)
	ndTrueLineEdit.setFixedWidth(150)
	ndTrueLineEdit.setText(settings["NdTrue"])
	ndTrueLineEdit.textChanged.connect(lambda t: drs.setSetting("NdTrue", float(t)))
	formLayout.addRow("NdTrue", ndTrueLineEdit)

	Sm147_149LineEdit = QtGui.QLineEdit(widget)
	Sm147_149LineEdit.setFixedWidth(150)
	Sm147_149LineEdit.setText(settings["Sm147_149"])
	Sm147_149LineEdit.textChanged.connect(lambda t: drs.setSetting("Sm147_149", float(t)))
	formLayout.addRow("Sm147_149", Sm147_149LineEdit)

	Sm144_149LineEdit = QtGui.QLineEdit(widget)
	Sm144_149LineEdit.setFixedWidth(150)
	Sm144_149LineEdit.setText(settings["Sm144_149"])
	Sm144_149LineEdit.textChanged.connect(lambda t: drs.setSetting("Sm144_149", float(t)))
	formLayout.addRow("Sm144_149", Sm144_149LineEdit)
	
	Sm148_149LineEdit = QtGui.QLineEdit(widget)
	Sm148_149LineEdit.setFixedWidth(150)
	Sm148_149LineEdit.setText(settings["Sm148_149"])
	Sm148_149LineEdit.textChanged.connect(lambda t: drs.setSetting("Sm148_149", float(t)))
	formLayout.addRow("Sm148_149", Sm148_149LineEdit)
	
	Sm150_149LineEdit = QtGui.QLineEdit(widget)
	Sm150_149LineEdit.setFixedWidth(150)
	Sm150_149LineEdit.setText(settings["Sm150_149"])
	Sm150_149LineEdit.textChanged.connect(lambda t: drs.setSetting("Sm150_149", float(t)))
	formLayout.addRow("Sm150_149", Sm150_149LineEdit)

	#ageLineEdit = QtGui.QLineEdit(widget)
	#ageLineEdit.setText(settings["Age"])
	#ageLineEdit.textChanged.connect(lambda t: drs.setSetting("Age", float(t)))
	#formLayout.addRow("Age", ageLineEdit)

	bsMethodComboBox = QtGui.QComboBox(widget)
	bsMethodComboBox.setFixedWidth(150)
	bsMethodComboBox.addItems(['Laser log', 'Gaps between samples', 'Cutoff threshold', 'Rate of change'])
	bsMethodComboBox.activated.connect(lambda v: drs.setSetting('BeamSecondsMethod', bsMethodComboBox.currentText))
	bsMethodComboBox.setCurrentText(drs.setting('BeamSecondsMethod'))
	formLayout.addRow('Beam seconds method', bsMethodComboBox)

	bsChannelComboBox = QtGui.QComboBox(widget)
	bsChannelComboBox.setFixedWidth(150)
	bsChannelComboBox.addItems(data.timeSeriesNames(data.Input))
	bsChannelComboBox.activated.connect(lambda v: drs.setSetting('BeamSecondsChannel', bsChannelComboBox.currentText))
	bsChannelComboBox.setCurrentText(drs.setting('BeamSecondsChannel'))
	formLayout.addRow('Beam seconds channel', bsChannelComboBox)

	bsValueLineEdit = QtGui.QLineEdit(widget)
	bsValueLineEdit.setFixedWidth(150)
	bsValueLineEdit.textEdited.connect(lambda v: drs.setSetting('BeamSecondsValue', float(v)))
	bsValueLineEdit.setText(drs.setting('BeamSecondsValue'))
	formLayout.addRow('Beam seconds value', bsValueLineEdit)

	plotWidget = Plot(widget)
	plotWidget.setFixedSize(500,400)
	formLayout.addRow('Down-hole fit', plotWidget)
	drs.setSetting('FitWidget', plotWidget)

	propCheckBox = QtGui.QCheckBox(widget)
	propCheckBox.setChecked(settings["PropagateError"])
	propCheckBox.toggled.connect(lambda t: drs.setSetting("PropagateError", bool(t)))
	formLayout.addRow("PropagateError", propCheckBox)

    # Restore settings
	try:
		settings = drs.settings()
		print('Restoring settings...')
		print(settings)
		rmComboBox.setCurrentText(settings["ReferenceMaterial"])
		startLineEdit.setText(str(settings['StartTrim']))
		endLineEdit.setText(str(settings['EndTrim']))
		bsMethodComboBox.setCurrentText(settings['BeamSecondsMethod'])
		bsChannelComboBox.setCurrentText(settings['BeamSecondsChannel'])
		bsValueLineEdit.setText(str(settings['BeamSecondsValue']))
	except KeyError:
		pass

	drs.setSettingsWidget(widget)
