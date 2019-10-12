#/ Type: DRS
#/ Name: Trace Elements Sum-Norm Multi-RM
#/ Authors: Joe Petrus and Bence Paul
#/ Description: Sum-normalizing trace elements data reduction scheme with multiple external RMs
#/ References: None
#/ Version: 1.0
#/ Contact: support@iolite-software.com

from iolite import QtGui
from iolite.Qt import Qt
import numpy as np
import re
#from sklearn.linear_model import HuberRegressor

oxide_factors = {
	"Ag2O": 1.0741,
	"Al2O3": 1.8895,
	"As2O3": 1.3203,
	"As2O5": 1.5339,
	"Au2O": 1.0406,
	"B2O3": 3.2202,
	"BaO": 1.1165,
	"BeO": 2.7758,
	"Bi2O5": 1.1914,
	"CO2": 3.6644,
	"CaO": 1.3992,
	"CdO": 1.1423,
	"Ce2O3": 1.1713,
	"CeO2": 1.2284,
	"CoO": 1.2715,
	"Cr2O3": 1.4615,
	"Cs2O": 1.0602,
	"CuO": 1.2518,
	"Dy2O3": 1.1477,
	"Er2O3": 1.1435,
	"Eu2O3": 1.1579,
	"FeO": 1.2865,
	"Fe2O3": 1.4297,
	"Ga2O3": 1.3442,
	"Gd2O3": 1.1526,
	"GeO2": 1.4408,
	"HfO2": 1.1793,
	"HgO": 1.0798,
	"Ho2O3": 1.1455,
	"In2O3": 1.2091,
	"IrO": 1.0832,
	"K2O": 1.2046,
	"La2O3": 1.1728,
	"Li2O": 2.1527,
	"Lu2O3": 1.1371,
	"MgO": 1.6582,
	"MnO": 1.2912,
	"MnO2": 1.5825,
	"MoO3": 1.5003,
	"N2O5": 3.8551,
	"Na2O": 1.348,
	"Nb2O5": 1.4305,
	"Nd2O3": 1.1664,
	"NiO": 1.2725,
	"OsO": 1.0841,
	"P2O5": 2.2916,
	"PbO": 1.0772,
	"PbO2": 1.1544,
	"PdO": 1.1504,
	"Pr2O3": 1.1703,
	"Pr6O11": 1.2082,
	"PtO": 1.082,
	"Rb2O": 1.0936,
	"ReO": 1.0859,
	"RhO": 1.5555,
	"RuO": 1.1583,
	"SO3": 2.4972,
	"Sb2O5": 1.3284,
	"Sc2O3": 1.5338,
	"SeO3": 1.6079,
	"SiO2": 2.1392,
	"Sm2O3": 1.1596,
	"SnO2": 1.2696,
	"SrO": 1.1826,
	"Ta2O5": 1.2211,
	"Tb2O3": 1.151,
	"Tb4O7": 1.1762,
	"TeO3": 1.3762,
	"ThO2": 1.1379,
	"TiO2": 1.6681,
	"Tl2O3": 1.1174,
	"Tm2O3": 1.1421,
	"UO2": 1.1344,
	"UO3": 1.2017,
	"U3O8": 1.1792,
	"V2O5": 1.7852,
	"WO3": 1.261,
	"Y2O3": 1.2699,
	"Yb2O3": 1.1387,
	"ZnO": 1.2448,
	"ZrO2":1.3508
}

def get_oxide_factor(element):
	prog = re.compile('^([A-Z][a-z]?)')
	matching_keys = [k for k in oxide_factors.keys() if prog.match(k).group(0) == element]
	if matching_keys:
		return oxide_factors[matching_keys[0]]
	else:
		print('Could not get oxide factor for element %s'%(element))
		return 0

def runDRS():
	drs.message("Starting DRS...")
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
	drs.baselineSubtract(data.selectionGroup("Baseline"), allInputChannels, None, 25, 50)
	
	commonProps = {'DRS': drs.name()}

	# Work out relative ablation yields (not yet time-resolved)
	drs.message('Working on relative ablation yields...')
	drs.progress(51)
	
	groupNames = settings["ExtraExternals"]
	masterGroupName = settings["MasterExternal"]
	masterGroup = data.selectionGroup(masterGroupName)
	
	cpsChannelNames = data.timeSeriesNames(data.Intermediate, {'DRSType': 'BaselineSubtracted'})
	cpsChannels = [data.timeSeries(c) for c in cpsChannelNames]
	
	M = np.empty( (len(cpsChannels), len(groupNames)) )
	M.fill(np.nan)
	
	for col, groupName in enumerate(groupNames):
		group = data.selectionGroup(groupName)
		for row, channel in enumerate(cpsChannels):
			channelElement = channel.property('Element')
			try:
				groupResult = data.groupResult(group, channel).value()
				groupRMValue = data.referenceMaterialData(groupName)[channelElement].valueInPPM()
				masterGroupResult = data.groupResult(masterGroup, channel).value()
				masterGroupRMValue = data.referenceMaterialData(masterGroupName)[channelElement].valueInPPM()
				M[row][col] = (groupResult/groupRMValue) / (masterGroupResult/masterGroupRMValue)
			except:
				pass
			
	ablationFactors = dict(zip(groupNames, np.nanmedian(M, axis=0)))
	print(ablationFactors)
	
	# Work out SQ concentrations using 
	drs.message('Calculating semi-quant concentrations...')
	drs.progress(60)
	
	for cpsChannel in cpsChannels:
		channelElement = cpsChannel.property('Element')
		groupsForChannel = [g for g in groupNames if channelElement in data.referenceMaterialData(g)]
		N = len(groupsForChannel)
		if N < 1:
			IoLog.warning('None of the reference materials have data for element %s'%(channelElement))
			continue
		
		ppmValues = [data.referenceMaterialData(g)[channelElement].valueInPPM() for g in groupsForChannel]
		cpsValues = np.empty( (len(cpsChannel.data()), len(groupsForChannel)) )
		
		for gi, g in enumerate(groupsForChannel):
			gSpline = data.spline(g, cpsChannel.name).data() / ablationFactors[g]
			cpsValues[:,gi] = gSpline
		
		def doFit(x):
			return x.dot(ppmValues)/x.dot(x)
		
		slopes = np.apply_along_axis(doFit, axis=1, arr=cpsValues)
		
		data.createTimeSeries('%s_slope'%(cpsChannel.name), data.Intermediate, None, slopes, commonProps)
		
		sq = cpsChannel.data() * slopes
		
		data.createTimeSeries('%s%i_ppm'%(cpsChannel.property('Element'), int(cpsChannel.property('Mass'))), data.Output, None, sq, commonProps)
		
	# Work out normalizing factor
	channels_for_norm = [c + '_ppm' for c in settings['Elements']]
	channels_sum = np.zeros(len(indexChannel.data()))
	use_oxides = settings['Oxides']
	for channel in channels_for_norm:
		# Need to add in oxide ability!!
		element = data.timeSeries(channel).property('Element')
		type_factor = get_oxide_factor(element) if use_oxides else 1
		#print('%s: %s: %f'%(channel, element, type_factor))
		channels_sum += data.timeSeries(channel).data()*type_factor
	
	factor = (1e6*settings['Value']/100)/channels_sum #(channels_sum/1e4)/settings['Value']
	data.createTimeSeries("NormalizationFactor", data.Intermediate, None, factor, commonProps)
	
	maskChannelData = data.timeSeries(settings['MaskChannel']).data()
	maskValue = float(settings['MaskValue'])
	mask = np.ones(len(maskChannelData))
	mask[maskChannelData < maskValue] = np.nan
	
	channels_to_adjust = [c for c in data.timeSeriesList(data.Output) if 'ppm' in c.name]
	print(channels_to_adjust)
	for channel in channels_to_adjust:
		print('Adjusting %s'%(channel.name))
		d = mask*channel.data()*factor
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
	rmComboBox.currentTextChanged.connect(lambda s: drs.setSetting("MasterExternal", str(s)))
	formLayout.addRow("Master reference material", rmComboBox)
	
	extraRMsList = QtGui.QListWidget(widget)
	extraRMsList.addItems(data.selectionGroupNames(data.ReferenceMaterial))
	extraRMsList.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
	def updateExtraGroups():
		groups = [wi.text() for wi in extraRMsList.selectedItems()]
		drs.setSetting('ExtraExternals', groups)
	extraRMsList.itemSelectionChanged.connect(updateExtraGroups)
	formLayout.addRow('Extra reference materials', extraRMsList)

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
	
	maskComboBox = QtGui.QComboBox(widget)
	maskComboBox.addItems(data.timeSeriesNames())
	maskComboBox.currentTextChanged.connect(lambda s: drs.setSetting("MaskChannel", str(s)))
	formLayout.addRow("Mask channel", maskComboBox)
	
	maskLineEdit = QtGui.QLineEdit(widget)
	maskLineEdit.textChanged.connect(lambda v: drs.setSetting("MaskValue", float(v)))
	formLayout.addRow("Mask value", maskLineEdit)
	
	
	# Restore settings
	try:
		settings = drs.settings()
		rmComboBox.setCurrentText(settings["MasterExternal"])
		extraRMsList.selectionModel().clear()
		for rm in settings["ExtraExternals"]:
			matches = extraRMsList.findItems(rm, Qt.MatchFixedString)
			if matches:
				matches[0].setSelected(True)
		elementsList.selectionModel().clear()
		for el in settings["Elements"]:
			matches = elementsList.findItems(el, Qt.MatchFixedString)
			if matches:
				matches[0].setSelected(True)
		oxideCheckBox.setChecked(settings["Oxides"])
		valueLineEdit.setText(str(float(settings["Value"])))
		maskComboBox.setCurrentText(settings["MaskChannel"])
		maskLineEdit.setText(str(float(settings["MaskValue"])))
	except KeyError:
		pass

	drs.setSettingsWidget(widget)
