#/ Type: DRS
#/ Name: LuHf_Geochron
#/ Authors: Kai Rankenburg & Bruno Ribeiro
#/ Description: Lu-Hf using NH3
#/ References: Ribeiro et al. (2023), https://doi.org/10.1016/j.epsl.2022.117969
#/ Version: 1.0
#/ Contact: support@iolite-software.com

from iolite import QtGui
import numpy as np
from iolite.types import Result


# Correlation coeffiient Associated Results
def calcRho1(s):
    Lu176_Hf177s = data.timeSeries("Lu176_Hf177").dataForSelection(s)
    Hf176_177s = data.timeSeries("Hf176_177").dataForSelection(s)
    rho1 = np.corrcoef(Hf176_177s, Lu176_Hf177s, rowvar=False)
    res1 = Result(rho1[0,1], 0.001, 'ratio', 'rho Hf176_177 v LuHf')
    return res1

def calcRho2(s):
    Lu176_Hf176s = data.timeSeries("Lu176_Hf176").dataForSelection(s)
    Hf177_176s = data.timeSeries("Hf177_176").dataForSelection(s)
    rho2 = np.corrcoef(Hf177_176s, Lu176_Hf176s, rowvar=False)
    res2 = Result(rho2[0, 1], 0.001, 'ratio', 'rho Hf177_176 v LuHf')
    return res2


def runDRS():

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
    Lu176_175ref = settings["Lu176_175ref"]
    Yb176_172ref = settings["Yb176_172ref"]
    Hf177_178ref = settings["Hf177_178ref"]
    Hf176_178ref = settings["Hf176_178ref"]
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

    allInputChannels = data.timeSeriesList(data.Input)
    blGrp = None

    if len(data.selectionGroupList(data.Baseline)) > 1:
        IoLog.error("There are multiple baseline groups. Lu-Hf DRS cannot proceed...")
        drs.message("Error. See Messages")
        drs.progress(100)
        drs.finished()
        return
    elif len(data.selectionGroupList(data.Baseline)) < 1:
        IoLog.error("No baselines. Please select some baselines. Lu-Hf DRS cannot proceed...")
        drs.message("Error. See Messages")
        drs.progress(100)
        drs.finished()
        return
    else:
        blGrp = data.selectionGroupList(data.Baseline)[0]

    if len(blGrp.selections()) < 1:
        IoLog.error("No baseline selections. Please select some baselines. Lu-Hf DRS cannot proceed...")
        drs.message("Error. See Messages")
        drs.progress(100)
        drs.finished()
        return

    drs.baselineSubtract(blGrp, allInputChannels, mask, 25, 75)

    Yb254 = data.timeSeriesList(data.Intermediate, {'Post-shift mass': '254'})[0].data()
    Lu175 = data.timeSeriesList(data.Intermediate, {'Post-shift mass': '175'})[0].data()
    #Lu191 = data.timeSeries("Lu191_CPS").data()
    Lu257 = data.timeSeriesList(data.Intermediate, {'Post-shift mass': '257'})[0].data()
    Hf258 = data.timeSeriesList(data.Intermediate, {'Post-shift mass': '258'})[0].data()
    Hf260 = data.timeSeriesList(data.Intermediate, {'Post-shift mass': '260'})[0].data()

    Lu176_Hf177 = (Lu175 * Lu176_175ref) / (Hf260 * Hf177_178ref)
    Lu176_Hf176 = (Lu175 * Lu176_175ref) / (Hf258)

    Hf176_177 = (Hf258-Lu257 * Lu176_175ref - Yb254 * Yb176_172ref) / Hf260 * Hf177_178ref
    Hf177_176 = (Hf260 * Hf177_178ref) / Hf258

    data.createTimeSeries("Lu176_Hf177", data.Output, indexChannel.time(), Lu176_Hf177)
    data.createTimeSeries("Hf176_177", data.Output, indexChannel.time(), Hf176_177)
    data.createTimeSeries("Lu176_Hf176", data.Output, indexChannel.time(), Lu176_Hf176)
    data.createTimeSeries("Hf177_176", data.Output, indexChannel.time(), Hf177_176)

    # Check that ref mat has 176Lu/177Hf value
    rm = data.referenceMaterialData(rmName)
    if not "176Lu/177Hf" in rm.keys():
        IoLog.error(f"There was no 176Lu/177Hf in the RM file \'{rmName}\'")
        drs.message("Error. See Messages")
        drs.progress(100)
        drs.finished()
        return

    # Assumes other values are present

    StdValue_Lu176_Hf177 = data.referenceMaterialData(rmName)["176Lu/177Hf"].value()
    StdValue_Hf176_177 = data.referenceMaterialData(rmName)["176Hf/177Hf"].value()
    StdValue_Lu176_Hf176 = StdValue_Lu176_Hf177 / StdValue_Hf176_177
    StdValue_Hf177_176 = 1 / StdValue_Hf176_177

    StdSpline_Lu176_Hf177 = data.spline(rmName, "Lu176_Hf177").data()
    StdCorr_Lu176_Hf177 = Lu176_Hf177 * StdValue_Lu176_Hf177 / StdSpline_Lu176_Hf177
    data.createTimeSeries("StdCorr_Lu176_Hf177", data.Output, indexChannel.time(), StdCorr_Lu176_Hf177)

    StdSpline_Hf176_177 = data.spline(rmName, "Hf176_177").data()
    StdCorr_Hf176_177 = Hf176_177 * StdValue_Hf176_177 / StdSpline_Hf176_177
    data.createTimeSeries("StdCorr_Hf176_177", data.Output, indexChannel.time(), StdCorr_Hf176_177)

    StdSpline_Lu176_Hf176 = data.spline(rmName, "Lu176_Hf176").data()
    StdCorr_Lu176_Hf176 = Lu176_Hf176 * StdValue_Lu176_Hf176 / StdSpline_Lu176_Hf176
    data.createTimeSeries("StdCorr_Lu176_Hf176", data.Output, indexChannel.time(), StdCorr_Lu176_Hf176)

    StdSpline_Hf177_176 = data.spline(rmName, "Hf177_176").data()
    StdCorr_Hf177_176 = Hf177_176 * StdValue_Hf177_176 / StdSpline_Hf177_176
    data.createTimeSeries("StdCorr_Hf177_176", data.Output, indexChannel.time(), StdCorr_Hf177_176)

    if propErrors:
        groups = [s for s in data.selectionGroupList() if s.type != data.Baseline]
        data.propagateErrors(groups, [data.timeSeries("StdCorr_Hf176_177")], data.timeSeries("Hf176_177"), rmName)

    # Register error correlations
    data.registerAssociatedResult('rho Hf176_177 v LuHf', calcRho1)
    data.registerAssociatedResult('rho Hf177_176 v LuHf', calcRho2)

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

    drs.setSetting("IndexChannel", "Lu175")
    drs.setSetting("ReferenceMaterial", "G_NIST610")
    drs.setSetting("MaskChannel", "Hf178")
    drs.setSetting("MaskCutoff", 1.0)
    drs.setSetting("MaskTrim", 0.0)
    drs.setSetting("Lu176_175ref", 0.02659)
    drs.setSetting("Yb176_172ref", 0.5845)
    drs.setSetting("Hf177_178ref", 0.682)
    drs.setSetting("Hf176_178ref", 0.192)
    drs.setSetting("PropagateError", False)

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

    Lu176_175refLineEdit = QtGui.QLineEdit(widget)
    Lu176_175refLineEdit.setText(settings["Lu176_175ref"])
    Lu176_175refLineEdit.textChanged.connect(lambda t: drs.setSetting("Lu176_175ref", float(t)))
    formLayout.addRow("Lu176_175ref", Lu176_175refLineEdit)

    Yb176_172refLineEdit = QtGui.QLineEdit(widget)
    Yb176_172refLineEdit.setText(settings["Yb176_172ref"])
    Yb176_172refLineEdit.textChanged.connect(lambda t: drs.setSetting("Yb176_172ref", float(t)))
    formLayout.addRow("Yb176_172ref", Yb176_172refLineEdit)

    Hf177_178refLineEdit = QtGui.QLineEdit(widget)
    Hf177_178refLineEdit.setText(settings["Hf177_178ref"])
    Hf177_178refLineEdit.textChanged.connect(lambda t: drs.setSetting("Hf177_178ref", float(t)))
    formLayout.addRow("Hf177_178ref", Hf177_178refLineEdit)

    Hf176_178refLineEdit = QtGui.QLineEdit(widget)
    Hf176_178refLineEdit.setText(settings["Hf176_178ref"])
    Hf176_178refLineEdit.textChanged.connect(lambda t: drs.setSetting("Hf176_178ref", float(t)))
    formLayout.addRow("Hf176_178ref", Hf176_178refLineEdit)

    propCheckBox = QtGui.QCheckBox(widget)
    propCheckBox.setChecked(settings["PropagateError"])
    propCheckBox.toggled.connect(lambda t: drs.setSetting("PropagateError", bool(t)))
    formLayout.addRow("PropagateError", propCheckBox)

    drs.setSettingsWidget(widget)
