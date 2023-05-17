# A python-based data reduction scheme for iolite 4 starts with some metadata
#/ Type: DRS
#/ Name: Pb_Isotopes_DRS_TEST
#/ Authors: BP
#/ Description: Test For Python Course
#/ References:
#/ Version:
#/ Contact:

from iolite import QtGui
from time import sleep
import numpy as np


def runDRS():
    drs.message("Starting Pb isotope TEST")
    drs.progress(0)

    # Get settings
    settings = drs.settings()
    print(settings)

    indexChannel = data.timeSeries(settings["IndexChannel"])
    maskChannel = data.timeSeries(settings["MaskChannel"])
    cutoff = settings["MaskCutoff"]
    trim = settings["MaskTrim"]
    rmName = settings["Primary RM"]

    # Create debug messages for the settings being used
    IoLog.debug("indexChannelName = %s" % indexChannel.name)
    IoLog.debug("maskChannelName = %s" % maskChannel.name)
    IoLog.debug("maskCutoff = %f" % cutoff)
    IoLog.debug("maskTrim = %f" % trim)
    IoLog.debug("Primary RM = %s" % rmName)

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

    drs.baselineSubtract(data.selectionGroup("Baselines"), allInputChannels, mask, 25, 35)

    # Calculate raw ratios:
    Pb208_CPS = data.timeSeriesByMass(data.Intermediate, 208., 0.1).data()
    Pb206_CPS = data.timeSeriesByMass(data.Intermediate, 206., 0.1).data()
    Pb86_raw = Pb208_CPS / Pb206_CPS

    try:
        StdValue_208_206 = data.referenceMaterialData(rmName)["208Pb/206Pb"].value()
    except Exception as e:
        IoLog.errorWithOrigin(f"Could not get the 208Pb/206Pb value from the RM file: {e}", "BP_DRS")
        drs.message("DRS did not finish. Please check Messages")
        drs.progress(100)
        drs.finished()
        return

    Pb_f = (np.log(StdValue_208_206 / (Pb208_CPS / Pb206_CPS))) / (np.log(207.976627 / 205.974440))

    data.createTimeSeries("Pb208/Pb206 raw", data.Intermediate, indexChannel.time(), Pb86_raw)
    data.createTimeSeries("Pb_f", data.Intermediate, indexChannel.time(), Pb_f)

    data.updateResults()
    StdSpline_86 = data.spline(rmName, "Pb_f").data()

    Pb86_corr = Pb86_raw * np.power((207.976627/205.974440), StdSpline_86)
    data.createTimeSeries("Pb86_corr", data.Output, indexChannel.time(), Pb86_corr)

    drs.message("Finished!")
    drs.progress(100)
    drs.finished()


def settingsWidget():
    widget = QtGui.QWidget()
    formLayout = QtGui.QFormLayout()
    widget.setLayout(formLayout)

    timeSeriesNames = data.timeSeriesNames(data.Input)
    defaultChannelName = ""
    if timeSeriesNames:
        defaultChannelName = timeSeriesNames[0]

    drs.setSetting("IndexChannel", defaultChannelName)
    drs.setSetting("MaskChannel", defaultChannelName)
    drs.setSetting("MaskCutoff", 100000.0)
    drs.setSetting("MaskTrim", 0.0)
    drs.setSetting("Primary RM", "G_NIST612")

    settings = drs.settings()

    indexComboBox = QtGui.QComboBox(widget)
    indexComboBox.addItems(data.timeSeriesNames(data.Input))
    indexComboBox.setCurrentText(settings["IndexChannel"])
    indexComboBox.textActivated.connect(lambda t: drs.setSetting("IndexChannel", t))

    maskComboBox = QtGui.QComboBox(widget)
    maskComboBox.addItems(data.timeSeriesNames(data.Input))
    maskComboBox.setCurrentText(settings["MaskChannel"])
    maskComboBox.textActivated.connect(lambda t: drs.setSetting("MaskChannel", t))

    maskLineEdit = QtGui.QLineEdit(widget)
    maskLineEdit.setText(settings["MaskCutoff"])
    maskLineEdit.textChanged.connect(lambda t: drs.setSetting("MaskCutoff", float(t)))

    maskTrimLineEdit = QtGui.QLineEdit(widget)
    maskTrimLineEdit.setText(settings["MaskTrim"])
    maskTrimLineEdit.textChanged.connect(lambda t: drs.setSetting("MaskTrim", float(t)))

    rmComboBox = QtGui.QComboBox(widget)
    rmComboBox.addItems(data.referenceMaterialNames())
    rmComboBox.setCurrentText(settings["Primary RM"])
    rmComboBox.textActivated.connect(lambda t: drs.setSetting("Primary RM", t))

    formLayout.addRow("Index channel", indexComboBox)
    formLayout.addRow("Mask channel", maskComboBox)
    formLayout.addRow("Mask cutoff", maskLineEdit)
    formLayout.addRow("Mask trim", maskTrimLineEdit)
    formLayout.addRow("Primary RM", rmComboBox)

    drs.setSettingsWidget(widget)
