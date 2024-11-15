# A python-based data reduction scheme for iolite 4 starts with some metadata
#/ Type: DRS
#/ Name: Pb Isotopes DRS (Beta)
#/ Authors: BP
#/ Description: A basic Pb Isotope DRS for iolite 4
#/ References:
#/ Version:
#/ Contact:

from iolite.QtGui import QWidget, QSpacerItem, QSizePolicy, QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox, QComboBox, QLineEdit
from iolite.QtGui import QDoubleSpinBox
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
    if settings["UseMask"]:
        drs.message("Making mask...")
        drs.progress(10)
        mask = drs.createMaskFromCutoff(maskChannel, cutoff, trim)
    
    else:
        mask = np.ones_like(maskChannel.data())

    # Interp onto index time and baseline subtract
    drs.message("Interpolating onto index time and baseline subtracting...")
    drs.progress(25)

    allInputChannels = data.timeSeriesList(data.Input)

    drs.baselineSubtract(data.selectionGroup("Baselines"), allInputChannels, mask, 25, 35)

    drs.message("Calculating raw ratios...")
    drs.progress(35)

    # Calculate raw ratios:
    Pb208_CPS = data.timeSeriesByMass(data.Intermediate, 208., 0.1).data()
    Pb207_CPS = data.timeSeriesByMass(data.Intermediate, 207., 0.1).data()
    Pb206_CPS = data.timeSeriesByMass(data.Intermediate, 206., 0.1).data()
    m204_CPS = data.timeSeriesByMass(data.Intermediate, 204., 0.1).data()
    Hg202_CPS = data.timeSeriesByMass(data.Intermediate, 202., 0.1).data()

    Pb86_raw = Pb208_CPS / Pb206_CPS * mask
    Pb76_raw = Pb207_CPS / Pb206_CPS * mask

    simpleHg204 = Hg202_CPS * 0.2301
    data.createTimeSeries("simpleHg204", data.Intermediate, indexChannel.time(), simpleHg204)

    Pb204_raw = m204_CPS - Hg202_CPS * 0.2299
    Pb84_raw = Pb208_CPS / Pb204_raw * mask
    Pb74_raw = Pb207_CPS / Pb204_raw * mask
    Pb64_raw = Pb206_CPS / Pb204_raw * mask

    data.createTimeSeries("Pb204_raw", data.Intermediate, indexChannel.time(), Pb204_raw)
    data.createTimeSeries("raw 208Pb/204Pb", data.Intermediate, indexChannel.time(), Pb84_raw)
    data.createTimeSeries("raw 207Pb/204Pb", data.Intermediate, indexChannel.time(), Pb74_raw)
    data.createTimeSeries("raw 206Pb/204Pb", data.Intermediate, indexChannel.time(), Pb64_raw)
    data.createTimeSeries("raw 208Pb/206Pb", data.Intermediate, indexChannel.time(), Pb86_raw)
    data.createTimeSeries("raw 207Pb/206Pb", data.Intermediate, indexChannel.time(), Pb76_raw)

    drs.message("Correcting for mass fractionation...")
    drs.progress(75)

    try:
        StdValue_208_206 = data.referenceMaterialData(rmName)["208Pb/206Pb"].value()
    except Exception as e:
        IoLog.errorWithOrigin(f"Could not get the 208Pb/206Pb value from the RM file: {e}", "Pb Isotope DRS")
        drs.message("DRS did not finish. Please check Messages")
        drs.progress(100)
        drs.finished()
        return

    Pb_f = (np.log(StdValue_208_206 / (Pb208_CPS / Pb206_CPS))) / (np.log(207.976627 / 205.974440))

    data.createTimeSeries("Pb_F", data.Intermediate, indexChannel.time(), Pb_f)

    data.updateResults()
    StdSpline_PbF = data.spline(rmName, "Pb_F").data()

    if settings["UseHgF"]:
        hg_factor = float(settings["HgFactor"])
        print(f'Using Hg Factor: {hg_factor}')
    else:
        hg_factor = 1.0

    HgF = StdSpline_PbF * hg_factor
    Hg204 = Hg202_CPS * 0.2301 / np.power((203.973481/201.970632), HgF)
    Pb204_corr = m204_CPS - Hg204
    data.createTimeSeries("Hg204", data.Intermediate, indexChannel.time(), Hg204)
    data.createTimeSeries("Pb204_corr", data.Intermediate, indexChannel.time(), Pb204_corr)

    Pb86_corr = Pb86_raw * np.power((207.976627/205.97446), StdSpline_PbF)
    Pb76_corr = Pb76_raw * np.power((206.97589/205.97446), StdSpline_PbF)
    Pb84_corr = (Pb208_CPS / Pb204_corr) * np.power((207.976627/203.97304), StdSpline_PbF)
    Pb74_corr = (Pb207_CPS / Pb204_corr) * np.power((206.97589/203.97304), StdSpline_PbF)
    Pb64_corr = (Pb206_CPS / Pb204_corr) * np.power((205.97446/203.97304), StdSpline_PbF)

    data.createTimeSeries("Final 208Pb/206Pb", data.Output, indexChannel.time(), Pb86_corr)
    data.createTimeSeries("Final 207Pb/206Pb", data.Output, indexChannel.time(), Pb76_corr)
    data.createTimeSeries("Final 208Pb/204Pb", data.Output, indexChannel.time(), Pb84_corr)
    data.createTimeSeries("Final 207Pb/204Pb", data.Output, indexChannel.time(), Pb74_corr)
    data.createTimeSeries("Final 206Pb/204Pb", data.Output, indexChannel.time(), Pb64_corr)

    # Normalise 204 ratios to the standard
    drs.message("Normalising 204 ratios to the standard...")
    drs.progress(90)

    try:
        StdValue_206_204 = data.referenceMaterialData(rmName)["206Pb/204Pb"].value()
        StdValue_207_204 = data.referenceMaterialData(rmName)["207Pb/204Pb"].value()
        StdValue_208_204 = data.referenceMaterialData(rmName)["208Pb/204Pb"].value()
    except Exception as e:
        IoLog.errorWithOrigin(f"Could not get one of the 20XPb/204Pb values from the RM file: {e}", "Pb Isotope DRS")
        drs.message("DRS did not finish. Please check Messages")
        drs.progress(100)
        drs.finished()
        return

    Pb64_corrFact = Pb64_corr / StdValue_206_204
    Pb74_corrFact = Pb74_corr / StdValue_207_204
    Pb84_corrFact = Pb84_corr / StdValue_208_204

    data.createTimeSeries("206Pb/204Pb_corrFact", data.Intermediate, indexChannel.time(), Pb64_corrFact)
    data.createTimeSeries("207Pb/204Pb_corrFact", data.Intermediate, indexChannel.time(), Pb74_corrFact)
    data.createTimeSeries("208Pb/204Pb_corrFact", data.Intermediate, indexChannel.time(), Pb84_corrFact)

    Pb64_corrSpline = data.spline(rmName, "206Pb/204Pb_corrFact").data()
    Pb74_corrSpline = data.spline(rmName, "207Pb/204Pb_corrFact").data()
    Pb84_corrSpline = data.spline(rmName, "208Pb/204Pb_corrFact").data()

    Pb64_norm = Pb64_corr / Pb64_corrSpline
    Pb74_norm = Pb74_corr / Pb74_corrSpline
    Pb84_norm = Pb84_corr / Pb84_corrSpline

    data.createTimeSeries("206Pb/204Pb_norm", data.Output, indexChannel.time(), Pb64_norm)
    data.createTimeSeries("207Pb/204Pb_norm", data.Output, indexChannel.time(), Pb74_norm)
    data.createTimeSeries("208Pb/204Pb_norm", data.Output, indexChannel.time(), Pb84_norm)

    # Finally, trying a separate PbF for ratio 204Pb ratio
    drs.message("Calculating PbF for 204Pb ratios...")
    drs.progress(95)

    PbF_64 = (np.log(StdValue_206_204 / (Pb206_CPS / Pb204_corr))) / (np.log(205.974440 / 203.97304))
    PbF_74 = (np.log(StdValue_207_204 / (Pb207_CPS / Pb204_corr))) / (np.log(206.97589 / 203.97304))
    PbF_84 = (np.log(StdValue_208_204 / (Pb208_CPS / Pb204_corr))) / (np.log(207.976627 / 203.97304))

    data.createTimeSeries("PbF_64", data.Intermediate, indexChannel.time(), PbF_64)
    data.createTimeSeries("PbF_74", data.Intermediate, indexChannel.time(), PbF_74)
    data.createTimeSeries("PbF_84", data.Intermediate, indexChannel.time(), PbF_84)

    PbF_64_std_spline = data.spline(rmName, "PbF_64").data()
    PbF_74_std_spline = data.spline(rmName, "PbF_74").data()
    PbF_84_std_spline = data.spline(rmName, "PbF_84").data()

    Pb64_w64PbF = (Pb206_CPS / Pb204_corr) * np.power((205.97446/203.97304), PbF_64_std_spline)
    Pb74_w74PbF = (Pb207_CPS / Pb204_corr) * np.power((206.97589/203.97304), PbF_74_std_spline)
    Pb84_w84PbF = (Pb208_CPS / Pb204_corr) * np.power((207.976627/203.97304), PbF_84_std_spline)

    data.createTimeSeries("Final 206Pb/204Pb w64PbF", data.Output, indexChannel.time(), Pb64_w64PbF)
    data.createTimeSeries("Final 207Pb/204Pb w74PbF", data.Output, indexChannel.time(), Pb74_w74PbF)
    data.createTimeSeries("Final 208Pb/204Pb w84PbF", data.Output, indexChannel.time(), Pb84_w84PbF)

    drs.message("Finished!")
    drs.progress(100)
    drs.finished()


def settingsWidget():
    widget = QWidget()
    hspacer_left = QSpacerItem(20, 40, QSizePolicy.Expanding, QSizePolicy.Expanding)
    hspacer_right = QSpacerItem(20, 40, QSizePolicy.Expanding, QSizePolicy.Expanding)
    mainLayout = QVBoxLayout()
    hLayout = QHBoxLayout()
    hLayout.addItem(hspacer_left)
    hLayout.addLayout(mainLayout)
    hLayout.addItem(hspacer_right)
    widget.setLayout(hLayout)

    timeSeriesNames = data.timeSeriesNames(data.Input)
    defaultChannelName = ""
    if timeSeriesNames:
        defaultChannelName = timeSeriesNames[0]

    drs.setSetting("Primary RM", "G_NIST612")
    drs.setSetting("IndexChannel", defaultChannelName)
    drs.setSetting("UseMask", True)
    drs.setSetting("MaskChannel", defaultChannelName)
    drs.setSetting("MaskCutoff", 100000.0)
    drs.setSetting("MaskTrim", 0.0)
    drs.setSetting("UseHgF", True)
    drs.setSetting("HgFactor", 1.0)

    settings = drs.settings()

    basicsGroupBox = QGroupBox(widget)
    basicsGroupBox.setTitle("Basic settings")
    basicsLayout = QFormLayout()

    rmComboBox = QComboBox(widget)
    rmComboBox.addItems(data.referenceMaterialNames())
    rmComboBox.setCurrentText(settings["Primary RM"])
    rmComboBox.textActivated.connect(lambda t: drs.setSetting("Primary RM", t))
    basicsLayout.addRow("Primary RM", rmComboBox)

    indexComboBox = QComboBox(widget)
    indexComboBox.addItems(data.timeSeriesNames(data.Input))
    indexComboBox.setCurrentText(settings["IndexChannel"])
    indexComboBox.textActivated.connect(lambda t: drs.setSetting("IndexChannel", t))
    basicsLayout.addRow("Index channel", indexComboBox)

    basicsGroupBox.setLayout(basicsLayout)

    maskGroupBox = QGroupBox(widget)
    maskGroupBox.setTitle("Mask settings")
    maskGroupBox.setCheckable(True)
    maskGroupBox.setChecked(True)
    maskGroupBox.toggled.connect(lambda t: drs.setSetting("UseMask", t))
    maskLayout = QFormLayout()

    maskComboBox = QComboBox(widget)
    maskComboBox.addItems(data.timeSeriesNames(data.Input))
    maskComboBox.setCurrentText(settings["MaskChannel"])
    maskComboBox.textActivated.connect(lambda t: drs.setSetting("MaskChannel", t))

    maskLineEdit = QLineEdit(widget)
    maskLineEdit.setText(settings["MaskCutoff"])
    maskLineEdit.textChanged.connect(lambda t: drs.setSetting("MaskCutoff", float(t)))

    maskTrimLineEdit = QLineEdit(widget)
    maskTrimLineEdit.setText(settings["MaskTrim"])
    maskTrimLineEdit.textChanged.connect(lambda t: drs.setSetting("MaskTrim", float(t)))

    maskLayout.addRow("Mask channel", maskComboBox)
    maskLayout.addRow("Mask cutoff", maskLineEdit)
    maskLayout.addRow("Mask trim", maskTrimLineEdit)
    maskGroupBox.setLayout(maskLayout)

    HgFGroupBox = QGroupBox(widget)
    HgFGroupBox.setTitle("HgF settings")
    HgFGroupBox.setCheckable(True)
    HgFGroupBox.setChecked(True)
    HgFGroupBox.toggled.connect(lambda t: drs.setSetting("UseHgF", t))
    HgFLayout = QFormLayout()
    
    HgFSpinBox = QDoubleSpinBox(widget)
    HgFSpinBox.setValue(settings["HgFactor"])
    HgFSpinBox.setRange(0.0, 5.0)
    HgFSpinBox.setSingleStep(0.1)
    HgFSpinBox.setDecimals(3)
    HgFSpinBox.valueChanged.connect(lambda t: drs.setSetting("HgFactor", t))
    HgFLayout.addRow("HgF Factor", HgFSpinBox)
    HgFGroupBox.setLayout(HgFLayout)

    vSpacer = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)

    mainLayout.addWidget(basicsGroupBox)
    mainLayout.addWidget(maskGroupBox)
    mainLayout.addWidget(HgFGroupBox)
    mainLayout.addItem(vSpacer)

    drs.setSettingsWidget(widget)
