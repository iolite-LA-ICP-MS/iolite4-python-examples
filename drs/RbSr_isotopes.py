# / Type: DRS
# / Name: Rb-Sr Isotopes
# / Authors: iolite Software
# / Description: A basic Rb-Sr isotope DRS, based on that described in Redaa et al. J. Anal. At. Spectrom., 2021, 36, 322
# / References: https://doi.org/10.1039/D0JA00299B
# / Version: 0.1
# / Contact: support@iolite-software.com

from iolite import QtGui
import numpy as np


def runDRS():

    drs.message("Starting Rb-Sr isotopes DRS...")
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

    # Create debug messages for the settings being used
    IoLog.debug("indexChannelName = %s" % indexChannel.name)
    IoLog.debug(
        "Masking data  = True" if maskOption else "Masking data  = False")
    IoLog.debug("maskChannelName = %s" % maskChannel.name)
    IoLog.debug("maskCutoff = %f" % cutoff)
    IoLog.debug("maskTrim = %f" % trim)

    # Setup index time
    drs.message("Setting up index time...")
    drs.progress(5)
    drs.setIndexChannel(indexChannel)

    # Setup the mask
    if maskOption:
        drs.message("Making mask...")
        drs.progress(10)
        mask = drs.createMaskFromCutoff(maskChannel, cutoff, trim)
        data.createTimeSeries('mask', data.Intermediate,
                              indexChannel.time(), mask)
    else:
        mask = np.ones_like(indexChannel.data())
        data.createTimeSeries('mask', data.Intermediate,
                              indexChannel.time(), mask)

    # Interp onto index time and baseline subtract
    drs.message("Interpolating onto index time and baseline subtracting...")
    drs.progress(25)

    allInputChannels = data.timeSeriesList(data.Input)
    blGrp = None

    if len(data.selectionGroupList(data.Baseline)) > 1:
        IoLog.error("There are multiple baseline groups. Rb-Sr DRS cannot proceed...")
        drs.message("Error. See Messages")
        drs.progress(100)
        drs.finished()
        return
    elif len(data.selectionGroupList(data.Baseline)) < 1:
        IoLog.error("No baselines. Please select some baselines. Rb-Sr DRS cannot proceed...")
        drs.message("Error. See Messages")
        drs.progress(100)
        drs.finished()
        return
    else:
        blGrp = data.selectionGroupList(data.Baseline)[0]

    if len(blGrp.selections()) < 1:
        IoLog.error("No baseline selections. Please select some baselines. Rb-Sr DRS cannot proceed...")
        drs.message("Error. See Messages")
        drs.progress(100)
        drs.finished()
        return

    for counter, channel in enumerate(allInputChannels):
        drs.message("Baseline subtracting %s" % channel.name)
        drs.progress(25 + 50*counter/len(allInputChannels))

        drs.baselineSubtract(blGrp, [allInputChannels[counter]], mask, 25, 75)
        cps_ch = data.timeSeries(channel.name + '_CPS')
        input_ch = data.timeSeries(channel.name)
        cps_ch.setProperty(
            'Pre-shift mass', input_ch.property('Pre-shift mass'))
        cps_ch.setProperty('Post-shift mass',
                           input_ch.property('Post-shift mass'))

    drs.message("Calculating raw ratios...")
    drs.progress(50)

    # Declare the channels used in the calculations:
    Rb85_CPS = data.timeSeriesList(
        data.Intermediate, {'Post-shift mass': '85'})[0].data()
    Sr86_102_CPS = data.timeSeriesList(
        data.Intermediate, {'Post-shift mass': '102'})[0].data()
    Sr87_103_CPS = data.timeSeriesList(
        data.Intermediate, {'Post-shift mass': '103'})[0].data()
    Sr88_104_CPS = data.timeSeriesList(
        data.Intermediate, {'Post-shift mass': '104'})[0].data()

    Rb87_CPS = Rb85_CPS * 0.38562
    Rb85_Sr86s_Raw = Rb85_CPS/Sr86_102_CPS
    Sr87s_Sr86s_Raw = Sr87_103_CPS/Sr86_102_CPS
    Rb87_Sr86s_Raw = Rb87_CPS/Sr86_102_CPS
    Sr87s_Rb87_Raw = Sr87_103_CPS/Rb87_CPS
    Sr88s_Sr86s_Raw = Sr88_104_CPS/Sr86_102_CPS

    # Gather up intermediate channels and add them as time series:
    int_channel_names = ['Rb87_CPS', 'Rb85_Sr86s_Raw', 'Sr87s_Sr86s_Raw',
                         'Rb87_Sr86s_Raw', 'Sr87s_Rb87_Raw', 'Sr88s_Sr86s_Raw']
    int_channels = [Rb87_CPS, Rb85_Sr86s_Raw, Sr87s_Sr86s_Raw,
                    Rb87_Sr86s_Raw, Sr87s_Rb87_Raw, Sr88s_Sr86s_Raw]
    for name, channel in zip(int_channel_names, int_channels):
        data.createTimeSeries(name, data.Intermediate,
                              indexChannel.time(), channel)

    drs.message("Correcting ratios...")
    drs.progress(80)

    StdSpline_Rb87_Sr86s = data.spline(rmName, "Rb87_Sr86s_Raw").data()
    try:
        StdValue_Rb87_Sr86 = data.referenceMaterialData(rmName)["87Rb/86Sr"].value()
    except KeyError:
        IoLog.error("There was no 87Rb/86Sr value in the " + rmName +
                      " datafile. Rb-Sr DRS cannot proceed.")
        drs.message("Error. See Messages")
        drs.progress(100)
        drs.finished()
        return

    print("StdSpline_Rb87_Sr86 mean = %f" % StdSpline_Rb87_Sr86s.mean())
    print("StdValue_Rb87_Sr86 = %f" % StdValue_Rb87_Sr86)

    StdCorr_Rb87_Sr86s = (Rb87_Sr86s_Raw) * StdValue_Rb87_Sr86 / StdSpline_Rb87_Sr86s
    data.createTimeSeries('StdCorr_Rb87_Sr86s', data.Output,
                          indexChannel.time(), StdCorr_Rb87_Sr86s)

    StdSpline_Sr87s_Sr86s = data.spline(rmName, "Sr87s_Sr86s_Raw").data()
    try:
        StdValue_Sr87_Sr86 = data.referenceMaterialData(rmName)["87Sr/86Sr"].value()
    except KeyError:
        IoLog.error("There was no 87Sr/86Sr value in the " + rmName +
                      " datafile. Rb-Sr DRS cannot proceed.")
        drs.message("Error. See Messages")
        drs.progress(100)
        drs.finished()
        return

    print("StdSpline_Sr87s_Sr86s mean = %f" % StdSpline_Sr87s_Sr86s.mean())
    print("StdValue_Sr87_Sr86 = %f" % StdValue_Sr87_Sr86)

    StdCorr_Sr87s_Sr86s = (Sr87s_Sr86s_Raw) * StdValue_Sr87_Sr86 / StdSpline_Sr87s_Sr86s
    data.createTimeSeries('StdCorr_Sr87s_Sr86s', data.Output,
                          indexChannel.time(), StdCorr_Sr87s_Sr86s)

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
    drs.setSetting("ReferenceMaterial", "A_MAD")
    drs.setSetting("Mask", False)
    drs.setSetting("MaskChannel", defaultChannelName)
    drs.setSetting("MaskCutoff", 0.1)
    drs.setSetting("MaskTrim", 0.0)

    settings = drs.settings()

    indexComboBox = QtGui.QComboBox(widget)
    indexComboBox.addItems(timeSeriesNames)
    indexComboBox.setCurrentText(settings["IndexChannel"])
    indexComboBox.currentTextChanged.connect(
        lambda t: drs.setSetting("IndexChannel", t))
    formLayout.addRow("Index channel", indexComboBox)

    rmComboBox = QtGui.QComboBox(widget)
    rmComboBox.addItems(rmNames)
    rmComboBox.setCurrentText(settings["ReferenceMaterial"])
    rmComboBox.currentTextChanged.connect(
        lambda t: drs.setSetting("ReferenceMaterial", t))
    formLayout.addRow("Reference material", rmComboBox)

    verticalSpacer = QtGui.QSpacerItem(
        20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
    formLayout.addItem(verticalSpacer)

    maskCheckBox = QtGui.QCheckBox(widget)
    maskCheckBox.setChecked(settings["Mask"])
    maskCheckBox.toggled.connect(lambda t: drs.setSetting("Mask", bool(t)))
    formLayout.addRow("Mask", maskCheckBox)

    maskComboBox = QtGui.QComboBox(widget)
    maskComboBox.addItems(data.timeSeriesNames(data.Input))
    maskComboBox.setCurrentText(settings["MaskChannel"])
    maskComboBox.currentTextChanged.connect(
        lambda t: drs.setSetting("MaskChannel", t))
    formLayout.addRow("Mask channel", maskComboBox)

    maskLineEdit = QtGui.QLineEdit(widget)
    maskLineEdit.setText(settings["MaskCutoff"])
    maskLineEdit.textChanged.connect(
        lambda t: drs.setSetting("MaskCutoff", float(t)))
    formLayout.addRow("Mask cutoff", maskLineEdit)

    maskTrimLineEdit = QtGui.QLineEdit(widget)
    maskTrimLineEdit.setText(settings["MaskTrim"])
    maskTrimLineEdit.textChanged.connect(
        lambda t: drs.setSetting("MaskTrim", float(t)))
    formLayout.addRow("Mask trim", maskTrimLineEdit)

    verticalSpacer2 = QtGui.QSpacerItem(
        20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
    formLayout.addItem(verticalSpacer2)

    drs.setSettingsWidget(widget)
