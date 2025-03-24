#/ Type: DRS
#/ Name: Ni Isotope Double Spike
#/ Authors: Bence Paul
#/ Description: A Ni isotope double-spike DRS, based on the Pt_DS DRS by John Creech
#/ References: isospike.org, https://doi.org/10.1111/j.1751-908X.2014.00276.x
#/ Version: 0.1
#/ Contact: support@iolite-software.com


from iolite.QtGui import QComboBox, QWidget, QFormLayout, QSpacerItem, QSizePolicy
from iolite.QtGui import QCheckBox, QLineEdit
import numpy as np
import pandas as pd
from IsoSpike_iolite4 import IsoSpike

def runDRS():
    drs.message("Starting Ni Double Spike DRS...")
    drs.progress(0)

    # Get DRS settings
    settings = drs.settings()
    indexChannel = data.timeSeries(settings["IndexChannel"])
    rmName = settings["ReferenceMaterial"]

    # Double-spike Inversion Parameters
    rationames = ['Ni60/Ni58', 'Ni61/Ni58', 'Ni62/Ni58']
    unmixedRatios = [0.3852, 0.0167, 0.0534]
    spikeRatios = [0.61283, 81.22346, 85.94008]
    logMassRatios = [0.03386, 0.05042, 0.06665]

    DSsettings = np.array([unmixedRatios, spikeRatios, logMassRatios])

    # Create debug messages for the settings being used
    IoLog.debug("indexChannelName = %s" % indexChannel.name)

    # Setup index time
    drs.message("Setting up index time...")
    drs.progress(5)

    # Interp onto index time and baseline subtract
    drs.message("Interpolating onto index time and baseline subtracting...")
    drs.progress(25)
    #
    # Setup the mask
    maskOption = settings["Mask"]
    maskChannel = data.timeSeries(settings["MaskChannel"])
    cutoff = settings["MaskCutoff"]
    trim = settings["MaskTrim"]
    if maskOption:
        drs.message("Making mask...")
        drs.progress(10)
        mask = drs.createMaskFromCutoff(maskChannel, cutoff, trim)
        data.createTimeSeries('mask', data.Intermediate, indexChannel.time(), mask)
    else:
        mask = np.ones_like(indexChannel.data())
        data.createTimeSeries('mask', data.Intermediate, indexChannel.time(), mask)


    # Do baseline subtraction
    allInputChannels = data.timeSeriesList(data.Input)
    if len(data.selectionGroupList(data.Baseline)) < 1:
        IoLog.error("Could not find a baseline group. Combined Sr DRS cannot proceed...")
        drs.message("Please select some baseline selections")
        drs.progress(100)
        drs.finished()
        return

    if len(data.selectionGroupList(data.Baseline)) > 1:
        IoLog.error("There are more than one baseline groups. Combined Sr DRS cannot proceed...")
        drs.message("DRS did not finish. Please check Messages")
        drs.progress(100)
        drs.finished()
        return
    elif len(data.selectionGroupList(data.Baseline)) < 1:
        IoLog.error("No baseline groups found. Have you selected baselines yet?")
        drs.message("DRS did not finish. Please check Messages")
        drs.progress(100)
        drs.finished()
        return
    else:
        blGrp = data.selectionGroupList(data.Baseline)[0]


    drs.baselineSubtract(blGrp, allInputChannels, mask, 15, 35)

    drs.message("Calculating raw ratios")
    drs.progress(50)

    # Reference baseline subtracted data
    Ni58_CPS = data.timeSeries("Ni58_CPS").data()
    Ni60_CPS = data.timeSeries("Ni60_CPS").data()
    Ni61_CPS = data.timeSeries("Ni61_CPS").data()
    Ni62_CPS = data.timeSeries("Ni62_CPS").data()
    x64_CPS = data.timeSeries("Ni64_CPS").data()  # This will be corrected for 64Zn interference below...

    # Correct for 64Zn interference
    Zn66_CPS = data.timeSeries("Zn66_CPS").data()
    Zn64 = 1.773 * Zn66_CPS
    Ni64_CPS = x64_CPS - Zn64

    # Calculate raw ratios
    Raw60_58 = Ni60_CPS / Ni58_CPS
    Raw61_58 = Ni61_CPS / Ni58_CPS
    Raw62_58 = Ni62_CPS / Ni58_CPS

    # Add raw ratios to intermediate channels
    data.createTimeSeries("Raw60_58", data.Intermediate, indexChannel.time(), Raw60_58)
    data.createTimeSeries("Raw61_58", data.Intermediate, indexChannel.time(), Raw61_58)
    data.createTimeSeries("Raw62_58", data.Intermediate, indexChannel.time(), Raw62_58)

    # Calculate total Ni volts
    Total_Ni_Volts = Ni58_CPS + Ni60_CPS + Ni61_CPS + Ni62_CPS + Ni64_CPS
    data.createTimeSeries("Total Ni Volts", data.Output, indexChannel.time(), Total_Ni_Volts)

    ####################################################################################################
    ### Everything in the DRS that needs to be customised is above this line; below is all abstracted.
    ####################################################################################################

    drs.message("Calling IsoSpike")
    drs.progress(70)
    #######
    #### Here is where we call IsoSpike, passing on the DS settings and the measured raw ratios.
    #######
    IsoSpike_results = IsoSpike(DSsettings, Raw60_58, Raw61_58, Raw62_58)

    drs.message("Calculating internally normalised values")
    drs.progress(85)

    #### Below is some post-processing stuff.
    ### We need to reference the results of the DS calculations so we can see it in Iolite.

    ones = np.ones(len(indexChannel.time())) # initialise an array of ones; some data didn't work if I just copy into arrays, but fine when multiplying by array of ones

    # Unpack results
    lam = ones*IsoSpike_results[:, 0]
    alpha = ones*IsoSpike_results[:, 1]
    beta = ones*IsoSpike_results[:, 2]
    DScorr_ratio1 = ones*IsoSpike_results[:, 3]
    DScorr_ratio2 = ones*IsoSpike_results[:, 4]
    DScorr_ratio3 = ones*IsoSpike_results[:, 5]
    deltaratio1 = ones*IsoSpike_results[:, 6]
    deltaratio2 = ones*IsoSpike_results[:, 7]
    deltaratio3 = ones*IsoSpike_results[:, 8]

    #create time series in Iolite
    data.createTimeSeries("lambda", data.Output,  indexChannel.time(), lam)
    data.createTimeSeries("alpha", data.Output,  indexChannel.time(), alpha)
    data.createTimeSeries("beta", data.Output,  indexChannel.time(), beta)
    data.createTimeSeries("DScorr"+rationames[0], data.Output,  indexChannel.time(), DScorr_ratio1)
    data.createTimeSeries("DScorr"+rationames[1], data.Output,  indexChannel.time(), DScorr_ratio2)
    data.createTimeSeries("DScorr"+rationames[2], data.Output,  indexChannel.time(), DScorr_ratio3)
    data.createTimeSeries("Delta"+rationames[0], data.Output,  indexChannel.time(), deltaratio1)
    data.createTimeSeries("Delta"+rationames[1], data.Output,  indexChannel.time(), deltaratio2)
    data.createTimeSeries("Delta"+rationames[2], data.Output,  indexChannel.time(), deltaratio3)

    #create splines on the DS corrected ratios for the standard
    try:
        data.selectionGroup(rmName)
    except:
        IoLog.error("The Ni Isotope DS DRS requires Ref Material selections to proceed.")
        drs.message("DRS did not finish. Please check Messages")
        drs.progress(100)
        drs.finished()
        return

    StdSpline1 = data.spline(rmName, 'DScorr'+rationames[0]).data()
    StdSpline2 = data.spline(rmName, 'DScorr'+rationames[1]).data()
    StdSpline3 = data.spline(rmName, 'DScorr'+rationames[2]).data()

    # Calculate internally normalised delta values
    DeltaInt_ratio1 = ((DScorr_ratio1/StdSpline1)-1) * 1000
    DeltaInt_ratio2 = ((DScorr_ratio2/StdSpline2)-1) * 1000
    DeltaInt_ratio3 = ((DScorr_ratio3/StdSpline3)-1) * 1000

    data.createTimeSeries("DeltaInt"+rationames[0], data.Output, indexChannel.time(), DeltaInt_ratio1)
    data.createTimeSeries("DeltaInt"+rationames[1], data.Output, indexChannel.time(), DeltaInt_ratio2)
    data.createTimeSeries("DeltaInt"+rationames[2], data.Output, indexChannel.time(), DeltaInt_ratio3)

    #finishing touches
    drs.message("Finished!")
    drs.progress(100)
    drs.finished()

def settingsWidget():
    widget = QWidget()
    formLayout = QFormLayout()
    widget.setLayout(formLayout)

    defaultChannelName = "Ni58"
    timeSeriesNames = data.timeSeriesNames(data.Input)
    if 'Ni58' in timeSeriesNames:
        defaultChannelName = 'Ni58'
    elif timeSeriesNames:
        defaultChannelName = timeSeriesNames[0]

    drs.setSetting("IndexChannel", defaultChannelName)
    drs.setSetting("ReferenceMaterial", "NIST986")
    drs.setSetting("Mask", True)
    drs.setSetting("MaskChannel", defaultChannelName)
    drs.setSetting("MaskCutoff", 0.05)
    drs.setSetting("MaskTrim", 0.0)

    settings = drs.settings()

    indexComboBox = QComboBox(widget)
    indexComboBox.addItems(timeSeriesNames)
    indexComboBox.setCurrentText(settings["IndexChannel"])
    indexComboBox.textActivated.connect(lambda t: drs.setSetting("IndexChannel", t))
    formLayout.addRow("Index channel", indexComboBox)

    def updateIndexCombo():
        timeSeriesNames = data.timeSeriesNames(data.Input)
        indexComboBox.clear()
        indexComboBox.addItems(timeSeriesNames)

    data.dataChanged.connect(updateIndexCombo)

    rmComboBox = QComboBox(widget)
    rmNames = data.referenceMaterialNames()
    rmComboBox.addItems(rmNames)
    if settings["ReferenceMaterial"] in rmNames:
        rmComboBox.setCurrentText(settings["ReferenceMaterial"])
    else:
        rmComboBox.setCurrentText(rmNames[0])
        drs.setSetting('ReferenceMaterial', rmNames[0])

    rmComboBox.textActivated.connect(lambda t: drs.setSetting("ReferenceMaterial", t))
    formLayout.addRow("Reference material", rmComboBox)

    def updateRMCombo():
        rmNames = data.selectionGroupNames(data.ReferenceMaterial)
        rmComboBox.clear()
        rmComboBox.addItems(rmNames)
        drs.setSetting("ReferenceMaterial", rmComboBox.currentText())

    data.selectionGroupsChanged.connect(updateRMCombo)

    verticalSpacer2 = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Minimum)
    formLayout.addItem(verticalSpacer2)

    # Set up Mask
    maskCheckBox = QCheckBox(widget)
    maskCheckBox.setChecked(settings["Mask"])
    maskCheckBox.toggled.connect(lambda t: drs.setSetting("Mask", bool(t)))
    formLayout.addRow("Mask", maskCheckBox)

    maskComboBox = QComboBox(widget)
    maskComboBox.addItems(data.timeSeriesNames(data.Input))
    maskComboBox.setCurrentText(settings["MaskChannel"])
    maskComboBox.currentTextChanged.connect(lambda t: drs.setSetting("MaskChannel", t))
    formLayout.addRow("Mask channel", maskComboBox)

    maskLineEdit = QLineEdit(widget)
    maskLineEdit.setText(settings["MaskCutoff"])
    maskLineEdit.textChanged.connect(lambda t: drs.setSetting("MaskCutoff", float(t)))
    formLayout.addRow("Mask cutoff", maskLineEdit)

    maskTrimLineEdit = QLineEdit(widget)
    maskTrimLineEdit.setText(settings["MaskTrim"])
    maskTrimLineEdit.textChanged.connect(lambda t: drs.setSetting("MaskTrim", float(t)))
    formLayout.addRow("Mask trim", maskTrimLineEdit)

    verticalSpacer3 = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Minimum)
    formLayout.addItem(verticalSpacer3)

    # Restore settings
    try:
        settings = drs.settings()
        print('Restoring settings...')
        print(settings)
        rmComboBox.setCurrentText(settings["ReferenceMaterial"])
    except KeyError:
        pass

    drs.setSettingsWidget(widget)
