#/ Type: DRS
#/ Name: Sr Isotopes (Combined)
#/ Authors: Bence Paul, Joe Petrus and author(s) of Sr_isotopes_Total_NIGL.ipf
#/ Description: A Sr isotopes DRS that corrects for REE and CaAr interferences
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
    Age = settings["Age"]
    RbBias = settings["RbBias"]
    CaArBias = settings["CaArBias"]
    propErrors = settings["PropagateError"]

    # Create debug messages for the settings being used
    IoLog.debug("indexChannelName = %s" % indexChannel.name)
    IoLog.debug("Masking data  = True" if maskOption else "Masking data  = False")
    IoLog.debug("maskChannelName = %s" % maskChannel.name)
    IoLog.debug("maskCutoff = %f" % cutoff)
    IoLog.debug("maskTrim = %f" % trim)
    IoLog.debug("Age = %f" % Age)
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
        IoLog.error("There are more than one baseline groups. Combined Sr DRS cannot proceed...")
        drs.message("DRS did not finish. Please check Messages")
        drs.progress(100)
        drs.finished()
        return
    else:
        blGrp = data.selectionGroupList(data.Baseline)[0]

    for counter, channel in enumerate(allInputChannels):
        drs.message("Baseline subtracting %s" % channel.name)
        drs.progress(25 + 50*counter/len(allInputChannels))

        drs.baselineSubtract(blGrp, [allInputChannels[counter]], mask, 25, 75)

    drs.message("Checking for half masses...")
    drs.progress(35)

    try:
        Yb86_5 = data.timeSeriesList(data.Intermediate, {'Mass': '86.5'})[0].data()
    except IndexError:
        IoLog.error("Could not find half masses. Combined Sr DRS could not proceed...")
        drs.message("DRS did not finish. Please check Messages")
        drs.progress(100)
        drs.finished()
        return

    #Y89 = data.timeSeriesList(data.Intermediate, {'Mass': '89'})[0].data()
    SrCaArYbLu88 = data.timeSeriesList(data.Intermediate, {'Mass': '88'})[0].data()
    #Lu87_5 = data.timeSeriesList(data.Intermediate, {'Mass': '87.5'})[0].data()
    SrRbYb87 = data.timeSeriesList(data.Intermediate, {'Mass': '87'})[0].data()
    Yb86_5 = data.timeSeriesList(data.Intermediate, {'Mass': '86.5'})[0].data()
    SrCaArYb86 = data.timeSeriesList(data.Intermediate, {'Mass': '86'})[0].data()
    RbYbEr85 = data.timeSeriesList(data.Intermediate, {'Mass': '85'})[0].data()
    #Tm84_5 = data.timeSeriesList(data.Intermediate, {'Mass': '84.5'})[0].data()
    SrCaArErYb84 = data.timeSeriesList(data.Intermediate, {'Mass': '84'})[0].data()
    Er83_5 = data.timeSeriesList(data.Intermediate, {'Mass': '83.5'})[0].data()
    CaArEr83 = data.timeSeriesList(data.Intermediate, {'Mass': '83'})[0].data()
    CaArErDy82 = data.timeSeriesList(data.Intermediate, {'Mass': '82'})[0].data()


    '''Note that this is not baseline subtracted as it is intended to be used as a
    reference of the raw counts
    '''
    Raw83 = data.timeSeriesList(data.Input, {'Mass': '83'})[0].data()

    drs.message("Subtracting interferences...")
    drs.progress(40)

    Sr88_86_reference = settings['Sr88_86_reference']
    Dy_Er = settings['Dy_Er']

    PFract = (np.log(Sr88_86_reference / (SrCaArYbLu88/SrCaArYb86))) / (np.log(87.9056/85.9093))*mask

    '''The following equations subtract CaAr using the signal on 82, itself corrected for REE2+ using the Er83_5 signal and canonical Dy/Er (because we don't measure Dy on a half mass).
    '''

    CaAr82 = CaArErDy82 - (Er83_5 * 1.601 / 22.869) / np.power((81.96460 / 83.46619), PFract) - ((Er83_5 * 1.601 / 22.869) / np.power((81.96460 / 83.46619), PFract) * Dy_Er * 28.260)

    PSrCaArYb86 = SrCaArYb86 - (CaAr82 * .004 / .647) / np.power((85.9160721 / 81.9210049), PFract) - (Yb86_5 * 21.68 / 16.103) / np.power((85.968195 / 86.46911), PFract)

    PSrCaArYbLu88 = SrCaArYbLu88 - ((CaAr82 * .187 / .647) / np.power((87.9149151 / 81.9210049), PFract)) - (Yb86_5 * 12.996 / 16.103) / np.power((87.97129 / 86.46911), PFract) - ((Yb86_5 * 12.996 / 16.103) / np.power((87.97129 / 86.46911), PFract) * Dy_Er * 2.59)

    # Use these CaAr and REE stripped Sr 86 and Sr 88 values to calculate a refined fractionation factor
    BetaSr = (np.log(Sr88_86_reference / (PSrCaArYbLu88 / PSrCaArYb86))) / (np.log(87.9056/85.9093))

    # Calculate Rb fractionation factor, with optional adjustment
    BetaRb = BetaSr + settings['RbBias']

    # Calculate Rb on mass 87
    Rb87_85_reference = settings['Rb87_85_reference']
    Rb87 = ((RbYbEr85-(Yb86_5 * 2.982 / 16.103) / np.power((84.967385 / 86.46911), BetaSr) - (Er83_5 * 14.910 / 22.869) / np.power((84.96774 / 83.46619), BetaSr)) * Rb87_85_reference) / np.power((86.90918 / 84.9118), BetaRb)

    # Subtract this amount of Rb 87 from the 87 beam
    Sr87 = (SrRbYb87 - Rb87 - (Yb86_5 * 32.026 / 16.103) / np.power((86.96943 / 86.46911), BetaSr))

    # Calculate unique CaAr fractionation factor relative to Sr fract.
    BetaCaAr = BetaSr + settings['CaArBias']

    # Calculate the amount of CaAr on mass 84
    CaArErYb84 = (CaAr82 * ((2.086/0.647) / np.power((83.917989 / 81.921122), BetaCaAr))) + (Er83_5 * 26.978 / 22.869) / np.power((83.96619 / 83.46619), BetaSr) + (Yb86_5 * 0.123 / 16.103) / np.power((83.96694 / 86.46911), BetaSr)

    # Subtract this amount from the 84 beam
    Sr84 = (SrCaArErYb84 - CaArErYb84)

    # Calculate a final CaAr and REE corrected 88 beam
    Sr88 = (SrCaArYbLu88 - ((CaAr82 * .187 / .647) / np.power((87.9149151 / 81.9210049), BetaSr)) - (Yb86_5 * 12.996 / 16.103) / np.power((87.97129 / 86.46911), BetaSr) - ((Yb86_5 * 12.996 / 16.103) / np.power((87.97129 / 86.46911), BetaSr) * Dy_Er * 2.59))                    

    # Calculate a final CaAr and REE corrected 86 beam
    Sr86 = (SrCaArYb86 - ((CaAr82 * .004 / .647) / np.power((85.9160721 / 81.9210049), BetaSr)) - (Yb86_5 * 21.68 /16.103) / np.power((85.968195 / 86.46911), BetaSr))

    # Gather up intermediate channels and add them as time series:
    int_channel_names = ['PFract', 'CaAr82', 'PSrCaArYb86', 'PSrCaArYbLu88',
            'BetaSr', 'BetaRb', 'Rb87', 'Sr87', 'BetaCaAr', 'CaArErYb84',
            'Sr84', 'Sr88', 'Sr86']
    int_channels = [PFract, CaAr82, PSrCaArYb86, PSrCaArYbLu88,
            BetaSr, BetaRb, Rb87, Sr87, BetaCaAr, CaArErYb84,
            Sr84, Sr88, Sr86]
    for name, channel in zip(int_channel_names, int_channels):
        data.createTimeSeries(name, data.Intermediate, indexChannel.time(), channel)

    drs.message("Calculating further interference corrections...")
    drs.progress(60)

    Sr8786_Uncorr = (SrRbYb87 / SrCaArYb86) * np.power((86.9089 / 85.9093), BetaSr) * mask
    Sr8786_Corr = (Sr87 / Sr86) * np.power((86.9089 / 85.9093), BetaSr) * mask
    Rb87Sr86_Corr = (Rb87 / Sr86) * np.power((86.9089 / 85.9093), BetaSr) * mask
    Sr8486_Uncorr = (SrCaArErYb84 / SrCaArYb86) * np.power((83.9134 / 85.9093), BetaSr) * mask
    Sr8486_Corr = (Sr84 / Sr86) * np.power((83.9134 / 85.9093), BetaSr) * mask
    Sr8488_Uncorr = (SrCaArErYb84 / SrCaArYbLu88) * np.power((83.9134 / 87.9056), BetaSr) * mask
    Sr8488_Corr = (Sr84 / Sr88) * np.power((83.9134 / 87.9056), BetaSr) * mask
    Rb87asPPM = (Rb87 / SrRbYb87) * 1000000 * mask
    CaArErYb84asPPM = (CaArErYb84 / SrCaArErYb84) * 100000 * mask
    TotalSrBeam = Sr88 + Sr84 + Sr86 + Sr87 * mask

    # The following equations subtract CaAr using the signal on 83, itself corrected for REE2+ using the Er83_5 signal.
    CaAr83 = CaArEr83 - (Er83_5 * 33.503 /22.869) / np.power((82.96514975 / 83.46619), PFract)
    # Use this preliminary fract in stripping CaAr and REE from Sr 86
    PSrCaArYb86_b = (SrCaArYb86 - ((CaAr83 * .004 / .135) / np.power((85.9160721 / 82.921150), PFract)) - (Yb86_5 * 21.68 / 16.103) / np.power((85.968195 / 86.46911), PFract))   
    # and Sr 88
    PSrCaArYbLu88_b = (SrCaArYbLu88 - ((CaAr83 * .187 / .135) / np.power((87.9149151 / 82.921150), PFract)) - (Yb86_5 * 12.996 / 16.103) / np.power((87.97129 / 86.46911), PFract) - ((Yb86_5 * 12.996 / 16.103) / np.power((87.97129 / 86.46911), PFract) * Dy_Er * 2.59))

    # Use these CaAr and REE stripped Sr 86 and Sr 88 values to calculate a refined fractionation factor
    BetaSr_b = (np.log(Sr88_86_reference/(PSrCaArYbLu88_b / PSrCaArYb86_b))) / (np.log(87.9056/85.9093))
    BetaRb_b = BetaSr_b + RbBias

    # Use the Rb fract to calculate the amount of Rb on mass 87
    Rb87_b = ((RbYbEr85 - (Yb86_5 * 2.982 / 16.103) / np.power((84.967385 / 86.46911), BetaSr_b) - (Er83_5 * 14.910 /22.869) / np.power((84.96774 / 83.46619), BetaSr_b)) * Rb87_85_reference) / np.power((86.90918 / 84.9118), BetaRb_b)
    # Subtract this amount of Rb 87 from the 87 beam
    Sr87_b = (SrRbYb87 - Rb87_b - (Yb86_5 * 32.026 /16.103) / np.power((86.96943 / 86.46911), BetaSr_b))
    # Allows for a modification of the CaAr fractionation factor relative to the Sr fract. We have never used anything but 1 (i.e. no modification)
    BetaCaAr_b = BetaSr_b + CaArBias
    # Calculate the amount of CaAr and REE on mass 84
    CaArErYb84_b = (CaAr83 * ((2.086/0.135) / np.power((83.917989 / 82.921150), BetaCaAr_b))) + (Er83_5 * 26.978 / 22.869) / np.power((83.96619 / 83.46619), BetaSr_b) + (Yb86_5 * 0.123 /16.103) / np.power((83.96694 / 86.46911), BetaSr_b)
    # Subtract this amount from the 84 beam
    Sr84_b = (SrCaArErYb84 - CaArErYb84_b)
    # Calculate a final CaAr and REE corrected 88 beam
    Sr88_b = (SrCaArYbLu88 - ((CaAr83 * .187 / .135) / np.power((87.9149151 / 82.921150), BetaSr_b)) - (Yb86_5 * 12.996 / 16.103) / np.power((87.97129 / 86.46911), BetaSr_b) - ((Yb86_5 * 12.996 / 16.103) / np.power((87.97129 / 86.46911), BetaSr_b) * Dy_Er * 2.59))
    # Calculate a final CaAr and REE corrected 86 beam
    Sr86_b = (SrCaArYb86 - ((CaAr83 * .004 / .135) / np.power((85.9160721 / 82.921150), BetaSr_b)) - (Yb86_5 * 21.68 / 16.103) / np.power((85.968195 / 86.46911), BetaSr_b))

    Sr8786_Uncorr_b = (SrRbYb87 / SrCaArYb86) * np.power((86.9089 / 85.9093), BetaSr_b) * mask
    Sr8786_Corr_b = (Sr87_b / Sr86_b) * np.power((86.9089 / 85.9093), BetaSr_b) * mask
    Rb87Sr86_Corr_b = (Rb87_b / Sr86_b) * np.power((86.9089 / 85.9093), BetaSr_b) * mask
    Sr8486_Uncorr_b = (SrCaArErYb84 / SrCaArYb86) * np.power((83.9134 / 85.9093), BetaSr_b) * mask
    Sr8486_Corr_b = (Sr84_b / Sr86_b) * np.power((83.9134 / 85.9093), BetaSr_b) * mask
    Sr8488_Uncorr_b = (SrCaArErYb84 / SrCaArYbLu88) * np.power((83.9134 / 87.9056), BetaSr_b) * mask
    Sr8488_Corr_b = (Sr84_b / Sr88_b) * np.power((83.9134 / 87.9056), BetaSr_b) * mask
    Rb87asPPM_b = (Rb87_b / SrRbYb87) * 1000000 * mask
    CaArErYb84asPPM_b = (CaArErYb84_b / SrCaArErYb84) * 100000 * mask
    TotalSrBeam_b = Sr88_b + Sr84_b + Sr86_b + Sr87_b * mask

    # Gather up intermediate channels and add them as time series:
    int_channel_names = ['Sr88', 'Sr86', 'Sr84', 'Rb87asPPM']
    int_channel_names += ['Sr88_b', 'Sr86_b', 'Sr84_b', 'Rb87asPPM_b']
    int_channel_names += ['Sr8786_Uncorr', 'Sr8786_Corr', 'Rb87Sr86_Corr', 'Sr8486_Uncorr']
    int_channel_names += ['Sr8486_Corr', 'Sr8488_Uncorr', 'Sr8488_Corr', 'Rb87asPPM']
    int_channel_names += ['CaArErYb84asPPM', 'TotalSrBeam', 'Sr8786_Uncorr_b', 'Sr8786_Corr_b']
    int_channel_names += ['Rb87Sr86_Corr_b', 'Sr8486_Uncorr_b', 'Sr8486_Corr_b', 'Sr8488_Uncorr_b']
    int_channel_names += ['Sr8488_Corr_b', 'Rb87asPPM_b', 'CaArErYb84asPPM_b', 'TotalSrBeam_b']

    int_channels = [Sr88, Sr86, Sr84, Rb87asPPM]
    int_channels += [Sr88_b, Sr86_b, Sr84_b, Rb87asPPM_b]
    int_channels += [Sr8786_Uncorr, Sr8786_Corr, Rb87Sr86_Corr, Sr8486_Uncorr]
    int_channels += [Sr8486_Corr, Sr8488_Uncorr, Sr8488_Corr, Rb87asPPM]
    int_channels += [CaArErYb84asPPM, TotalSrBeam, Sr8786_Uncorr_b, Sr8786_Corr_b]
    int_channels += [Rb87Sr86_Corr_b, Sr8486_Uncorr_b, Sr8486_Corr_b, Sr8488_Uncorr_b]
    int_channels += [Sr8488_Corr_b, Rb87asPPM_b, CaArErYb84asPPM_b, TotalSrBeam_b]

    for name, channel in zip(int_channel_names, int_channels):
        data.createTimeSeries(name, data.Intermediate, indexChannel.time(), channel)

    drs.message("Calculating reference material corrected results...")
    drs.progress(70)

    # Now check if there are selections for the reference standard, and if so, generate standard-normalised ratios
    try:
        StdSpline_Sr87_86 = data.spline(rmName, "Sr8786_Corr").data()
        StdValue_Sr87_86 = data.referenceMaterialData(rmName)["87Sr/86Sr"].value()
    except:
        IoLog.error("The Combined Sr DRS requires Ref Material selections to proceed.")
        drs.message("DRS did not finish. Please check Messages")
        drs.progress(100)
        drs.finished()
        return

    StdCorr_Sr8786 = Sr8786_Corr * StdValue_Sr87_86 / StdSpline_Sr87_86

    # Now generate age-corrected values (using the observed Rb/Sr ratio)
    # NOTE: The line below used Sr8786_Corr. Changed to use
    # Std Corrected data (BP, 2021-02-01)
    Sr8786_AgeCorr = StdCorr_Sr8786 - (Rb87Sr86_Corr * (0.000013972 ** Age) - 1)

    try:
        StdSpline_Sr8786_b = data.spline(rmName, "Sr8786_Corr_b").data()
    except:
        IoLog.error("The Combined Sr DRS requires Ref Material selections to proceed.")
        drs.message("DRS did not finish. Please check Messages")
        drs.progress(100)
        drs.finished()
        return

    StdCorr_Sr8786_b = Sr8786_Corr_b * StdValue_Sr87_86 / StdSpline_Sr8786_b
    # Now generate age-corrected values (using the observed Rb/Sr ratio)
    Sr8786_AgeCorr_b = Sr8786_Corr_b - (Rb87Sr86_Corr_b * (0.000013972 ** Age) - 1)

    output_channels_names = ['StdCorr_Sr8786', 'Sr8786_AgeCorr', 'StdCorr_Sr8786_b', 'Sr8786_AgeCorr_b']
    output_channels = [StdCorr_Sr8786, Sr8786_AgeCorr, StdCorr_Sr8786_b, Sr8786_AgeCorr_b]
    for name, channel in zip(output_channels_names, output_channels):
        data.createTimeSeries(name, data.Output, indexChannel.time(), channel)

    if propErrors:
        drs.message("Propagating errors...")
        drs.progress(90)

        groups = [s for s in data.selectionGroupList() if s.type != data.Baseline]
        data.propagateErrors(groups, [data.timeSeries("StdCorr_Sr87_86_b")], data.timeSeries("Sr8786_Corr_b"), rmName)


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
    if 'TotalBeam' in timeSeriesNames:
        defaultChannelName = 'TotalBeam'
    elif timeSeriesNames:
        defaultChannelName = timeSeriesNames[0]

    rmNames = data.selectionGroupNames(data.ReferenceMaterial)

    drs.setSetting("IndexChannel", defaultChannelName)
    drs.setSetting("ReferenceMaterial", "CO3_shell")
    drs.setSetting("Mask", True)
    drs.setSetting("MaskChannel", defaultChannelName)
    drs.setSetting("MaskCutoff", 0.05)
    drs.setSetting("MaskTrim", 0.0)
    drs.setSetting("Age", 0.)
    drs.setSetting("Dy_Er", 1.8)
    drs.setSetting("RbBias", 0.)
    drs.setSetting("CaArBias", 0.)
    drs.setSetting("Sr88_86_reference", 8.37520938)  #Konter & Storm (2014)
    drs.setSetting("Rb87_85_reference", 0.385710)    #Konter & Storm (2014)

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

    ageLineEdit = QtGui.QLineEdit(widget)
    ageLineEdit.setText(settings["Age"])
    ageLineEdit.textChanged.connect(lambda t: drs.setSetting("Age", float(t)))
    formLayout.addRow("Age", ageLineEdit)

    rbBiasLineEdit = QtGui.QLineEdit(widget)
    rbBiasLineEdit.setText(settings["RbBias"])
    rbBiasLineEdit.textChanged.connect(lambda t: drs.setSetting("RbBias", float(t)))
    formLayout.addRow("Rb Bias Adjustment", rbBiasLineEdit)

    caArBiasLineEdit = QtGui.QLineEdit(widget)
    caArBiasLineEdit.setText(settings["CaArBias"])
    caArBiasLineEdit.textChanged.connect(lambda t: drs.setSetting("CaArBias", float(t)))
    formLayout.addRow("CaAr Bias Adjustment", caArBiasLineEdit)

    dyErLineEdit = QtGui.QLineEdit(widget)
    dyErLineEdit.setText(settings["Dy_Er"])
    dyErLineEdit.textChanged.connect(lambda t: drs.setSetting("Dy_Er", float(t)))
    formLayout.addRow("Dy/Er ratio", dyErLineEdit)

    sr88_86_refLineEdit = QtGui.QLineEdit(widget)
    sr88_86_refLineEdit.setText(settings["Sr88_86_reference"])
    sr88_86_refLineEdit.textChanged.connect(lambda t: drs.setSetting("Sr88_86_reference", float(t)))
    formLayout.addRow("Reference Sr88/Sr86 value", sr88_86_refLineEdit)

    rb87_85_refLineEdit = QtGui.QLineEdit(widget)
    rb87_85_refLineEdit.setText(settings["Rb87_85_reference"])
    rb87_85_refLineEdit.textChanged.connect(lambda t: drs.setSetting("Rb87_85_reference", float(t)))
    formLayout.addRow("Reference Rb87/Rb85 value", rb87_85_refLineEdit)

    propCheckBox = QtGui.QCheckBox(widget)
    propCheckBox.setChecked(settings["PropagateError"])
    propCheckBox.toggled.connect(lambda t: drs.setSetting("PropagateError", bool(t)))
    formLayout.addRow("Propagate Errors?", propCheckBox)

    drs.setSettingsWidget(widget)
