# A python-based data reduction scheme for iolite 4 starts with some metadata
#/ Type: DRS
#/ Name: RbSr_Neoma
#/ Authors: Grant Craig and Bence Paul
#/ Description: RbSr DRS for Neoma MSMS, includes approx age calc
#/ References: None
#/ Version: 1.0
#/ Contact: grant.craig@thermofisher.com

import math
import numpy as np
from scipy import constants

from iolite import QtGui
from iolite.types import Result

'''
This function is for calculating associated results (i.e. results
for a selection that are not directly tied to the values of a single channel)
This function calculates an age from the individual data points within a selection
'''
def calcSelectionAge(sel):
    result = Result()
    try:
        Rb87_CPS = data.timeSeriesList(
            data.Intermediate, {'Post-shift mass': '87'})[0]
        Sr86_CPS = data.timeSeriesList(
            data.Intermediate, {'Post-shift mass': '105'})[0]
        Sr87_CPS = data.timeSeriesList(
            data.Intermediate, {'Post-shift mass': '106'})[0]

        #Sr87_Sr86 = data.timeSeries('StdCorr_Sr87s_Sr86s')
        Sr87_Sr86 = data.timeSeries('StdCorr_Sr87_Sr86_MBC')
        #Rb87_Sr86 = data.timeSeries('StdCorr_Rb87_Sr86s')
        Rb87_Sr86 = data.timeSeries('StdCorr_Rb87_Sr86_MBC')
    except RuntimeError:
        print("Could not get the channels needed for plotting iso in Rb Sr")

    # Resistor values for Rb87, Sr87 and Sr86
    R_Rb87 = 1.0e+11
    R_Sr87 = 1.0e+13
    R_Sr86 = 1.0e+13
    # Boltzman constant
    K = constants.k
    # Temperature (in K)
    TK = 308.15
    # integration time (in s)
    it = 1.
    # Volts to CPS factor
    V2CPS = 6.2414e+07

    # If getting from settings... TODO
    #settings = drs.settings()
    #R_Rb87 = settings['R_Rb87']
    #R_Sr87 = settings['R_Sr87']
    #R_Sr86 = settings['R_Sr86']
    #TK = settings['T_K']
    #it = settings['IntegTime']
    #V2CPS = settings['Volts2CPS']
    #l = settings['RbSrLambda']

    # Calculate point uncertainties for Rb87, Sr87, and Sr86
    # Amplifier noise on Rb87 (in CPS)
    s_amp_Rb87 = math.sqrt(4. * K * TK * R_Rb87 / it) * V2CPS
    # Total noise on Rb87 (as RSE)
    s_Rb87 = np.sqrt((s_amp_Rb87 / Rb87_CPS.dataForSelection(sel)) ** 2 +
                    (1 / np.sqrt(Rb87_CPS.dataForSelection(sel) + it))**2)

    # Sr87
    # Amplifier noise on Sr87 (in CPS) NOTE: Different form to basic Johnson
    # eqn. See Grant Craig's spreadsheet
    s_amp_Sr87 = math.sqrt(4. * K * TK * R_Sr87 / it) / R_Sr87 * R_Rb87 * V2CPS
    # Total noise on Sr87 (as RSE)
    s_Sr87 = np.sqrt((s_amp_Sr87 / Sr87_CPS.dataForSelection(sel)) ** 2 +
                    (1 / np.sqrt(Sr87_CPS.dataForSelection(sel) + it))**2)

    # Sr86
    # Amplifier noise on Sr86 (in CPS)
    s_amp_Sr86 = math.sqrt(4. * K * TK * R_Sr86 / it) / R_Sr87 * R_Rb87 * V2CPS
    # Total noise on Sr86 (as RSE)
    s_Sr86 = np.sqrt((s_amp_Sr86 / Sr86_CPS.dataForSelection(sel)) ** 2 +
                    (1 / np.sqrt(Sr86_CPS.dataForSelection(sel) + it))**2)

    s_Rb87_86 = np.sqrt(s_Rb87 ** 2. + s_Sr86 ** 2.) * Rb87_Sr86.dataForSelection(sel) * 2.
    s_Sr87_86 = np.sqrt(s_Sr87 ** 2. + s_Sr86 ** 2.) * Sr87_Sr86.dataForSelection(sel) * 2.

    rho = np.corrcoef(Rb87_Sr86.dataForSelection(sel),
                        Sr87_Sr86.dataForSelection(sel))[0,1]

    fit = iolite_helpers.fitLine(Rb87_Sr86.dataForSelection(sel), s_Rb87_86,
                                    Sr87_Sr86.dataForSelection(sel), s_Sr87_86, rho)

    # TODO: get value from settings:
    #
    #lambda
    l = 1.3981E-11

    t = np.log(fit['m'] + 1) /  l / 1E6
    s_t = 2. * np.log(fit['sigma_m'] + 1) /  l / 1E6

    print(f'{sel.name} Age: {t:.2f} ± {s_t:.2f}')

    result.setValue(t)
    result.setUncertainty(s_t)
    return result


def runDRS():
    """
    This method will be called by iolite when the user clicks
    Crunch Data in the DRS window or as part of a processing
    template. It should transform 'input' data into 'output'
    data using the provided settings.

    DRS progress can be updated via the 'message' and 'progress'
    signals. These will be displayed in the iolite interface.

    When finished, the 'finished' signal should be emitted.

    As an example, we will do baseline subtraction of all
    input channels using a DRS helper function.
    """

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

    drs.baselineSubtract(blGrp, allInputChannels, mask, 25, 75)
    for ch in data.timeSeriesList(data.Input):
        cps_ch = data.timeSeries(ch.name + '_CPS')
        cps_ch.setProperty(
            'Pre-shift mass', ch.property('Pre-shift mass'))
        cps_ch.setProperty('Post-shift mass',
                            ch.property('Post-shift mass'))

    drs.message("Calculating raw ratios...")
    drs.progress(50)


    # Declare the channels used in the calculations:

    '''
    Probably better to use this approx, with post-shift mass, but
    leaving it as it is to avoid confusion with Grant's calculations
    '''
    # Rb87_CPS = data.timeSeriesList(
        # data.Intermediate, {'Post-shift mass': '87'})[0].data()
    # Sr86_105_CPS = data.timeSeriesList(
        # data.Intermediate, {'Post-shift mass': '105'})[0].data()
    # Sr87_106_CPS = data.timeSeriesList(
        # data.Intermediate, {'Post-shift mass': '106'})[0].data()
    # Sr88_107_CPS = data.timeSeriesList(
        # data.Intermediate, {'Post-shift mass': '107'})[0].data()

    #Rb85_CPS = data.timeSeriesList(
        #data.Intermediate, {'Element': 'Rb', 'Mass' : '85'})[0].data()
    Rb87_CPS = data.timeSeriesList(
        data.Intermediate, {'Element': 'Rb', 'Mass' : '87'})[0].data()
    #Sr86_CPS = data.timeSeriesList(
        #data.Intermediate, {'Element': 'Sr', 'Mass' : '86'})[0].data()
    Sr88_CPS = data.timeSeriesList(
        data.Intermediate, {'Element': 'Sr', 'Mass' : '88'})[0].data()
    Sr84F_CPS = data.timeSeriesList(
        data.Intermediate, {'Element': 'F', 'Mass' : '84'})[0].data()
    Sr86F_CPS = data.timeSeriesList(
        data.Intermediate, {'Element': 'F', 'Mass' : '86'})[0].data()
    Sr87F_CPS = data.timeSeriesList(
        data.Intermediate, {'Element': 'F', 'Mass' : '87'})[0].data()
    Sr88F_CPS = data.timeSeriesList(
        data.Intermediate, {'Element': 'F', 'Mass' : '88'})[0].data()

    #Rb87_CPS = Rb85_CPS * 0.38562
    Rb87c_CPS = Rb87_CPS - Sr88_CPS * Sr87F_CPS/Sr88F_CPS #use either this line, or line above, depending on cup configuration.
    #Rb85_Sr86s_Raw = Rb85_CPS/Sr86F_CPS
    Sr87_Sr86_Raw = Sr87F_CPS/Sr86F_CPS
    Rb87_Sr86_Raw = Rb87c_CPS/Sr86F_CPS
    Sr87_Rb87_Raw = Sr87F_CPS/Rb87c_CPS
    Sr88_Sr86_Raw = Sr88F_CPS/Sr86F_CPS
    Sr84_Sr86_Raw = Sr84F_CPS/Sr86F_CPS
    Beta = np.log(8.37520938/(Sr88F_CPS/Sr86F_CPS))/np.log(87.9056/85.9093)
    Sr87m_Sr86m_Raw = Sr87F_CPS/Sr86F_CPS * pow(86.9089/85.9093, Beta)
    Rb87_Sr86m_Raw = Rb87c_CPS/Sr86F_CPS * pow(86.9089/85.9093, Beta)
    Sr84m_Sr86m_Raw = Sr84F_CPS/Sr86F_CPS * pow(83.9134/85.9093, Beta)

    # Gather up intermediate channels and add them as time series:
    int_channel_names = ['Rb87c_CPS', 'Sr87_Sr86_Raw', 'Rb87_Sr86_Raw', 'Sr87_Rb87_Raw', 'Sr88_Sr86_Raw', 'Sr84_Sr86_Raw', 'Beta', 'Sr87m_Sr86m_Raw', 'Rb87_Sr86m_Raw', 'Sr84m_Sr86m_Raw']
    int_channels = [Rb87c_CPS, Sr87_Sr86_Raw, Rb87_Sr86_Raw, Sr87_Rb87_Raw, Sr88_Sr86_Raw, Sr84_Sr86_Raw, Beta, Sr87m_Sr86m_Raw, Rb87_Sr86m_Raw, Sr84m_Sr86m_Raw]
    for name, channel in zip(int_channel_names, int_channels):
        data.createTimeSeries(name, data.Intermediate,
                              indexChannel.time(), channel, {'DRSSettings': settings})

    drs.message("Correcting ratios...")
    drs.progress(80)

    StdSpline_Rb87_Sr86 = data.spline(rmName, "Rb87_Sr86_Raw").data()
    try:
        StdValue_Rb87_Sr86 = data.referenceMaterialData(rmName)["87Rb/86Sr"].value()
    except KeyError:
        IoLog.error("There was no 87Rb/86Sr value in the " + rmName +
                      " datafile. Rb-Sr DRS cannot proceed.")
        drs.message("Error. See Messages")
        drs.progress(100)
        drs.finished()
        return

    print("StdSpline_Rb87_Sr86 mean = %f" % StdSpline_Rb87_Sr86.mean())
    print("StdValue_Rb87_Sr86 = %f" % StdValue_Rb87_Sr86)

    StdCorr_Rb87_Sr86 = (Rb87_Sr86_Raw) * StdValue_Rb87_Sr86 / StdSpline_Rb87_Sr86
    data.createTimeSeries('StdCorr_Rb87_Sr86', data.Output,
                          indexChannel.time(), StdCorr_Rb87_Sr86)

    StdSpline_Sr87_Sr86 = data.spline(rmName, "Sr87_Sr86_Raw").data()
    try:
        StdValue_Sr87_Sr86 = data.referenceMaterialData(rmName)["87Sr/86Sr"].value()
    except KeyError:
        IoLog.error("There was no 87Sr/86Sr value in the " + rmName +
                      " datafile. Rb-Sr DRS cannot proceed.")
        drs.message("Error. See Messages")
        drs.progress(100)
        drs.finished()
        return

    print("StdSpline_Sr87_Sr86 mean = %f" % StdSpline_Sr87_Sr86.mean())
    print("StdValue_Sr87_Sr86 = %f" % StdValue_Sr87_Sr86)

    StdCorr_Sr87_Sr86 = (Sr87_Sr86_Raw) * StdValue_Sr87_Sr86 / StdSpline_Sr87_Sr86
    data.createTimeSeries('StdCorr_Sr87_Sr86', data.Output,
                          indexChannel.time(), StdCorr_Sr87_Sr86)

    StdSpline_Rb87_Sr86m = data.spline(rmName, "Rb87_Sr86m_Raw").data()
    try:
        StdValue_Rb87_Sr86m = data.referenceMaterialData(rmName)["87Rb/86Sr"].value()
    except KeyError:
        IoLog.error("There was no 87Rb/86Sr value in the " + rmName +
                      " datafile. Rb-Sr DRS cannot proceed.")
        drs.message("Error. See Messages")
        drs.progress(100)
        drs.finished()
        return

    print("StdSpline_Rb87_Sr86m mean = %f" % StdSpline_Rb87_Sr86m.mean())
    print("StdValue_Rb87_Sr86m = %f" % StdValue_Rb87_Sr86m)

    StdCorr_Rb87_Sr86_MBC = (Rb87_Sr86m_Raw) * StdValue_Rb87_Sr86m / StdSpline_Rb87_Sr86m
    data.createTimeSeries('StdCorr_Rb87_Sr86_MBC', data.Output,
                          indexChannel.time(), StdCorr_Rb87_Sr86_MBC)

    StdSpline_Sr87m_Sr86m = data.spline(rmName, "Sr87m_Sr86m_Raw").data()
    try:
        StdValue_Sr87m_Sr86m = data.referenceMaterialData(rmName)["87Sr/86Sr"].value()
    except KeyError:
        IoLog.error("There was no 87Sr/86Sr value in the " + rmName +
                      " datafile. Rb-Sr DRS cannot proceed.")
        drs.message("Error. See Messages")
        drs.progress(100)
        drs.finished()
        return

    print("StdSpline_Sr87m_Sr86m mean = %f" % StdSpline_Sr87m_Sr86m.mean())
    print("StdValue_Sr87m_Sr86m = %f" % StdValue_Sr87m_Sr86m)

    StdCorr_Sr87_Sr86_MBC = (Sr87m_Sr86m_Raw) * StdValue_Sr87m_Sr86m / StdSpline_Sr87m_Sr86m
    data.createTimeSeries('StdCorr_Sr87_Sr86_MBC', data.Output,
                          indexChannel.time(), StdCorr_Sr87_Sr86_MBC)

    print("Registering associated results...")
    data.registerAssociatedResult("Selection_Age", calcSelectionAge)

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
    drs.setSetting("R_Rb87", 1.0e+11)
    drs.setSetting("R_Sr87", 1.0e+13)
    drs.setSetting("R_Sr86", 1.0e+13)
    drs.setSetting("T_K", 308.15)
    drs.setSetting("IntegTime", 1.)
    drs.setSetting("Volts2CPS", 6.2414E+07)
    drs.setSetting("RbSrLambda", 1.3981E-11)

    settings = drs.settings()

    indexComboBox = QtGui.QComboBox(widget)
    indexComboBox.addItems(timeSeriesNames)
    indexComboBox.setCurrentText(settings["IndexChannel"])
    indexComboBox.currentTextChanged.connect(
        lambda t: drs.setSetting("IndexChannel", t))
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

    verticalSpacer = QtGui.QSpacerItem(
        20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
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
        20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
    formLayout.addItem(verticalSpacer2)

    # Resistance values for individual datapoint uncertainties
    R_Rb87 = QtGui.QLineEdit(widget)
    R_Rb87.setText(settings["R_Rb87"])
    R_Rb87.textChanged.connect(
        lambda t: drs.setSetting("R_Rb87", float(t)))
    formLayout.addRow("Rb87 Resistance (Ω)", R_Rb87)

    R_Sr87 = QtGui.QLineEdit(widget)
    R_Sr87.setText(settings["R_Sr87"])
    R_Sr87.textChanged.connect(
        lambda t: drs.setSetting("R_Sr87", float(t)))
    formLayout.addRow("Sr87 Resistance (Ω)", R_Sr87)

    R_Sr86 = QtGui.QLineEdit(widget)
    R_Sr86.setText(settings["R_Sr86"])
    R_Sr86.textChanged.connect(
        lambda t: drs.setSetting("R_Sr86", float(t)))
    formLayout.addRow("Sr86 Resistance (Ω)", R_Sr86)

    T_K_LineEdit = QtGui.QLineEdit(widget)
    T_K_LineEdit.setText(settings["T_K"])
    T_K_LineEdit.textChanged.connect(
        lambda t: drs.setSetting("T_K", float(t)))
    formLayout.addRow("Temperature (K)", T_K_LineEdit)

    integTimeLineEdit = QtGui.QLineEdit(widget)
    integTimeLineEdit.setText(settings["IntegTime"])
    integTimeLineEdit.textChanged.connect(
        lambda t: drs.setSetting("IntegTime", float(t)))
    formLayout.addRow("Measurement integration time (s)", integTimeLineEdit)

    v2CPSLineEdit = QtGui.QLineEdit(widget)
    v2CPSLineEdit.setText(settings["Volts2CPS"])
    v2CPSLineEdit.textChanged.connect(
        lambda t: drs.setSetting("Volts2CPS", float(t)))
    formLayout.addRow("Volts to CPS factor", v2CPSLineEdit)

    lambdaLineEdit = QtGui.QLineEdit(widget)
    lambdaLineEdit.setText(settings["RbSrLambda"])
    lambdaLineEdit.textChanged.connect(
        lambda t: drs.setSetting("RbSrLambda", float(t)))
    formLayout.addRow("Rb-Sr λ", lambdaLineEdit)

    drs.setSettingsWidget(widget)
