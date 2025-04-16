# A python-based data reduction scheme for iolite 4 starts with some metadata
#/ Type: DRS
#/ Name: S_Isotopes
#/ Authors: B Paul & Joe Petrus
#/ Description: A sulfur isotope example DRS
#/ References:
#/ Version:
#/ Contact:

from iolite import QtGui
from iolite.types import Result
from time import sleep
import numpy as np

'''
This function calculates the final, total uncertainty based on the Gilbert el al. (2014) paper:
https://doi.org/10.1039/C4JA00011K

The forumla (from Eq  2 in the paper) is:

total_uncert^2 = sigma^2 / n + u_RM^2 + u_D^2

where total_uncert is the total uncertainty
sigma is the RSD of the sample (in ‰)
n is the number of sweeps in the sample
u_RM is the uncertainty of δ34S (in ‰) for the reference material
and u_D is the uncertainty of the mass bias correction, defined in the paper
as the relative standard error of the mean of the drift corrected RM analyses.
'''

def totalSUncertainty(sel):
    result = Result()

    try:
        corr_3432 = data.timeSeries("Corrected 34S/32S")
        vals = corr_3432.dataForSelection(sel)
        sigma = np.std(vals)
        n = len(vals)
        u_RM = float(sel.property("u_RM"))
        u_D = float(sel.property("u_D"))
        total_uncert = np.sqrt(sigma**2. / n + u_RM**2. + u_D**2.)
        result.setValue(total_uncert)
        return result

    except:
        print(f'Could not calculate total uncertainty for selection {sel.name}')
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
    """

    drs.message("Starting baseline subtract DRS...")
    drs.progress(0)

    # Get settings
    settings = drs.settings()
    print(settings)

    indexChannel = data.timeSeries(settings["IndexChannel"])
    maskChannel = data.timeSeries(settings["MaskChannel"])
    cutoff = settings["MaskCutoff"]
    trim = settings["MaskTrim"]

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

    bs_list = data.selectionGroupList(data.Baseline)
    if(len(bs_list) == 0):
        drs.message("No baseline selection group found. Please create one before continuing.")
        drs.progress(100)
        drs.finished()
        return
    if(len(bs_list) > 1):
        drs.message("Multiple baseline selection groups found. Please ensure there is only one baseline group.")
        drs.progress(100)
        drs.finished()
        return

    bs_grp = bs_list[0]

    allInputChannels = data.timeSeriesList(data.Input)
    for counter, channel in enumerate(allInputChannels):
        drs.baselineSubtract(bs_grp, allInputChannels, mask, 25, 40)

    drs.message("Calculating ratios...")
    drs.progress(45)

    s32_CPS = data.timeSeriesList(
        data.Intermediate, {'Pre-shift mass': '32'})[0].data()
    s34_CPS = data.timeSeriesList(
        data.Intermediate, {'Pre-shift mass': '34'})[0].data()
    s33_measured = False
    try:
        s33_CPS = data.timeSeriesList(
            data.Intermediate, {'Pre-shift mass': '33'})[0].data()
        s33_measured = True
    except:
        pass

    raw_34_32 = s34_CPS / s32_CPS * mask
    if s33_measured:
        raw_33_32 = s33_CPS / s32_CPS * mask

    data.createTimeSeries("Raw 34S/32S", data.Intermediate, indexChannel.time(), raw_34_32)
    if s33_measured:
        data.createTimeSeries("Raw 33S/32S", data.Intermediate, indexChannel.time(), raw_33_32)

    drs.message("Calculating mass bias...")
    drs.progress(50)

    # Check that we have a reference material group with selections
    rmName = settings["ReferenceMaterial"]
    try:
        data.selectionGroup(rmName)
    except:
        IoLog.error(f"Could not find the {rmName} selection group. Please create one before continuing.")
        drs.message("DRS did not finish. Please check Messages")
        drs.progress(100)
        return

    # Now calculate the mass bias for the 34S/32S ratio
    try:
        StdSpline_S34_S32 = data.spline(rmName, "Raw 34S/32S").data()
    except:
        IoLog.error("Could not find the spline for the 34S/32S for the reference material selection group.")
        drs.message("DRS did not finish. Please check Messages")
        drs.progress(100)
        drs.finished()
        return

    try:
        StdValue_S34_S32 = data.referenceMaterialData(rmName)["34S/32S"].value()
    except:
        IoLog.error("Could not find the \"34S/32S\" value for the reference material selection group. Please make sure it is in the RM file (exactly as written here).")
        drs.message("DRS did not finish. Please check Messages")
        drs.progress(100)
        drs.finished()
        return

    print(f'StdValue_S34_S32: {StdValue_S34_S32}')

    mb_34_32 = StdValue_S34_S32 / StdSpline_S34_S32

    # Now repeat for the 33S/32S ratio
    if s33_measured:
        try:
            StdSpline_S33_S32 = data.spline(rmName, "Raw 33S/32S").data()
        except:
            IoLog.error("Could not find the spline for the 33S/32S for the reference material selection group.")
            drs.message("DRS did not finish. Please check Messages")
            drs.progress(100)
            drs.finished()
            return

        try:
            StdValue_S33_S32 = data.referenceMaterialData(rmName)["33S/32S"].value()
        except:
            IoLog.error("Could not find the \"33S/32S\" value for the reference material selection group. Please make sure it is in the RM file (exactly as written here).")
            drs.message("DRS did not finish. Please check Messages")
            drs.progress(100)
            drs.finished()
            return

        print(f'StdValue_S33_S32: {StdValue_S33_S32}')

        mb_33_32 = StdValue_S33_S32 / StdSpline_S33_S32

    # Now make them available for viewing
    data.createTimeSeries("34S/32S Mass Bias", data.Intermediate, indexChannel.time(), mb_34_32)
    if s33_measured:
        data.createTimeSeries("33S/32S Mass Bias", data.Intermediate, indexChannel.time(), mb_33_32)

    # Now calculate the corrected ratios
    drs.message("Calculating corrected ratios...")
    drs.progress(60)

    corrected_34_32 = raw_34_32 * mb_34_32
    if s33_measured:
        corrected_33_32 = raw_33_32 * mb_33_32

    data.createTimeSeries(
        "Corrected 34S/32S", data.Intermediate,
        indexChannel.time(), corrected_34_32)

    if s33_measured:
        data.createTimeSeries(
            "Corrected 33S/32S", data.Intermediate,
            indexChannel.time(), corrected_33_32)

    # Now calculate the del channels
    drs.message("Calculating del values...")
    drs.progress(70)

    del_34_32 = (corrected_34_32 / 0.0441626 - 1) * 1000.
    if s33_measured:
        del_33_32 = (corrected_33_32 / 0.00787724 - 1) * 1000.

    data.createTimeSeries("del 34S/32S", data.Output,
                          indexChannel.time(), del_34_32)

    if s33_measured:
        data.createTimeSeries("del 33S/32S", data.Output,
                              indexChannel.time(), del_33_32)

    # Now calculate the delta channel
    if s33_measured:
        drs.message("Calculating delta values...")
        drs.progress(80)

        delta_33S = del_33_32 - 1000. * (np.power((del_34_32/1000. + 1.), 0.515) - 1.)
        data.createTimeSeries("delta 33S", data.Output, indexChannel.time(), delta_33S)

    # Now calculate the total uncertainty as an associated result (see totalSUncertainty function at top of this script)
    # Will store some of the parameters needed as selection properties
    rm_res = data.referenceMaterialData(rmName)["34S/32S"]
    # This is the relative standard error of the reference material
    u_RM = rm_res.uncertainty() / rm_res.value() * (rm_res.value())
    rm_del_vals = data.frame([data.selectionGroup(rmName)], ['Name'], [data.timeSeries('del 34S/32S')], ['mean', '2SE'])['del 34S/32S mean'].values
    # This is the relative standard error of the mean
    u_D = np.std(rm_del_vals) / (np.sqrt(len(rm_del_vals)) - 1.)

    print(f'u_RM: {u_RM}')
    print(f'u_D: {u_D}')

    for grp in data.selectionGroupList():
        if grp.type == data.Baseline:
            continue
        for sel in grp.selections():
            sel.setProperty("u_RM", u_RM)
            sel.setProperty("u_D", u_D)

    data.registerAssociatedResult("S_Iso_TotalUncert", totalSUncertainty)

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

    drs.setSetting("IndexChannel", defaultChannelName)
    drs.setSetting("MaskChannel", defaultChannelName)
    drs.setSetting("MaskCutoff", 500000.0)
    drs.setSetting("MaskTrim", 0.0)
    drs.setSetting("ReferenceMaterial", "My_S_RM")

    settings = drs.settings()

    indexComboBox = QtGui.QComboBox(widget)
    indexComboBox.addItems(data.timeSeriesNames(data.Input))
    indexComboBox.setCurrentText(settings["IndexChannel"])
    indexComboBox.currentTextChanged.connect(lambda t: drs.setSetting("IndexChannel", t))

    def updateIndexCombo():
        timeSeriesNames = data.timeSeriesNames(data.Input)
        indexComboBox.clear()
        indexComboBox.addItems(timeSeriesNames)

    data.dataChanged.connect(updateIndexCombo)

    rmComboBox = QtGui.QComboBox(widget)
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
        if len(rmNames) == 1:
            print("Only one RM found, setting to: ", rmNames[0])
            rmComboBox.setCurrentText(rmNames[0])
        else:
             print("Multiple RMs found, setting to: ", len(rmNames))

        drs.setSetting("ReferenceMaterial", rmComboBox.currentText)

    data.selectionGroupsChanged.connect(updateRMCombo)

    updateRMCombo()

    verticalSpacer1 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
    formLayout.addItem(verticalSpacer1)

    maskComboBox = QtGui.QComboBox(widget)
    maskComboBox.addItems(data.timeSeriesNames(data.Input))
    maskComboBox.setCurrentText(settings["MaskChannel"])
    maskComboBox.currentTextChanged.connect(lambda t: drs.setSetting("MaskChannel", t))

    maskLineEdit = QtGui.QLineEdit(widget)
    maskLineEdit.setText(settings["MaskCutoff"])
    maskLineEdit.textChanged.connect(lambda t: drs.setSetting("MaskCutoff", float(t)))

    maskTrimLineEdit = QtGui.QLineEdit(widget)
    maskTrimLineEdit.setText(settings["MaskTrim"])
    maskTrimLineEdit.textChanged.connect(lambda t: drs.setSetting("MaskTrim", float(t)))

    formLayout.addRow("Index channel", indexComboBox)
    formLayout.addRow("Mask channel", maskComboBox)
    formLayout.addRow("Mask cutoff", maskLineEdit)
    formLayout.addRow("Mask trim", maskTrimLineEdit)

    drs.setSettingsWidget(widget)
