#/ Type: DRS
#/ Name: Sr Isotopes (Combined)
#/ Authors: Bence Paul, Joe Petrus and author(s) of Sr_isotopes_Total_NIGL.ipf
#/ Description: A Sr isotopes DRS that corrects for REE and CaAr interferences
#/ References: None
#/ Version: 2.0
#/ Contact: support@iolite-software.com

from iolite import QtGui
from iolite import QtCore
from iolite.Qt import Qt, QColor
from iolite.ui import CommonUIPyInterface as CUI
from iolite.ui import IolitePlotPyInterface as Plot
from iolite.ui import IolitePlotSettingsDialog as PlotSettings
from iolite.QtGui import QAction
from iolite.types import Result

from scipy.optimize import curve_fit
import numpy as np
import itertools
from functools import partial


'''
A plot to handle plotting of CaPO correction
If this is the first time the DRS script has been run (e.g.
when you first click on the DRS in iolite) the plot will be
initiated.
It also adds a Settings menu item to the context menu
'''
def showSettings():
    d = PlotSettings(PLOT)
    d.exec_()

try:
    PLOT
except:
    PLOT = Plot()
    PLOT.setAttribute(Qt.WA_DeleteOnClose)
    PLOT.setFixedSize(500,400)

    def showSettings():
        d = PlotSettings(PLOT)
        d.exec_()

    settingsAction = QAction(PLOT.contextMenu())
    settingsAction.setText('Settings')
    settingsAction.triggered.connect(showSettings)
    PLOT.contextMenu().addAction(settingsAction)


'''
A menu to handle which Reference Materials (RMs) to use in the CaPO correction
'''
class ReferenceMaterialsMenu(QtGui.QMenu):

    rmsChanged = QtCore.Signal(list)

    def __init__(self, parent):
        super().__init__(parent)
        self.rmsForCaPOCorrection = []

        for rm in data.referenceMaterialNames():
            a = QtGui.QWidgetAction(self)
            cb = QtGui.QCheckBox(rm, self)
            cb.setStyleSheet('QCheckBox { padding-left: 5px; margin: 3px; }')
            a.setDefaultWidget(cb)
            self.addAction(a)
            cb.clicked.connect(partial(self.updateChannels, rm))

        self.aboutToShow.connect(self.updateMenu)

    def updateMenu(self):
         for a in self.actions():
            cb = a.defaultWidget()
            cb.setChecked(False)
            try:
                if cb.text in self.rmsForCaPOCorrection:
                    cb.setChecked(True)
            except Exception as e:
                print(e)

    def updateChannels(self, rmName, b):
        if b:
            self.rmsForCaPOCorrection = list(set(self.rmsForCaPOCorrection + [rmName]))
        else:
            self.rmsForCaPOCorrection = list(filter(lambda rm: rm != rmName, self.rmsForCaPOCorrection))

        self.rmsChanged.emit(self.rmsForCaPOCorrection)

'''
Defining the colors to use in the CaPO plot here. By default only has 5 colors:
'''
PLOT_COLORS = [
    QColor(239, 71, 111),   # Carmine Red
    QColor(255, 209, 102),  # Orange Yellow
    QColor(6, 214, 160),    # Green
    QColor(17, 138, 178),   # Blue
    QColor(7, 59, 76),      # Midnight Green
]


'''
This function is for calculating associated results (i.e. results
for a selection that are not directly tied to the values of a single channel)
This particular example calculates the correlation between the corrected
87Sr/86Sr ratios and the Rb87 interference for each selection.
'''
def getSelRatioIntensityCorr(sel):
    result = Result()

    try:
        Sr8786_Corr = data.timeSeries("Sr8786_Corr")
        Rb87AsPPM = data.timeSeries("Rb87asPPM")
    except RuntimeError:
        return result

    array_1 = Sr8786_Corr.dataForSelection(sel)
    array_2 = Rb87AsPPM.dataForSelection(sel)

    result.setValue(np.corrcoef(array_1, array_2)[0,1])
    return result



def runDRS():

    drs.message("Starting Sr isotopes DRS...")
    drs.progress(0)

    # Get settings
    settings = drs.settings()
    print(settings)

    if settings["IndexChannel"] not in data.timeSeriesNames(data.Input):
        IoLog.errorWithOrigin("Index channel not found. DRS cannot continue", "Sr Isotopes CaAr REE DRS")
        drs.message("DRS did not finish. Please check Messages")
        drs.progress(100)
        drs.finished()
        return

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

    drs.baselineSubtract(blGrp, allInputChannels, mask, 15, 35)

    drs.message("Checking for half masses...")
    drs.progress(35)

    try:
        Yb86_5 = data.timeSeriesByMass(data.Intermediate, 86.5, 0.1).data()
    except IndexError:
        IoLog.error("Could not find half masses. Combined Sr DRS could not proceed...")
        drs.message("DRS did not finish. Please check Messages")
        drs.progress(100)
        drs.finished()
        return

    #Y89 = data.timeSeriesList(data.Intermediate, {'Mass': '89'})[0].data()
    SrCaArYbLu88 = data.timeSeriesByMass(data.Intermediate, 88, 0.1).data()
    #Lu87_5 = data.timeSeriesList(data.Intermediate, {'Mass': '87.5'})[0].data()
    SrRbYb87 = data.timeSeriesByMass(data.Intermediate, 87, 0.1).data()
    Yb86_5 = data.timeSeriesByMass(data.Intermediate, 86.5, 0.1).data()
    SrCaArYb86 = data.timeSeriesByMass(data.Intermediate, 86, 0.1).data()
    RbYbEr85 = data.timeSeriesByMass(data.Intermediate, 85, 0.1).data()
    #Tm84_5 = data.timeSeriesByMass(data.Intermediate, 84.5, 0.1).data()
    SrCaArErYb84 = data.timeSeriesByMass(data.Intermediate, 84, 0.1).data()
    Er83_5 = data.timeSeriesByMass(data.Intermediate, 83.5, 0.1).data()
    CaArEr83 = data.timeSeriesByMass(data.Intermediate, 83, 0.1).data()
    CaArErDy82 = data.timeSeriesByMass(data.Intermediate, 82, 0.1).data()

    '''
    Note that the following is not baseline subtracted as it is intended to be used as a
    reference of the raw counts
    '''
    Raw83 = data.timeSeriesByMass(data.Input, 83, 0.1).data()

    drs.message("Subtracting interferences...")
    drs.progress(40)

    Sr88_86_reference = settings['Sr88_86_reference']
    Dy_Er = settings['Dy_Er']
    Lu_Yb = settings['Lu_Yb']

    PFract = (np.log(Sr88_86_reference / (SrCaArYbLu88/SrCaArYb86))) / (np.log(87.9056/85.9093))*mask

    # The following equations subtract CaAr using the signal on 82, itself corrected for REE2+ using the
    # Er83_5 signal and canonical Dy/Er (because we don't measure Dy on a half mass).
    # The default value of Dy_Er is 1.8
    Er82 = (Er83_5 * 1.601 / 22.869) / np.power((81.96460 / 83.46619), PFract)
    Dy82 = Er82 / (1.601 / 100.) * Dy_Er * (28.260 / 100.)
    CaAr82 = CaArErDy82 - Er82 - Dy82

    # Now correct mass 86 for CaAr using CaAr82, and correct for 172Yb2+ using 173Yb2+ (measured on mass 86.5)
    CaAr86 = (CaAr82 * .004 / .647) / np.power((85.9160721 / 81.9210049), PFract)
    Yb86 = (Yb86_5 * 21.68 / 16.103) / np.power((85.968195 / 86.46911), PFract)
    PSr86 = SrCaArYb86 - CaAr86 - Yb86

    # Now correct mass 88 for CaAr using CaAr82, and correct for 176Yb2+ using 173Yb2+ (measured on mass 86.5)
    # and 176Lu2+ using 173Yb2+ and applying the same Dy_Er ratio
    CaAr88 = (CaAr82 * .187 / .647) / np.power((87.9149151 / 81.9210049), PFract)
    Yb88 = (Yb86_5 * 12.996 / 16.103) / np.power((87.97129 / 86.46911), PFract)
    Lu88 = Yb88 / (12.996 / 100.) * Lu_Yb * (2.59 / 100.)
    PSr88 = SrCaArYbLu88 - CaAr88 - Yb88 - Lu88

    # Use these CaAr and REE stripped Sr 86 and Sr 88 values to calculate a refined fractionation factor
    BetaSr = (np.log(Sr88_86_reference / (PSr88 / PSr86))) / (np.log(87.9056/85.9093))

    # You might notice that a lot of the following equations look the same as those above,
    # but just using `BetaSr` fractionation factor instead of `PFract` fractionation factor

    # Calculate Rb fractionation factor, with optional adjustment
    BetaRb = BetaSr + settings['RbBias']

    # Correct mass 85 for 170Yb2+ using 173Yb2+ (measured on mass 86.5) and 170Er2+ using 167Er2+ (measured on mass 83.5)
    Yb85 =  Yb86_5 * 2.982 / 16.103 / np.power((84.967385 / 86.46911), BetaSr)
    Er85 = Er83_5 * 14.910 / 22.869 / np.power((84.96774 / 83.46619), BetaSr)
    Rb85_c = RbYbEr85 - Yb85 - Er85
    # Calculate Rb on mass 87
    Rb87_85_reference = settings['Rb87_85_reference']
    Rb87 = Rb85_c * Rb87_85_reference / np.power((86.90918 / 84.9118), BetaRb)

    # Subtract this Rb87 amount from the 87 beam, along with 174Yb2+ (calculated from
    # 173Yb2++, on mass 86.5), to get Sr87
    Yb87 = Yb86_5 * 32.026 / 16.103 / np.power((86.96943 / 86.46911), BetaSr)
    Sr87_c = SrRbYb87 - Rb87 - Yb87

    # Calculate unique CaAr fractionation factor relative to Sr fract.
    BetaCaAr = BetaSr + settings['CaArBias']

    # Then determine Sr84 by removing CaAr84, Yb84 (168Yb2+) (from Yb86.5), and Er84 (168Er2+) (from Er83.5)
    Yb84 = Yb86_5 * 0.123 / 16.103 / np.power((83.96694 / 86.46911), BetaSr)
    Er84 = Er83_5 * 26.978 / 22.869 / np.power((83.96619 / 83.46619), BetaSr)
    CaAr84 = CaAr82 * 2.086/0.647 / np.power((83.917989 / 81.921122), BetaCaAr)
    Sr84_c = SrCaArErYb84 - CaAr84 - Er84 - Yb84

    # Then determine final Sr88 by removing CaAr88, Yb88 (176Yb2+) (from Yb86.5), and Lu88 (176Lu2+) (from Yb86.5)
    # and applying the Dy_Er factor * 2.59
    CaAr88 = CaAr82 * .187 / .647 / np.power((87.9149151 / 81.9210049), BetaSr)
    Yb88 = Yb86_5 * 12.996 / 16.103 / np.power((87.97129 / 86.46911), BetaSr)
    Lu88 = Yb88 / (12.996 / 100.) * Lu_Yb * (2.59 / 100.)
    Sr88_c = SrCaArYbLu88 - CaAr88 - Yb88 - Lu88

    # Calculate a final CaAr and REE corrected 86 beam
    CaAr86 = CaAr82 * .004 / .647 / np.power((85.9160721 / 81.9210049), BetaSr)
    Yb86 = Yb86_5 * 21.68 / 16.103 / np.power((85.968195 / 86.46911), BetaSr)
    Sr86_c = SrCaArYb86 - CaAr86 - Yb86

    # Gather up intermediate channels and add them as time series:
    # NOTE: adding _c to some channel names to avoid clashes with input channel names
    int_channel_names = ['PFract', 'CaAr82', 'PSr86', 'PSr88', 'BetaSr', 'BetaRb', 'BetaCaAr',
                         'Sr84_c', 'Sr86_c', 'Sr87_c', 'Sr88_c', 'Yb84', 'Yb85', 'Yb86', 'Yb87', 'Yb88',
                         'Er84', 'Er85', 'CaAr84', 'CaAr86', 'CaAr88', 'Lu88', 'Rb85_c', 'Rb87', 'Er82']

    int_channels = [PFract, CaAr82, PSr86, PSr88, BetaSr, BetaRb, BetaCaAr,
                         Sr84_c, Sr86_c, Sr87_c, Sr88_c, Yb84, Yb85, Yb86, Yb87, Yb88,
                         Er84, Er85, CaAr84, CaAr86, CaAr88, Lu88, Rb85_c, Rb87, Er82]

    for name, channel in zip(int_channel_names, int_channels):
        data.createTimeSeries(name, data.Intermediate, indexChannel.time(), channel)

    drs.message("Calculating further interference corrections...")
    drs.progress(60)

    Sr8786_Uncorr = (SrRbYb87 / SrCaArYb86) * np.power((86.9089 / 85.9093), BetaSr) * mask
    Sr8786_Corr = (Sr87_c / Sr86_c) * np.power((86.9089 / 85.9093), BetaSr) * mask
    Rb87Sr86_Corr = (Rb87 / Sr86_c) * np.power((86.9089 / 85.9093), BetaSr)
    Sr8486_Uncorr = (SrCaArErYb84 / SrCaArYb86) * np.power((83.9134 / 85.9093), BetaSr)
    Sr8486_Corr = (Sr84_c / Sr86_c) * np.power((83.9134 / 85.9093), BetaSr)
    Sr8488_Uncorr = (SrCaArErYb84 / SrCaArYbLu88) * np.power((83.9134 / 87.9056), BetaSr)
    Sr8488_Corr = (Sr84_c / Sr88_c) * np.power((83.9134 / 87.9056), BetaSr)
    Rb87asPPM = (Rb87 / SrRbYb87) * 1000000
    CaArErYb84asPPM = (CaAr84 + Er84 + Yb84) / SrCaArErYb84 * 100000
    TotalSrBeam = Sr88_c + Sr84_c + Sr86_c + Sr87_c

    # Gather up intermediate channels and add them as time series (again):
    int_channel_names = ['Sr8786_Uncorr', 'Sr8786_Corr', 'Rb87Sr86_Corr', 'Sr8486_Uncorr', 'Sr8486_Corr',
                        'Sr8488_Uncorr', 'Sr8488_Corr', 'Rb87asPPM', 'CaArErYb84asPPM', 'TotalSrBeam']
    int_channels = [Sr8786_Uncorr, Sr8786_Corr, Rb87Sr86_Corr, Sr8486_Uncorr, Sr8486_Corr,
                        Sr8488_Uncorr, Sr8488_Corr, Rb87asPPM, CaArErYb84asPPM, TotalSrBeam]

    for name, channel in zip(int_channel_names, int_channels):
        data.createTimeSeries(name, data.Intermediate, indexChannel.time(), channel)

    # The following equations subtract CaAr using the signal on 83, itself corrected for REE2+ using the Er83_5 signal.
    Er83 = Er83_5 * 33.503 / 22.869 / np.power((82.96514975 / 83.46619), PFract)
    CaAr83 = CaArEr83 - Er83

    # Use the preliminary fract in stripping CaAr and REE from Sr 86
    CaAr86_b = CaAr83 * .004 / .135 / np.power((85.9160721 / 82.921150), PFract)
    Yb86_b = Yb86_5 * 21.68 / 16.103 / np.power((85.968195 / 86.46911), PFract)
    PSrCaArYb86_b = SrCaArYb86 - CaAr86_b - Yb86_b

    # and Sr 88
    CaAr88_b = CaAr83 * .187 / .135 / np.power((87.9149151 / 82.921150), PFract)
    Yb88_b = Yb86_5 * 12.996 / 16.103 / np.power((87.97129 / 86.46911), PFract)
    Lu88_b = Yb88_b / (12.996 / 100.) * Lu_Yb * (2.59 / 100.)
    PSrCaArYbLu88_b = SrCaArYbLu88 - CaAr88_b - Yb88_b - Lu88_b

    # Use these CaAr and REE stripped Sr 86 and Sr 88 values to calculate a refined fractionation factor
    BetaSr_b = (np.log(Sr88_86_reference/(PSrCaArYbLu88_b / PSrCaArYb86_b))) / (np.log(87.9056/85.9093))
    BetaRb_b = BetaSr_b + RbBias

    # Use the Rb fract to calculate the amount of Rb on mass 87
    # Correct mass 85 for 170Yb2+ using 173Yb2+ (measured on mass 86.5) and 170Er2+ using 167Er2+ (measured on mass 83.5)
    Yb85_b =  Yb86_5 * 2.982 / 16.103 / np.power((84.967385 / 86.46911), BetaSr_b)
    Er85_b = Er83_5 * 14.910 / 22.869 / np.power((84.96774 / 83.46619), BetaSr_b)
    Rb85_b = RbYbEr85 - Yb85_b - Er85_b
    # Calculate Rb on mass 87
    Rb87_b = Rb85_b * Rb87_85_reference / np.power((86.90918 / 84.9118), BetaRb_b)

    # Subtract this Rb87 amount from the 87 beam, along with 174Yb2+ (calculated from
    # 173Yb2++, on mass 86.5), to get Sr87
    Yb87_b = Yb86_5 * 32.026 / 16.103 / np.power((86.96943 / 86.46911), BetaSr_b)
    Sr87_b = SrRbYb87 - Rb87_b - Yb87_b

    # Allows for a modification of the CaAr fractionation factor relative to the Sr fract. We have never used anything but 1 (i.e. no modification)
    BetaCaAr_b = BetaSr_b + CaArBias

    # Then determine Sr84 by removing CaAr84, Yb84 (168Yb2+) (from Yb86.5), and Er84 (168Er2+) (from Er83.5)
    Yb84_b = Yb86_5 * 0.123 / 16.103 / np.power((83.96694 / 86.46911), BetaSr_b)
    Er84_b = Er83_5 * 26.978 / 22.869 / np.power((83.96619 / 83.46619), BetaSr_b)
    CaAr84_b = CaAr83 * 2.086/0.135 / np.power((83.917989 / 82.921150), BetaCaAr_b)
    CaArErYb84_b = CaAr84_b + Er84_b + Yb84_b
    Sr84_b = SrCaArErYb84 - CaAr84_b - Er84_b - Yb84_b

    # Then determine final Sr88 by removing CaAr88, Yb88 (176Yb2+) (from Yb86.5), and Lu88 (176Lu2+) (from Yb86.5)
    # and applying the Dy_Er factor * 2.59
    CaAr88_b = CaAr83 * .187 / .135 / np.power((87.9149151 / 82.921150), BetaSr_b)
    Yb88_b = Yb86_5 * 12.996 / 16.103 / np.power((87.97129 / 86.46911), BetaSr_b)
    Lu88_b = Yb88_b / (12.996 / 100.) * Lu_Yb * (2.59 / 100.)
    Sr88_b = SrCaArYbLu88 - CaAr88_b - Yb88_b - Lu88_b

    # Calculate a final CaAr and REE corrected 86 beam
    CaAr86_b = CaAr83 * .004 / .135 / np.power((85.9160721 / 82.921150), BetaSr_b)
    Yb86_b = Yb86_5 * 21.68 / 16.103 / np.power((85.968195 / 86.46911), BetaSr_b)
    Sr86_b = SrCaArYb86 - CaAr86_b - Yb86_b

    Sr8786_Uncorr_b = (SrRbYb87 / SrCaArYb86) * np.power((86.9089 / 85.9093), BetaSr_b)
    Sr8786_Corr_b = (Sr87_b / Sr86_b) * np.power((86.9089 / 85.9093), BetaSr_b)
    Rb87Sr86_Corr_b = (Rb87_b / Sr86_b) * np.power((86.9089 / 85.9093), BetaSr_b)
    Sr8486_Uncorr_b = (SrCaArErYb84 / SrCaArYb86) * np.power((83.9134 / 85.9093), BetaSr_b)
    Sr8486_Corr_b = (Sr84_b / Sr86_b) * np.power((83.9134 / 85.9093), BetaSr_b)
    Sr8488_Uncorr_b = (SrCaArErYb84 / SrCaArYbLu88) * np.power((83.9134 / 87.9056), BetaSr_b)
    Sr8488_Corr_b = (Sr84_b / Sr88_b) * np.power((83.9134 / 87.9056), BetaSr_b)
    Rb87asPPM_b = (Rb87_b / SrRbYb87) * 1000000
    CaArErYb84asPPM_b = (CaArErYb84_b / SrCaArErYb84) * 100000
    TotalSrBeam_b = Sr88_b + Sr84_b + Sr86_b + Sr87_b

    # Gather up intermediate channels and add them as time series:
    int_channel_names = ['Er83', 'CaAr83', 'Sr84_b', 'Sr86_b', 'Sr87_b', 'Sr88_b', 'BetaSr_b', 'Rb87_b',
                        'Sr8786_Uncorr_b', 'Sr8786_Corr_b', 'Rb87Sr86_Corr_b',
                        'Sr8486_Uncorr_b', 'Sr8486_Corr_b', 'Sr8488_Uncorr_b',
                        'Sr8488_Corr_b', 'Rb87asPPM_b', 'CaArErYb84asPPM_b', 'TotalSrBeam_b']

    int_channels = [Er83, CaAr83, Sr84_b, Sr86_b, Sr87_b, Sr88_b, BetaSr_b, Rb87_b,
                    Sr8786_Uncorr_b, Sr8786_Corr_b, Rb87Sr86_Corr_b,
                    Sr8486_Uncorr_b, Sr8486_Corr_b, Sr8488_Uncorr_b,
                    Sr8488_Corr_b, Rb87asPPM_b, CaArErYb84asPPM_b, TotalSrBeam_b]

    for name, channel in zip(int_channel_names, int_channels):
        data.createTimeSeries(name, data.Intermediate, indexChannel.time(), channel)

    drs.message("Calculating reference material corrected results...")
    drs.progress(70)
    
    data.updateResults()
    
    # Now check if there are selections for the reference standard, and if so, generate standard-normalised ratios
    try:
        StdSpline_Sr87_86 = data.spline(rmName, "Sr8786_Corr").data()
    except:
        IoLog.errorWithOrigin("The Combined Sr DRS requires Ref Material selections to proceed.", "BP_Sr_DRS")
        drs.message("DRS did not finish. Please check Messages")
        drs.progress(100)
        drs.finished()
        return
    # And check that we can get the RM 87/86 value    
    try:
        StdValue_Sr87_86 = data.referenceMaterialData(rmName)["87Sr/86Sr"].value()
    except Exception as e:
        IoLog.errorWithOrigin(f"Could not get the 87Sr_86Sr value from the RM file: {e}", "BP_DRS")
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


    # CaPO Plot and Correction
    if len(settings['CaPO_RMs']) > 0:
        print("Here are the RMs for CaPO correction:", settings['CaPO_RMs'])
        PLOT.clearGraphs()

        # Get Sr signal and deviations for chosen RMs:
        Sr8786 = data.timeSeries('StdCorr_Sr8786')
        totalSr = data.timeSeries('TotalSrBeam')

        #collect up all deviations and Sr intensities here
        devs_main = []
        sigs_main = []

        for i, rm in enumerate(settings['CaPO_RMs']):
            # Also have list for this group only. For plotting below
            devs = []
            sigs = []
            try:
                sg = data.selectionGroup(rm)
            except RuntimeError as err:
                IoLog.warning(f"Could not get the selections for group {rm}: {err}")
                continue

            try:
                true = data.referenceMaterialData(rm)['87Sr/86Sr'].value()
            except KeyError as err:
                IoLog.warning(f"Could not get the \'87Sr/86Sr\' value for {rm}: {err}")
                continue

            for sel in sg.selections():
                meas = data.result(sel, Sr8786).value()
                devs.append(meas/true)
                devs_main.append(meas/true)
                sigs.append(data.result(sel, totalSr).value())
                sigs_main.append(data.result(sel, totalSr).value())

            g = PLOT.addGraph()
            g.setName(rm)
            g.setLineStyle('lsNone')
            g.setScatterStyle('ssDisc', 6.0, PLOT_COLORS[i])
            g.setData(np.array(sigs), np.array(devs))

        # Now add fit to all data
        def fitFunc(t, a, b):
            return a + b*t

        # Fit a curve to the raw downhole data
        sigs_array = np.array(sigs_main)
        devs_array = np.array(devs_main)
        params, cov = curve_fit(fitFunc, sigs_array, devs_array, ftol=1e-5)
        print(f"Here are the fit params: Slope: {params[1]}, Intercept: {params[0]}")

        # Show the fit +/- 10% of the spread in signals
        sigs_range = sigs_array.max() - sigs_array.min()
        fit_x = np.linspace(
            sigs_array.min() - sigs_range * 0.1,
            sigs_array.max() + sigs_range * 0.1 ,
            10
        )
        fit_y = params[1]*fit_x + params[0]

        g = PLOT.addGraph()
        g.setName("fit")
        g.setData(fit_x, fit_y)

        PLOT.left().setLabel('87Sr/86Sr Deviation (meas/true)')
        PLOT.bottom().setLabel('Total Sr (V)')
        PLOT.setToolsVisible(False)
        PLOT.rescaleAxes()
        PLOT.replot()

    data.registerAssociatedResult("Corr_Sr8786_Sr88_V", getSelRatioIntensityCorr)

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
    drs.setSetting("Lu_Yb", 7.0)
    drs.setSetting("RbBias", 0.)
    drs.setSetting("CaArBias", 0.)
    drs.setSetting("Sr88_86_reference", 8.37520938)  #Konter & Storm (2014)
    drs.setSetting("Rb87_85_reference", 0.385710)    #Konter & Storm (2014)
    drs.setSetting("CaPO_RMs", [])

    drs.setSetting("PropagateError", False)

    settings = drs.settings()

    indexComboBox = QtGui.QComboBox(widget)
    indexComboBox.addItems(timeSeriesNames)
    indexComboBox.setCurrentText(settings["IndexChannel"])
    indexComboBox.textActivated.connect(lambda t: drs.setSetting("IndexChannel", t))
    formLayout.addRow("Index channel", indexComboBox)

    rmComboBox = QtGui.QComboBox(widget)
    rmComboBox.addItems(rmNames)
    if settings["ReferenceMaterial"] in rmNames:
        rmComboBox.setCurrentText(settings["ReferenceMaterial"])
        drs.setSetting("ReferenceMaterial", settings["ReferenceMaterial"])

    rmComboBox.textActivated.connect(lambda t: drs.setSetting("ReferenceMaterial", t))
    formLayout.addRow("Reference material", rmComboBox)

    verticalSpacer = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
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

    verticalSpacer2 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
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

    luYbLineEdit = QtGui.QLineEdit(widget)
    luYbLineEdit.setText(settings["Lu_Yb"])
    luYbLineEdit.textChanged.connect(lambda t: drs.setSetting("Lu_Yb", float(t)))
    formLayout.addRow("Lu/Yb ratio", luYbLineEdit)

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

    setExtButton = QtGui.QToolButton(widget)
    rmMenu = ReferenceMaterialsMenu(setExtButton)
    rmMenu.rmsChanged.connect(lambda l: drs.setSetting("CaPO_RMs", l))
    setExtButton.setMenu(rmMenu)
    setExtButton.setPopupMode(QtGui.QToolButton.InstantPopup)
    setExtButton.setIcon(CUI().icon('trophy'))
    setExtButton.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
    formLayout.addRow("RMs for CaPO correction:", setExtButton)

    propCheckBox = QtGui.QCheckBox(widget)
    propCheckBox.setChecked(settings["PropagateError"])
    propCheckBox.toggled.connect(lambda t: drs.setSetting("PropagateError", bool(t)))
    formLayout.addRow("Propagate Errors?", propCheckBox)

    # Plot for CaPO correction
    formLayout.addRow('CaPO Correction fit', PLOT)

    # Restore settings
    try:
        settings = drs.settings()
        print('Restoring settings...')
        print(settings)
        rmComboBox.setCurrentText(settings["ReferenceMaterial"])
    except KeyError:
        pass

    drs.setSettingsWidget(widget)
