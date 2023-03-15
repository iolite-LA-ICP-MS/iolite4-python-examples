#/ Type: DRS
#/ Name: Sr Isotopes Universal
#/ Authors: Bence Paul, Joe Petrus, Graham Hagen-Peter, author(s) of Sr_isotopes_Total_NIGL.ipf, and various authors of previous Iolite Sr isotope DRS
#/ Description: A Sr isotopes DRS that can correct for different combinations of interferences
#/ References: Mulder et al. (2023) Geostandards and Geolanalytical Research
#/ Version: 2.0
#/ Contact: support@iolite-software.com

from iolite import QtGui
from iolite import QtCore
from iolite.Qt import Qt, QColor
from iolite.ui import CommonUIPyInterface as CUI
from iolite.ui import IolitePlotPyInterface as Plot
from iolite.ui import IolitePlotSettingsDialog as PlotSettings
from iolite.QtGui import QAction, QPen
from iolite.types import Result

from scipy.optimize import curve_fit
import numpy as np
import itertools
import math
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

    # Add annotation to show fit parameters
    ann = PLOT.annotate('', 0.01, 0.951, 'ptAxisRectRatio', Qt.AlignLeft | Qt.AlignBottom)
    ann.visible = False

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
Defining the colors to use in the CaPO plot here. By default only has 10 colors:
'''
PLOT_COLORS = [
    QColor(239, 71, 111),   # Carmine Red
    QColor(255, 209, 102),  # Orange Yellow
    QColor(6, 214, 160),    # Green
    QColor(17, 138, 178),   # Blue
    QColor(7, 59, 76),      # Midnight Green
    QColor(204, 245, 172),  # Light Green
    QColor(128, 138, 159),  # Roman Silver
    QColor(61, 43, 86),     # Russian Purple
    QColor(247, 169, 168),  # Pastel Pink
    QColor(207, 212, 197),  # Bone
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
        IoLog.errorWithOrigin("Index channel not found. DRS cannot continue", "Sr Isotopes Universal DRS")
        drs.message("DRS did not finish. Please check Messages")
        drs.progress(100)
        drs.finished()
        return

    indexChannel = data.timeSeries(settings["IndexChannel"])
    Sr_rmName = settings["ReferenceMaterial"]
    maskOption = settings["Mask"]
    maskChannel = data.timeSeries(settings["MaskChannel"])
    cutoff = settings["MaskCutoff"]
    trim = settings["MaskTrim"]

    Rb_Beta_adjust = settings["Rb_Beta_adjust"]
    RbBias = settings["RbBias"]
    Ref_Rb_Beta = settings["ReferenceMaterialRb"]
    Rb_Sr_elemental = settings["Rb_Sr_elemental"]
    Rb_Sr_fract = settings["Rb_Sr_fract"]
    Ref_Rb_Sr_elemental = settings["ReferenceMaterialRbSr"]
    REE_subtract = settings["REE_subtract"]
    REEBias = settings["REEBias"]
    Dy_Er = settings["Dy_Er"]
    Lu_Yb = settings["Lu_Yb"]
    CaAr_CaCa_subtract = settings["CaAr_CaCa_subtract"]
    CaAr_83 = settings["CaAr_83"]
    CaArBias = settings["CaArBias"]
    NaNi_CaAlO_subtract = settings["NaNi_CaAlO_subtract"]
    NaNiBias = settings["NaNiBias"]
    ProportionCaAr = settings["ProportionCaAr"]
    ProportionNaNi = settings["ProportionNaNi"]
    Age = settings["Age"]
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
    #IoLog.debug("Corrections selected: " + corrections)

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

    '''
    The following channels names include the REE interferences but
    we won't be correcting for REE unless the user selects CaAr+REE.
    The names are common to both correction profiles
    '''

#########################
    try:
        total88 = data.timeSeriesList(data.Intermediate, {'Mass': '88'})[0].data()
        total87 = data.timeSeriesList(data.Intermediate, {'Mass': '87'})[0].data()
        total86 = data.timeSeriesList(data.Intermediate, {'Mass': '86'})[0].data()
        total85 = data.timeSeriesList(data.Intermediate, {'Mass': '85'})[0].data()
        total84 = data.timeSeriesList(data.Intermediate, {'Mass': '84'})[0].data()
    
    except:
        IoLog.error("This DRS requires data for mass 88, 87, 86, 85, and 84 for the basic calculations.")
        drs.message("DRS did not finish. Please check Messages")
        drs.progress(100)
        drs.finished()
        return

    try:
        Y89 = data.timeSeriesList(data.Intermediate, {'Mass': '89'})[0].data()
    except:
        pass

    try:
        total87_5 = data.timeSeriesList(data.Intermediate, {'Mass': '87.5'})[0].data()
    except:
        pass

    try:
        total86_5 = data.timeSeriesList(data.Intermediate, {'Mass': '86.5'})[0].data()
    except:
        pass

    try:
        total84_5 = data.timeSeriesList(data.Intermediate, {'Mass': '84.5'})[0].data()
    except:
        pass

    try:
        total83_5 = data.timeSeriesList(data.Intermediate, {'Mass': '83.5'})[0].data()
    except:
        pass

    try:
        total83 = data.timeSeriesList(data.Intermediate, {'Mass': '83'})[0].data()
    except:
        pass

    try:
        total82_5 = data.timeSeriesList(data.Intermediate, {'Mass': '82.5'})[0].data()
    except:
        pass

    try:
        total82 = data.timeSeriesList(data.Intermediate, {'Mass': '82'})[0].data()
    except:  
        pass

    try:
        total81_5 = data.timeSeriesList(data.Intermediate, {'Mass': '81.5'})[0].data()
    except:
        pass  

    if ProportionCaAr < 0 or ProportionCaAr >1:
        IoLog.error("The proportion CaAr must be between 0 and 1.")
        drs.message("DRS did not finish. Please check Messages")
        drs.progress(100)
        drs.finished()
        return

    if ProportionNaNi < 0 or ProportionNaNi >1:
        IoLog.error("The proportion NaNi must be between 0 and 1.")
        drs.message("DRS did not finish. Please check Messages")
        drs.progress(100)
        drs.finished()
        return


    # UNIVERSAL DRS interference subtraction combinations here

    #True masses from NIST; Monoatomic isotopice abundances from IUPAC subcommittee for Isotopic Abundance Measurements (1999) compiled by U. of Alberta; Polyatomic isotopologue abundances calculated from Scientific Instrument Services "Isotope Distribution Calculator" (sisweb.com/mstools/isotope.htm)
    Sr84mass = 83.9134191
    Sr86mass = 85.9092606
    Sr87mass = 86.9088775
    Sr88mass = 87.9056125
    Rb85mass = 84.9117897379
    Rb87mass = 86.9091805310


    Sr88_86_reference = 8.37520938 #Reciprocal of 86Sr/88Sr = 0.1194 from Konter & Storm (2014, Chem. Geol)

    Rb87_85_reference = 0.38571 

    PFract = (np.log(Sr88_86_reference/(total88/total86)))/(np.log(Sr88mass/Sr86mass))##calculate preliminary fractionation factor

    drs.message("Subtracting interferences...")
    drs.progress(40)

    ##REE plus/minus Ca-Ar-CaCa and/or NaNi-CaAlO subtractions

    if REE_subtract:
        residual_86_REE = total86 - ((total86_5 * 21.83 /16.103) / np.power((85.968195 / 86.46911), PFract*REEBias))

        try:
            residual_88_REE = total88 - ((total86_5 * 12.76 /16.103) / np.power((87.971284 / 86.46911), PFract*REEBias))-((total87_5 * 2.59 /97.41) / np.power((87.97134 / 87.47039), PFract*REEBias))
        except:
            IoLog.error("There is no 87.5 channel, so 176Lu++ will be calculated using the monitored 173Yb++ on 86.5 and a user-defined Lu/Yb. Click OK to proceed.")
            residual_88_REE = total88 - ((total86_5 * 12.76 /16.103) / np.power((87.971284 / 86.46911), PFract*REEBias))- ((total86_5 * 12.76 /16.103) / np.power((87.971284 / 86.46911), PFract*REEBias)*(1/0.1276) * Lu_Yb * 0.0259)

        PFract_REE = (np.log(Sr88_86_reference/(residual_88_REE/residual_86_REE)))/(np.log(Sr88mass/Sr86mass))

        try:
            residual_82_REE = total82 - (total83_5 * 1.61 /22.93) / np.power((81.96460 / 83.46619), PFract_REE*REEBias) - ((total81_5 * 28.18 /24.90) / np.power((81.96459 / 81.46436), PFract_REE*REEBias))
        except:
            try:
                IoLog.error("There is no 81.5 channel, so 164Dy++ will be calculated using the monitored 167Er++ on 83.5 and a user-defined Dy/Er. Click OK to proceed.")
                residual_82_REE = total82 - (total83_5 * 1.61 /22.93) / np.power((81.96460 / 83.46619), PFract_REE*REEBias) - ((total83_5 * 1.61 /22.93) / np.power((81.96460 / 83.46619), PFract_REE*REEBias)*(1/0.0161) * Dy_Er * 0.2818)
            except:
                pass
            pass

        try:
            residual_83_REE = total83 -  (total83_5 * 33.61 /22.93) / np.power((82.96514975 / 83.46619), PFract_REE*REEBias)
        except:
            pass


        residual_84_REE = total84 - ((total83_5 * 26.78 /22.93)/ np.power((83.96619 / 83.46619), PFract_REE*REEBias)-(total86_5 * 0.13 /16.13) / np.power((83.96694 / 86.46911), PFract_REE*REEBias))
        residual_85_REE = total85 - (total86_5 * 3.04 / 16.13) / np.power((84.967385 / 86.46911), PFract_REE*REEBias) - (total83_5 * 14.93 /22.93) / np.power((84.96774 / 83.46619), PFract_REE*REEBias)
        residual_87_REE = total87 - ((total86_5 * 31.83 /16.13) / np.power((86.96943 / 86.46911), PFract_REE*REEBias))

        Er_Sr_ppm = (total83_5/0.2293)/(total88/0.8258) * 1e6
        Yb_Sr_ppm = (total86_5/0.1613)/(total88/0.8258) * 1e6


        if CaAr_CaCa_subtract:
            if CaAr_83:
                residual_86_CaAr = residual_86_REE - (((ProportionCaAr)*residual_83_REE * 0.0039 /0.139) / np.power((85.91607 / 82.92115), PFract_REE*CaArBias)) - (((1-ProportionCaAr)*residual_83_REE * 0.035 /0.272) / np.power((85.91411 / 82.92136), PFract_REE*CaArBias))
                residual_88_CaAr = residual_88_REE - (((ProportionCaAr)*residual_83_REE * 0.189 /0.139) / np.power((87.91491 / 82.92115), PFract_REE*CaArBias)) - (((1-ProportionCaAr)*residual_83_REE * 0.412 /0.272) / np.power((87.91512 / 82.92136), PFract_REE*CaArBias))

                PFract_CaAr = (np.log(Sr88_86_reference/(residual_88_CaAr/residual_86_CaAr)))/(np.log(Sr88mass/Sr86mass))

                try:
                    residual_82_CaAr = residual_82_REE -  (((ProportionCaAr)*residual_83_REE * 0.649 /0.139) / np.power((81.921 / 82.92115), PFract_CaAr*CaArBias)) - (((1-ProportionCaAr)*residual_83_REE * 1.260 /0.272) / np.power((81.92121 / 82.92136), PFract_CaAr*CaArBias))
                except:
                    pass

                residual_84_CaAr = residual_84_REE - (((ProportionCaAr)*residual_83_REE * 2.078 /0.139) / np.power((83.92008 / 82.92115), PFract_CaAr*CaArBias)) - (((1-ProportionCaAr)*residual_83_REE * 4.048 /0.272) / np.power((83.91808 / 82.92136), PFract_CaAr*CaArBias))
                residual_85_CaAr = residual_85_REE - (((1-ProportionCaAr)*residual_83_REE * 0.0019 /0.272) / np.power((84.91739 / 82.92136), PFract_CaAr*CaArBias))
                residual_87_CaAr = residual_87_REE - (((1-ProportionCaAr)*residual_83_REE * 0.0056 /0.272) / np.power((86.91426 / 82.92136), PFract_CaAr*CaArBias))

            else:
                residual_86_CaAr = residual_86_REE - (((ProportionCaAr)*residual_82_REE * 0.0039 /0.649) / np.power((85.91607 / 81.921), PFract_REE*CaArBias)) - (((1-ProportionCaAr)*residual_82_REE * 0.035 /1.260) / np.power((85.91411 / 81.92121), PFract_REE*CaArBias))
                residual_88_CaAr = residual_88_REE - (((ProportionCaAr)*residual_82_REE * 0.189 /0.649) / np.power((87.91491 / 81.921), PFract_REE*CaArBias)) - (((1-ProportionCaAr)*residual_82_REE * 0.412 /1.260) / np.power((87.91512 / 81.92121), PFract_REE*CaArBias))

                PFract_CaAr = (np.log(Sr88_86_reference/(residual_88_CaAr/residual_86_CaAr)))/(np.log(Sr88mass/Sr86mass))

                try:
                    residual_83_CaAr = residual_83_REE -  (((ProportionCaAr)*residual_82_REE * 0.139 /0.649) / np.power((82.92115 / 81.921), PFract_CaAr*CaArBias)) - (((1-ProportionCaAr)*residual_82_REE * 0.272 /1.260) / np.power((82.92136 / 81.92121), PFract_CaAr*CaArBias))
                except:
                    pass

                residual_84_CaAr = residual_84_REE - (((ProportionCaAr)*residual_82_REE * 2.078 /0.649) / np.power((83.92008 / 81.921), PFract_CaAr*CaArBias)) - (((1-ProportionCaAr)*residual_82_REE * 4.048 /1.260) / np.power((83.91808 / 81.92121), PFract_CaAr*CaArBias))
                residual_85_CaAr = residual_85_REE - (((1-ProportionCaAr)*residual_82_REE * 0.0019 /1.260) / np.power((84.91739 / 81.92121), PFract_CaAr*CaArBias))
                residual_87_CaAr = residual_87_REE - (((1-ProportionCaAr)*residual_82_REE * 0.0056 /1.260) / np.power((86.91426 / 81.92121), PFract_CaAr*CaArBias))

        else:
            PFract_CaAr = PFract_REE

            residual_86_CaAr = residual_86_REE
            residual_88_CaAr = residual_88_REE

            try:
                residual_83_CaAr = residual_83_REE
            except:
                pass

            residual_84_CaAr = residual_84_REE
            residual_85_CaAr = residual_85_REE
            residual_87_CaAr = residual_87_REE

        if NaNi_CaAlO_subtract:
            if CaAr_83:
                IoLog.error("Using mass 83 for CaAr-CaCa peak stripping assumes that there is no NaNi-CaAlO on mass 83. Use mass 82 for CaAr-CaCa peak stripping if you want to do these corrections independently.")
                drs.message("DRS did not finish. Please check Messages")
                drs.progress(100)
                drs.finished()
                return

            residual_86_NaNi = residual_86_CaAr - (((1-ProportionNaNi)*residual_83_CaAr * 0.139 /96.701) / np.power((85.93523 / 82.93905), PFract_CaAr*NaNiBias))
            residual_88_NaNi = residual_88_CaAr - (((1-ProportionNaNi)*residual_83_CaAr * 0.001 /96.701) / np.power((87.93616 / 82.93905), PFract_CaAr*NaNiBias))
            PFract_NaNi = (np.log(Sr88_86_reference/(residual_88_NaNi/residual_86_NaNi)))/(np.log(Sr88mass/Sr86mass))

            residual_84_NaNi = residual_84_CaAr - (((ProportionNaNi)*residual_83_CaAr * 1.130 /26.100) / np.power((83.92083 / 82.92056), PFract_NaNi*NaNiBias)) - (((1-ProportionNaNi)*residual_83_CaAr * 0.039 /96.701) / np.power((83.94326 / 82.93905), PFract_NaNi*NaNiBias))
            residual_85_NaNi = residual_85_CaAr - (((ProportionNaNi)*residual_83_CaAr * 3.590 /26.100) / np.power((84.91812 / 82.92056), PFract_NaNi*NaNiBias)) - (((1-ProportionNaNi)*residual_83_CaAr * 0.841 /96.701) / np.power((84.93508 / 82.93905), PFract_NaNi*NaNiBias))
            residual_87_NaNi = residual_87_CaAr - (((ProportionNaNi)*residual_83_CaAr * 0.910 /26.100) / np.power((86.91774 / 82.92056), PFract_NaNi*NaNiBias)) - (((1-ProportionNaNi)*residual_83_CaAr * 2.082 /96.701) / np.power((86.93195 / 82.93905), PFract_NaNi*NaNiBias))
        else:
            residual_84_NaNi = residual_84_CaAr
            residual_85_NaNi = residual_85_CaAr
            residual_86_NaNi = residual_86_CaAr
            residual_87_NaNi = residual_87_CaAr
            residual_88_NaNi = residual_88_CaAr

    ##CaAr-CaCa plus/minus Na-Ni-CaAlO subtractions
    else:

        if CaAr_CaCa_subtract:

            if CaAr_83:

                residual_86_CaAr = total86 - (((ProportionCaAr)*total83 * 0.0039 /0.139) / np.power((85.91607 / 82.92115), PFract*CaArBias)) - (((1-ProportionCaAr)*total83 * 0.035 /0.272) / np.power((85.91411 / 82.92136), PFract*CaArBias))
                residual_88_CaAr = total88 - (((ProportionCaAr)*total83 * 0.189 /0.139) / np.power((87.91491 / 82.92115), PFract*CaArBias)) - (((1-ProportionCaAr)*total83 * 0.412 /0.272) / np.power((87.91512 / 82.92136), PFract*CaArBias))

                PFract_CaAr = (np.log(Sr88_86_reference/(residual_88_CaAr/residual_86_CaAr)))/(np.log(Sr88mass/Sr86mass))

                try:
                    residual_82_CaAr = total82 -  (((ProportionCaAr)*total83 * 0.649 /0.139) / np.power((81.921 / 82.92115), PFract_CaAr*CaArBias)) - (((1-ProportionCaAr)*total83 * 1.260 /0.272) / np.power((81.92121 / 82.92136), PFract_CaAr*CaArBias))
                except:
                    pass

                residual_84_CaAr = total84 - (((ProportionCaAr)*total83 * 2.078 /0.139) / np.power((83.92008 / 82.92115), PFract_CaAr*CaArBias)) - (((1-ProportionCaAr)*total83 * 4.048 /0.272) / np.power((83.91808 / 82.92136), PFract_CaAr*CaArBias))
                residual_85_CaAr = total85 - (((1-ProportionCaAr)*total83 * 0.0019 /0.272) / np.power((84.91739 / 82.92136), PFract_CaAr*CaArBias))
                residual_87_CaAr = total87 - (((1-ProportionCaAr)*total83 * 0.0056 /0.272) / np.power((86.91426 / 82.92136), PFract_CaAr*CaArBias))

            else:

                residual_86_CaAr = total86 - (((ProportionCaAr)*total82 * 0.0039 /0.649) / np.power((85.91607 / 81.921), PFract*CaArBias)) - (((1-ProportionCaAr)*total82 * 0.035 /1.260) / np.power((85.91411 / 81.92121), PFract*CaArBias))
                residual_88_CaAr = total88 - (((ProportionCaAr)*total82 * 0.189 /0.649) / np.power((87.91491 / 81.921), PFract*CaArBias)) - (((1-ProportionCaAr)*total82 * 0.412 /1.260) / np.power((87.91512 / 81.92121), PFract*CaArBias))

                PFract_CaAr = (np.log(Sr88_86_reference/(residual_88_CaAr/residual_86_CaAr)))/(np.log(Sr88mass/Sr86mass))

                try:
                    residual_83_CaAr = total83 -  (((ProportionCaAr)*total82 * 0.139 /0.649) / np.power((82.92115 / 81.921), PFract_CaAr*CaArBias)) - (((1-ProportionCaAr)*total82 * 0.272 /1.260) / np.power((82.92136 / 81.92121), PFract_CaAr*CaArBias))
                except:
                    pass

                residual_84_CaAr = total84 - (((ProportionCaAr)*total82 * 2.078 /0.649) / np.power((83.92008 / 81.921), PFract_CaAr*CaArBias)) - (((1-ProportionCaAr)*total82 * 4.048 /1.260) / np.power((83.91808 / 81.92121), PFract_CaAr*CaArBias))
                residual_85_CaAr = total85 - (((1-ProportionCaAr)*total82 * 0.0019 /1.260) / np.power((84.91739 / 81.92121), PFract_CaAr*CaArBias))
                residual_87_CaAr = total87 - (((1-ProportionCaAr)*total82 * 0.0056 /1.260) / np.power((86.91426 / 81.92121), PFract_CaAr*CaArBias))

        else:
            PFract_CaAr = PFract

            residual_86_CaAr = total86
            residual_88_CaAr = total88

            try:
                residual_83_CaAr = total83
            except:
                pass

            residual_84_CaAr = total84
            residual_85_CaAr = total85
            residual_87_CaAr = total87

        if NaNi_CaAlO_subtract:

            if CaAr_83:
                IoLog.error("Using mass 83 for CaAr-CaCa peak stripping assumes that there is no NaNi-CaAlO on mass 83. Use mass 82 for CaAr-CaCa peak stripping if you want to do these corrections independently.")
                drs.message("DRS did not finish. Please check Messages")
                drs.progress(100)
                drs.finished()
                return

            residual_86_NaNi = residual_86_CaAr - (((1-ProportionNaNi)*residual_83_CaAr * 0.139 /96.701) / np.power((85.93523 / 82.93905), PFract_CaAr*NaNiBias))
            residual_88_NaNi = residual_88_CaAr - (((1-ProportionNaNi)*residual_83_CaAr * 0.001 /96.701) / np.power((87.93616 / 82.93905), PFract_CaAr*NaNiBias))
            PFract_NaNi = (np.log(Sr88_86_reference/(residual_88_NaNi/residual_86_NaNi)))/(np.log(Sr88mass/Sr86mass))

            residual_84_NaNi = residual_84_CaAr - (((ProportionNaNi)*residual_83_CaAr * 1.130 /26.100) / np.power((83.92083 / 82.92056), PFract_NaNi*NaNiBias)) - (((1-ProportionNaNi)*residual_83_CaAr * 0.039 /96.701) / np.power((83.94326 / 82.93905), PFract_NaNi*NaNiBias))
            residual_85_NaNi = residual_85_CaAr - (((ProportionNaNi)*residual_83_CaAr * 3.590 /26.100) / np.power((84.91812 / 82.92056), PFract_NaNi*NaNiBias)) - (((1-ProportionNaNi)*residual_83_CaAr * 0.841 /96.701) / np.power((84.93508 / 82.93905), PFract_NaNi*NaNiBias))
            residual_87_NaNi = residual_87_CaAr - (((ProportionNaNi)*residual_83_CaAr * 0.910 /26.100) / np.power((86.91774 / 82.92056), PFract_NaNi*NaNiBias)) - (((1-ProportionNaNi)*residual_83_CaAr * 2.082 /96.701) / np.power((86.93195 / 82.93905), PFract_NaNi*NaNiBias))
        else:
            residual_84_NaNi = residual_84_CaAr
            residual_85_NaNi = residual_85_CaAr
            residual_86_NaNi = residual_86_CaAr
            residual_87_NaNi = residual_87_CaAr
            residual_88_NaNi = residual_88_CaAr

    #final=risidualNaNi
    Sr84_corr = residual_84_NaNi
    Rb85_corr = residual_85_NaNi
    Sr86_corr = residual_86_NaNi
    SrRb87_corr = residual_87_NaNi
    Sr88_corr = residual_88_NaNi

    BetaSr = (np.log(Sr88_86_reference/(Sr88_corr/Sr86_corr)))/(np.log(Sr88mass/Sr86mass))
    Rb87_corr = (Rb85_corr * Rb87_85_reference) / np.power((Rb87mass / Rb85mass), BetaSr)
    Sr87_corr = (SrRb87_corr - Rb87_corr)

    Sr84_86_Corr = (Sr84_corr / Sr86_corr) * np.power((Sr84mass / Sr86mass), BetaSr)
    Sr84_88_Corr = (Sr84_corr / Sr88_corr) * np.power((Sr84mass / Sr88mass), BetaSr)

    Sr87_86_Corr = (Sr87_corr / Sr86_corr) * np.power((Sr87mass / Sr86mass), BetaSr)
    Rb87_Sr86_Corr = (Rb87_corr / Sr86_corr) * np.power((Rb87mass / Sr86mass), BetaSr)

    # Gather up intermediate channels and add them as time series:
    int_channel_names = ['Sr84_corr','Rb85_corr','Sr86_corr','Sr87_corr','Sr88_corr']
    int_channel_names += ['Sr87_86_Corr','BetaSr']

    int_channels = [Sr84_corr,Rb85_corr,Sr86_corr,Sr87_corr,Sr88_corr]
    int_channels += [Sr87_86_Corr,BetaSr]

    for name, channel in zip(int_channel_names, int_channels):
        data.createTimeSeries(name, data.Intermediate, indexChannel.time(), channel)

    drs.message("Calculating reference material corrected results...")
    drs.progress(70)

    try:
        StdSpline_Sr87_86 = data.spline(Sr_rmName, "Sr87_86_Corr").data()
        StdValue_Sr87_86 = data.referenceMaterialData(Sr_rmName)["87Sr/86Sr"].value()
    except:
        IoLog.error("The Combined Sr DRS requires Sr isotope Ref Material selections to proceed.")
        drs.message("DRS did not finish. Please check Messages")
        drs.progress(100)
        drs.finished()
        return

    StdCorr_Sr87_86 = Sr87_86_Corr * StdValue_Sr87_86 / StdSpline_Sr87_86

    # Calculating Sr isotope ratio with adjusted Rb mass-bias factor (if "Rb_Beta_adjust" is selected)
    if Rb_Beta_adjust:
        try:
            RefRbValue_Sr87_86 = data.referenceMaterialData(Ref_Rb_Beta)["87Sr/86Sr"].value()
            BetaRb = np.log(((SrRb87_corr)-(RefRbValue_Sr87_86 /np.power((Sr87mass / Sr86mass), BetaSr)*Sr86_corr))/(Rb85_corr*Rb87_85_reference))/np.log(Rb85mass/Rb87mass)
            data.createTimeSeries('BetaRb', data.Intermediate, indexChannel.time(), BetaRb)
            StdSpline_BetaRb = data.spline(Ref_Rb_Beta, "BetaRb").data()
        
        except:
            IoLog.error("The Combined Sr DRS requires Rb Beta Ref Material (with reference 87Sr/86Sr) selections to proceed.")
            drs.message("DRS did not finish. Please check Messages")
            drs.progress(100)
            drs.finished()
            return

        Rb87_CorrRb = (Rb85_corr * Rb87_85_reference) / np.power((Rb87mass / Rb85mass), (StdSpline_BetaRb))
        Sr87_CorrRb = (SrRb87_corr - Rb87_CorrRb)
        Sr87_86_CorrRb = (Sr87_CorrRb / Sr86_corr) * np.power((Sr87mass / Sr86mass), BetaSr)
        Rb87_Sr86_CorrRb = (Rb87_CorrRb / Sr86_corr) * np.power((Rb87mass / Sr86mass), BetaSr)

        data.createTimeSeries('Sr87_86_CorrRb', data.Intermediate, indexChannel.time(), Sr87_86_CorrRb)

        StdSplineRb_Sr87_86 = data.spline(Sr_rmName, "Sr87_86_CorrRb").data()
        StdCorrRb_Sr87_86 = Sr87_86_CorrRb * StdValue_Sr87_86 / StdSplineRb_Sr87_86

    else:
        BetaRb = BetaSr * RbBias
        data.createTimeSeries('BetaRb', data.Intermediate, indexChannel.time(), BetaRb)

        Rb87_CorrRb = (Rb85_corr * Rb87_85_reference) / np.power((Rb87mass / Rb85mass), (BetaRb))
        Sr87_CorrRb = (SrRb87_corr - Rb87_CorrRb)
        Sr87_86_CorrRb = (Sr87_CorrRb / Sr86_corr) * np.power((Sr87mass / Sr86mass), BetaSr)
        Rb87_Sr86_CorrRb = (Rb87_CorrRb / Sr86_corr) * np.power((Rb87mass / Sr86mass), BetaSr)

        data.createTimeSeries('Sr87_86_CorrRb', data.Intermediate, indexChannel.time(), Sr87_86_CorrRb)

        StdSplineRb_Sr87_86 = data.spline(Sr_rmName, "Sr87_86_CorrRb").data()
        StdCorrRb_Sr87_86 = Sr87_86_CorrRb * StdValue_Sr87_86 / StdSplineRb_Sr87_86

    totalSrBeam = Sr84_corr + Sr86_corr + Sr87_CorrRb + Sr88_corr
    data.createTimeSeries('totalSrBeam', data.Intermediate, indexChannel.time(), totalSrBeam)

    # Calculating fractionation-corrected Rb/Sr (if "Rb_Sr_elemental" is selcted)
    if Rb_Sr_elemental:
        data.createTimeSeries('Rb87_Sr86_CorrRb', data.Intermediate, indexChannel.time(), Rb87_Sr86_CorrRb)

        try:
            Sr_conc = data.referenceMaterialData(Ref_Rb_Sr_elemental)["Sr"].value()
            Rb_conc = data.referenceMaterialData(Ref_Rb_Sr_elemental)["Rb"].value()
            StdSpline_RbSr = data.spline(Ref_Rb_Sr_elemental, "Rb87_Sr86_CorrRb").data()
        except:
            IoLog.error("The Combined Sr DRS requires Rb/Sr Ref Material selections to proceed.")
            drs.message("DRS did not finish. Please check Messages")
            drs.progress(100)
            drs.finished()
            return

        Rb87_Sr86_final = Rb87_Sr86_CorrRb * ((Rb_conc/Sr_conc) * (27.83 / 9.86)) / StdSpline_RbSr

    else:
        Rb87_Sr86_final = Rb87_Sr86_CorrRb * Rb_Sr_fract

    # Calculating age-corrected 87Sr/86Sr
    Sr87_86_AgeCorr = StdCorrRb_Sr87_86 - Rb87_Sr86_final * (math.exp(0.000013972 * Age)-1)# 87Rb decay constant from Villa et al. (2015, Geoch. et Cosm. Acta)

    # Output channels
    try:
        output_channels_names = ['Er_Sr_ppm','Yb_Sr_ppm','Sr84_86_Corr','Sr84_88_Corr','Rb87_Sr86_final','StdCorr_Sr87_86','StdCorrRb_Sr87_86','Sr87_86_AgeCorr']
        output_channels = [Er_Sr_ppm,Yb_Sr_ppm,Sr84_86_Corr,Sr84_88_Corr,Rb87_Sr86_final,StdCorr_Sr87_86,StdCorrRb_Sr87_86,Sr87_86_AgeCorr]
        for name, channel in zip(output_channels_names, output_channels):
            data.createTimeSeries(name, data.Output, indexChannel.time(), channel)

    except:

        output_channels_names = ['Sr84_86_Corr','Sr84_88_Corr','Rb87_Sr86_final','StdCorr_Sr87_86','StdCorrRb_Sr87_86','Sr87_86_AgeCorr']
        output_channels = [Sr84_86_Corr,Sr84_88_Corr,Rb87_Sr86_final,StdCorr_Sr87_86,StdCorrRb_Sr87_86,Sr87_86_AgeCorr]
        for name, channel in zip(output_channels_names, output_channels):
            data.createTimeSeries(name, data.Output, indexChannel.time(), channel)

    if propErrors:
        drs.message("Propagating errors...")
        drs.progress(90)

        groups = [s for s in data.selectionGroupList() if s.type != data.Baseline]
        data.propagateErrors(groups, [data.timeSeries("StdCorrRb_Sr87_86")], data.timeSeries("StdCorr_Sr87_86"), Sr_rmName)

##############################

    # CaPO Plot and Correction
    if len(settings['CaPO_RMs']) > 0:
        print("Here are the RMs for CaPO correction:", settings['CaPO_RMs'])
        PLOT.clearGraphs()

        # Get Sr signal and deviations for chosen RMs:
        Sr8786 = data.timeSeries('StdCorrRb_Sr87_86')
        totalSr = data.timeSeries('totalSrBeam')

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
                print(f"Could not get the selections for group {rm} so it could not be used in the CaPO correction: {err}")
                IoLog.warning(f"Could not get the selections for group {rm} so it could not be used in the CaPO correction: {err}")
                continue

            try:
                true = data.referenceMaterialData(rm)['87Sr/86Sr'].value()
            except KeyError as err:
                print(f"Could not get the \'87Sr/86Sr\' value for {rm} so it could not be used in the CaPO correction: {err}")
                IoLog.warning(f"Could not get the \'87Sr/86Sr\' value for {rm} so it could not be used in the CaPO correction: {err}")
                continue

            print(rm)
            print(true)

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

        sigs_array = np.array(sigs_main)
        devs_array = np.array(devs_main)
        #
        # Show the fit +/- 10% of the spread in signals
        sigs_range = sigs_array.max() - sigs_array.min()
        fit_x = np.linspace(
            sigs_array.min() - sigs_range * 0.1,
            sigs_array.max() + sigs_range * 0.1,
            50
        )

        # Now add fit to all data
        if settings["CaPOEqType"] == 'Linear':
            def fitFunc(x, a, b):
                return a + b*x

            params, cov = curve_fit(fitFunc, sigs_array, devs_array, ftol=1e-5)

            print(f"Here are the fit params: Slope: {params[1]}, Intercept: {params[0]}")
            fit_y = params[1]*fit_x + params[0]

        '''
           If user has chosen exp decay model, it might not find a fit.
           If so, we'll revert to linear fit
        '''
        if settings["CaPOEqType"] == 'Exponential Decay':
            def f(x, a, b, c):
                return a + b * np.exp(-c * x)

            try:
                params, cov = curve_fit(f, sigs_array, devs_array)
                print(f"Here are the fit params: {params}")
                fit_y = params[0] + params[1] * np.exp(-params[2] * fit_x)

            except RuntimeError as e:
                if 'Optimal parameters not found: Number of calls' in str(e):

                    IoLog.warning('Could not find fit to data with Exp model. Switching to Linear model.')

                    settings['CaPOEqType'] = 'Linear'

                    def fitFunc(x, a, b):
                        return a + b*x

                    params, cov = curve_fit(fitFunc, sigs_array, devs_array, ftol=1e-5)

                    print(f"Here are the fit params: Slope: {params[1]}, Intercept: {params[0]}")
                    fit_y = params[1]*fit_x + params[0]
                else:
                    raise

        g = PLOT.addGraph()
        g.setName("fit")
        fit_pen = QPen(QColor('dark grey'))
        fit_pen.setStyle(Qt.DashLine);
        g.pen = fit_pen
        g.setData(fit_x, fit_y)

        ann.visible = True
        ann.text = f'''
            <p style="color:black;">
            <b>Fit Parameters:</b><br>
            Slope: {params[1]:.7f} <br>
            Intercept: {params[0]:.6f}</p>'''

        PLOT.left().label = '87Sr/86Sr Deviation (meas/true)'
        PLOT.bottom().label = 'Total Sr (V)'
        PLOT.setToolsVisible(False)
        PLOT.rescaleAxes()
        PLOT.replot()

        if settings['CaPOEqType'] == 'Linear':
            CaPO_corrAmt = totalSrBeam * params[1] + params[0]
        else:
            CaPO_corrAmt = params[0] + params[1] * np.exp(-params[2] * totalSrBeam)

        data.createTimeSeries('CaPO_corrAmt', data.Output, indexChannel.time(), CaPO_corrAmt)

        CaPOCorr_Sr8786 = StdCorrRb_Sr87_86 / CaPO_corrAmt
        data.createTimeSeries('CaPOCorr_Sr8786', data.Output, indexChannel.time(), CaPOCorr_Sr8786)

    else:
        print("No CaPO correction invoked (no selected RMs)")
        PLOT.clearGraphs()
        ann.visible = False
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

    drs.setSetting("IndexChannel", defaultChannelName)
    drs.setSetting("ReferenceMaterial", "CO3_shell")
    drs.setSetting("Mask", True)
    drs.setSetting("MaskChannel", defaultChannelName)
    drs.setSetting("MaskCutoff", 0.05)
    drs.setSetting("MaskTrim", 0.0)
    drs.setSetting("Rb_Beta_adjust", False)
    drs.setSetting("RbBias", 1.)
    drs.setSetting("ReferenceMaterialRb", "G_BCR2G")
    drs.setSetting("Rb_Sr_elemental", False)
    drs.setSetting("Rb_Sr_fract", 1.)
    drs.setSetting("ReferenceMaterialRbSr", "G_BCR2G")
    drs.setSetting("Age", 0.)
    drs.setSetting("Dy_Er", 1.5)
    drs.setSetting("Lu_Yb", 0.15)
    drs.setSetting("ProportionCaAr", 1.)
    drs.setSetting("ProportionNaNi", 1.)
    drs.setSetting("Sr88_86_reference", 8.37520938)  #Konter & Storm (2014)
    drs.setSetting("Rb87_85_reference", 0.385710)    #Konter & Storm (2014)
    drs.setSetting("REE_subtract", False)
    drs.setSetting("REEBias", 1.)
    drs.setSetting("CaAr_CaCa_subtract", False)
    drs.setSetting("CaAr_83", False)
    drs.setSetting("CaArBias", 1.)
    drs.setSetting("NaNi_CaAlO_subtract", False)
    drs.setSetting("NaNiBias", 1.)
    drs.setSetting("CaPOEqType", 'Linear')
    drs.setSetting("CaPO_RMs", [])
    drs.setSetting("PropagateError", False)

    settings = drs.settings()

    DRS_message = QtGui.QLabel("This DRS requires data at least 86.5 and 83.5 half-masses for the REE subtraction (also 87.5 and 81.5 for direct 176Lu++ and 164Dy++ correction), mass 82 or 83 for the CaAr-CaCa subtraction, and mass 83 for the NaNi-CaAlO subtraction.\nIf you did not collect data for these channels, do not check the corresponding boxes. Otherwise you will get an error message.")
    DRS_message.setStyleSheet('color:yellow')
    formLayout.addRow(DRS_message)

    verticalSpacer = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
    formLayout.addItem(verticalSpacer)

    indexComboBox = QtGui.QComboBox(widget)
    indexComboBox.addItems(timeSeriesNames)
    indexComboBox.setCurrentText(settings["IndexChannel"])
    indexComboBox.textActivated.connect(lambda t: drs.setSetting("IndexChannel", t))
    formLayout.addRow("Index channel", indexComboBox)

    def updateIndexCombo():
        timeSeriesNames = data.timeSeriesNames(data.Input)
        indexComboBox.clear()
        indexComboBox.addItems(timeSeriesNames)

    data.dataChanged.connect(updateIndexCombo)

    rmComboBox = QtGui.QComboBox(widget)
    rmNames = data.selectionGroupNames(data.ReferenceMaterial)
    rmComboBox.addItems(rmNames)
    if settings["ReferenceMaterial"] in rmNames:
        rmComboBox.setCurrentText(settings["ReferenceMaterial"])
        drs.setSetting("ReferenceMaterial", settings["ReferenceMaterial"])

    rmComboBox.textActivated.connect(lambda t: drs.setSetting("ReferenceMaterial", t))
    formLayout.addRow("Reference material", rmComboBox)

    def updateRMCombo():
        rmNames = data.selectionGroupNames(data.ReferenceMaterial)
        rmComboBox.clear()
        rmComboBox.addItems(rmNames)

    data.selectionGroupsChanged.connect(updateRMCombo)

    verticalSpacer2 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
    formLayout.addItem(verticalSpacer2)

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

    verticalSpacer3 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
    formLayout.addItem(verticalSpacer3)

    REECheckBox = QtGui.QCheckBox(widget)
    REECheckBox.setChecked(settings["REE_subtract"])
    REECheckBox.toggled.connect(lambda t: drs.setSetting("REE_subtract", bool(t)))
    formLayout.addRow("REE subtraction?", REECheckBox)

    dyErLineEdit = QtGui.QLineEdit(widget)
    dyErLineEdit.setText(settings["Dy_Er"])
    dyErLineEdit.textChanged.connect(lambda t: drs.setSetting("Dy_Er", float(t)))
    formLayout.addRow("Dy/Er ratio (default is approximately chondritic)", dyErLineEdit)

    LuYbLineEdit = QtGui.QLineEdit(widget)
    LuYbLineEdit.setText(settings["Lu_Yb"])
    LuYbLineEdit.textChanged.connect(lambda t: drs.setSetting("Lu_Yb", float(t)))
    formLayout.addRow("Lu/Yb ratio (default is approximately chondritic)", LuYbLineEdit)

    REEBiasLineEdit = QtGui.QLineEdit(widget)
    REEBiasLineEdit.setText(settings["REEBias"])
    REEBiasLineEdit.textChanged.connect(lambda t: drs.setSetting("REEBias", float(t)))
    formLayout.addRow("Scale REE Beta (1 = BetaSr)", REEBiasLineEdit)

    verticalSpacer4 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
    formLayout.addItem(verticalSpacer4)

    CaArCheckBox = QtGui.QCheckBox(widget)
    CaArCheckBox.setChecked(settings["CaAr_CaCa_subtract"])
    CaArCheckBox.toggled.connect(lambda t: drs.setSetting("CaAr_CaCa_subtract", bool(t)))
    formLayout.addRow("CaAr-CaCa subtraction?", CaArCheckBox)

    CaAr83CheckBox = QtGui.QCheckBox(widget)
    CaAr83CheckBox.setChecked(settings["CaAr_83"])
    CaAr83CheckBox.toggled.connect(lambda t: drs.setSetting("CaAr_83", bool(t)))
    formLayout.addRow("Use 83(CaAr,CaCa) as monitor? (uses 82 as default otherwise)", CaAr83CheckBox)

    CaArLineEdit = QtGui.QLineEdit(widget)
    CaArLineEdit.setText(settings["ProportionCaAr"])
    CaArLineEdit.textChanged.connect(lambda t: drs.setSetting("ProportionCaAr", float(t)))
    formLayout.addRow("Proportion CaAr (0 to 1)", CaArLineEdit)

    CaArBiasLineEdit = QtGui.QLineEdit(widget)
    CaArBiasLineEdit.setText(settings["CaArBias"])
    CaArBiasLineEdit.textChanged.connect(lambda t: drs.setSetting("CaArBias", float(t)))
    formLayout.addRow("Scale CaAr-CaCa Beta (1 = BetaSr)", CaArBiasLineEdit)

    verticalSpacer5 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
    formLayout.addItem(verticalSpacer5)

    NaNiCheckBox = QtGui.QCheckBox(widget)
    NaNiCheckBox.setChecked(settings["NaNi_CaAlO_subtract"])
    NaNiCheckBox.toggled.connect(lambda t: drs.setSetting("NaNi_CaAlO_subtract", bool(t)))
    formLayout.addRow("NaNi-CaAlO subtraction?", NaNiCheckBox)

    NaNiLineEdit = QtGui.QLineEdit(widget)
    NaNiLineEdit.setText(settings["ProportionNaNi"])
    NaNiLineEdit.textChanged.connect(lambda t: drs.setSetting("ProportionNaNi", float(t)))
    formLayout.addRow("Proportion NaNi (0 to 1)", NaNiLineEdit)

    NaNiBiasLineEdit = QtGui.QLineEdit(widget)
    NaNiBiasLineEdit.setText(settings["NaNiBias"])
    NaNiBiasLineEdit.textChanged.connect(lambda t: drs.setSetting("NaNiBias", float(t)))
    formLayout.addRow("Scale NaNi-CaAlO Beta (1 = BetaSr)", NaNiBiasLineEdit)

    verticalSpacer6 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
    formLayout.addItem(verticalSpacer6)

    RbBetaCheckBox = QtGui.QCheckBox(widget)
    RbBetaCheckBox.setChecked(settings["Rb_Beta_adjust"])
    RbBetaCheckBox.toggled.connect(lambda t: drs.setSetting("Rb_Beta_adjust", bool(t)))
    formLayout.addRow("Rb Beta adjust?", RbBetaCheckBox)

    RbrmComboBox = QtGui.QComboBox(widget)
    RbrmComboBox.addItems(rmNames)
    RbrmComboBox.setCurrentText(settings["ReferenceMaterialRb"])
    RbrmComboBox.currentTextChanged.connect(lambda t: drs.setSetting("ReferenceMaterialRb", t))
    formLayout.addRow("Reference material Rb beta", RbrmComboBox)

    RbBiasLineEdit = QtGui.QLineEdit(widget)
    RbBiasLineEdit.setText(settings["RbBias"])
    RbBiasLineEdit.textChanged.connect(lambda t: drs.setSetting("RbBias", float(t)))
    formLayout.addRow("Scale Rb Beta (if Rb Beta adjust unchecked; 1 = BetaSr)", RbBiasLineEdit)

    verticalSpacer7 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
    formLayout.addItem(verticalSpacer7)

    RbSrCheckBox = QtGui.QCheckBox(widget)
    RbSrCheckBox.setChecked(settings["Rb_Sr_elemental"])
    RbSrCheckBox.toggled.connect(lambda t: drs.setSetting("Rb_Sr_elemental", bool(t)))
    formLayout.addRow("Rb-Sr fractionation?", RbSrCheckBox)

    RbSrrmComboBox = QtGui.QComboBox(widget)
    RbSrrmComboBox.addItems(rmNames)
    RbSrrmComboBox.setCurrentText(settings["ReferenceMaterialRbSr"])
    RbSrrmComboBox.currentTextChanged.connect(lambda t: drs.setSetting("ReferenceMaterialRbSr", t))
    formLayout.addRow("Reference material Rb/Sr fractionation", RbSrrmComboBox)

    RbSrLineEdit = QtGui.QLineEdit(widget)
    RbSrLineEdit.setText(settings["Rb_Sr_fract"])
    RbSrLineEdit.textChanged.connect(lambda t: drs.setSetting("Rb_Sr_fract", float(t)))
    formLayout.addRow("Scale Rb/Sr (if no Rb/Sr standard)", RbSrLineEdit)


    verticalSpacer8 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
    formLayout.addItem(verticalSpacer8)


    ageLineEdit = QtGui.QLineEdit(widget)
    ageLineEdit.setText(settings["Age"])
    ageLineEdit.textChanged.connect(lambda t: drs.setSetting("Age", float(t)))
    formLayout.addRow("Age (Ma)", ageLineEdit)


    verticalSpacer9 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
    formLayout.addItem(verticalSpacer9)

    # CaPO controls
    capoEqType = QtGui.QComboBox(widget)
    capoEqType.addItems(['Linear', 'Exponential Decay'])
    capoEqType.setCurrentText(settings["CaPOEqType"])
    capoEqType.currentTextChanged.connect(lambda t: drs.setSetting("CaPOEqType", t))
    formLayout.addRow("CaPO Fit Eqn", capoEqType)

    setExtButton = QtGui.QToolButton(widget)
    rmMenu = ReferenceMaterialsMenu(setExtButton)
    rmMenu.rmsChanged.connect(lambda l: drs.setSetting("CaPO_RMs", l))
    setExtButton.setMenu(rmMenu)
    setExtButton.setPopupMode(QtGui.QToolButton.InstantPopup)
    setExtButton.setIcon(CUI().icon('trophy'))
    setExtButton.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
    formLayout.addRow("RMs for CaPO correction:", setExtButton)

    # Plot for CaPO correction
    formLayout.addRow('CaPO Correction fit', PLOT)

    # Prop errors controls
    propCheckBox = QtGui.QCheckBox(widget)
    propCheckBox.setChecked(settings["PropagateError"])
    propCheckBox.toggled.connect(lambda t: drs.setSetting("PropagateError", bool(t)))
    formLayout.addRow("Propagate Errors?", propCheckBox)


    # Restore settings
    try:
        settings = drs.settings()
        print('Restoring settings...')
        print(settings)
        rmComboBox.setCurrentText(settings["ReferenceMaterial"])
    except KeyError:
        pass

    drs.setSettingsWidget(widget)
