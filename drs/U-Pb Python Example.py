#/ Type: DRS
#/ Name: U-Pb Python Example
#/ Authors: Joe Petrus and Bence Paul
#/ Description: Simple U-Pb with downhole fractionation corrections
#/ References: Paton et al., 2010 G3
#/ Version: 1.0
#/ Contact: support@iolite-software.com

from iolite import QtGui
from iolite.Qt import Qt
import numpy as np
from scipy.optimize import curve_fit

def downholeFunc(x, a, b, c):
    return a + b * np.exp(-c * x)

def runDRS():
    drs.message("Starting baseline subtract DRS...")
    drs.progress(0)

    # Get settings
    settings = drs.settings()
    print(settings)

    indexChannel = data.timeSeries('U238')

    # Setup index time
    drs.message("Setting up index time...")
    drs.progress(5)
    drs.setIndexChannel(indexChannel)

    # Interp onto index time and baseline subtract
    drs.message("Interpolating onto index time and baseline subtracting...")
    drs.progress(25)

    allInputChannels = data.timeSeriesList(data.Input)
    
    commonProps = {'DRS': drs.name()}

    mask = drs.createMaskFromLaserLog(0)
    drs.baselineSubtract(data.selectionGroup('Baseline'), data.timeSeriesList(data.Input), mask, 25, 50)
    drs.createBeamSecondsFromLaserLog()
    
    U238 = data.timeSeries('U238').data()
    Pb206 = data.timeSeries('Pb206').data()
    beamSeconds = data.timeSeries('BeamSeconds').data()
    
    # Create raw ratio
    Pb206_U238 = Pb206/U238
    Pb206_U238_ts = data.createTimeSeries('Pb206/U238', data.Intermediate, indexChannel.time(), Pb206_U238)
    
    # DHF correct it
    DHF = data.compileDownhole(data.selectionGroup(settings['ReferenceMaterial']), Pb206_U238_ts)
    DHFnans = np.isnan(DHF[1])
    DHFt = DHF[0][~DHFnans]
    DHFx = DHF[1][~DHFnans]
    popt, pcov = curve_fit(downholeFunc, DHFt, DHFx)
    Pb206_U238_dc = Pb206_U238/( 1 + (popt[1]/popt[0])*np.exp(-popt[2]*beamSeconds) )
    data.createTimeSeries('DC Pb206/U238', data.Intermediate, indexChannel.time(), Pb206_U238_dc)

    # Calibrate it
    rm = data.referenceMaterialData(settings["ReferenceMaterial"])
    rm6_38 = rm['206Pb/238U'].value()
    rm6_38_spline = data.spline(settings['ReferenceMaterial'], 'DC Pb206/U238').data()
    
    Pb206_U238_final = (rm6_38/rm6_38_spline)*Pb206_U238_dc
    data.createTimeSeries('Final Pb206/U238', data.Output, indexChannel.time(), Pb206_U238_final)
   

    drs.message("Finished!")
    drs.progress(100)
    drs.finished()
    

def settingsWidget():
    widget = QtGui.QWidget()
    formLayout = QtGui.QFormLayout()
    formLayout.setFieldGrowthPolicy(QtGui.QFormLayout.FieldsStayAtSizeHint)
    formLayout.setFormAlignment(Qt.AlignHCenter | Qt.AlignTop)
    widget.setLayout(formLayout)

    rmComboBox = QtGui.QComboBox(widget)
    rmComboBox.addItems(data.referenceMaterialNames())
    rmComboBox.currentTextChanged.connect(lambda s: drs.setSetting("ReferenceMaterial", str(s)))
    formLayout.addRow("Reference material", rmComboBox)
    
    # Restore settings
    try:
        settings = drs.settings()
        rmComboBox.setCurrentText(settings["External"])
    except KeyError:
        pass

    drs.setSettingsWidget(widget)
