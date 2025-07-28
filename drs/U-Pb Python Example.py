#/ Type: DRS
#/ Name: U-Pb Python Example
#/ Authors: Joe Petrus and Bence Paul
#/ Description: Simple U-Pb with downhole fractionation corrections
#/ References: Paton et al., 2010 G3
#/ Version: 1.0
#/ Contact: support@iolite-software.com

from iolite import QtGui
from iolite.Qt import Qt, QColor
from iolite.ui import IolitePlotPyInterface as Plot
import numpy as np
from scipy.optimize import curve_fit, leastsq

# Constants
l238 = 1.55125e-10
l235 = 9.8485e-10
l232 = 0.49475e-10
k = 137.818

# Lookup tables for 7/6 age 
lut = np.logspace(0, 23, 1000, base=2.7)
lu76 = (1/k)*(np.exp(l235*lut) - 1)/(np.exp(l238*lut) - 1)

def downholeFunc_Lin(t, a, b):
    return a*t + b

def downholeFunc_Exp(t, a, b, c):
    return a + b*np.exp(-c * t)

def downholeFunc_Exp_Lin(t, a, b, c, d):
    return a + b*t + c*np.exp(-d * t)


def age638(r):
    return np.log(r + 1)/l238


def age735(r):
    return np.log(r + 1)/l235


def age832(r):
    return np.log(r + 1)/l232
    

def age76(r):
    return np.interp(r, lu76, lut)
    

def runDRS():
    drs.message("Starting baseline subtract DRS...")
    drs.progress(0)

    # Get settings
    settings = drs.settings()
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
    beamSeconds = data.timeSeries('BeamSeconds').data()

    timeStep = indexChannel.time()[1] - indexChannel.time()[0]
    startTrimSec = settings["StartTrim"]
    startTrimIndex = int(round(startTrimSec/timeStep))
    endTrimSec = settings["EndTrim"]
    endTrimIndex = int(round(endTrimSec/timeStep))

    ratios = [
        {
            'data': lambda: data.timeSeries('Pb206_CPS').data()/data.timeSeries('U238_CPS').data(),
            'dhfc': True,
            'age': age638,
            'name': 'Pb206/U238',
            'rmName': '206Pb/238U',
            'plotName': 'Plot68'  
        },
        {
            'data': lambda: k*data.timeSeries('Pb207_CPS').data()/data.timeSeries('U238_CPS').data(),
            'dhfc': True,
            'age': age735,
            'name': 'Pb207/U235',
            'rmName': '207Pb/235U',
            'plotName': 'Plot75'
        },
        {
            'data': lambda: data.timeSeries('Pb208_CPS').data()/data.timeSeries('Th232_CPS').data(),
            'dhfc': True,
            'age': age832,
            'name': 'Pb208/Th232',
            'rmName': '208Pb/232Th',
            'plotName': 'Plot82'
        },
        {
            'data': lambda: data.timeSeries('Pb207_CPS').data()/data.timeSeries('Pb206_CPS').data(),
            'dhfc': False,
            'age': age76,
            'name': 'Pb207/Pb206',
            'rmName': '207Pb/206Pb',
        }
    ]

    # Clear previous plots
    settings['FitsWidget'].clear()

    def processRatio(ratio):
        print('Processing ratio: ' + ratio['name'])
        rawRatio = ratio['data']()
        rawRatio[np.isinf(rawRatio)] = np.nan
        ts = data.createTimeSeries(ratio['name'], data.Intermediate, indexChannel.time(), rawRatio, commonProps)

        ratioToCalibrate = rawRatio
        ratioToCalibrateName = ratio['name']

        if ratio['dhfc']:
            print(f'... doing down-hole correction for {ratio["name"]}')
            DHF = data.compileDownhole(data.selectionGroup(settings['ReferenceMaterial']), ts)
            DHFnans = np.isnan(DHF[1])
            DHFt = DHF[0][~DHFnans]
            DHFr = DHF[1][~DHFnans]

            if startTrimIndex != 0:
                DHFt = DHFt[startTrimIndex:]
                DHFr = DHFr[startTrimIndex:]
    
            if endTrimIndex != 0:
                DHFt = DHFt[:-endTrimIndex]
                DHFr = DHFr[:-endTrimIndex]
            
            print('... fitting down-hole function')
            print(f'... fit model: {settings["Model"]}')
            try:
                if settings['Model'] == 'Linear':
                    params, cov = curve_fit(downholeFunc_Lin, DHFt, DHFr, ftol=1e-5, maxfev = 3000)
                elif settings['Model'] == 'Exponential':
                    params, cov = curve_fit(downholeFunc_Exp, DHFt, DHFr, ftol=1e-5, maxfev = 3000)
                elif settings['Model'] == 'Exponential + Linear':
                    params, cov = curve_fit(downholeFunc_Exp_Lin, DHFt, DHFr, ftol=1e-5, maxfev = 3000)
            except RuntimeError as err:
                raise RuntimeError(f'Could not fit down-hole function: {err}')

            #print(f'... fit parameters: {params}')
            #print(f'... covariance: {cov}')
            param_name_list = ['a', 'b', 'c', 'd']
            uncertainties_list = []
            for i, param in enumerate(params):
                sigma = np.sqrt(cov[i, i])
                uncertainties_list.append(f'{param_name_list[i]}: {param} +/- {sigma:.4e} ({sigma/param:.4%})')
                print(f'... {param_name_list[i]}: {param} +/- {sigma:.4e} ({sigma/param:.4%})')

            if settings['Model'] == 'Linear':
                dc = rawRatio / (1 + (params[0] / params[1]) * beamSeconds)
            elif settings['Model'] == 'Exponential':
                dc = rawRatio / (1 + (params[1] / params[0]) * np.exp(-params[2] * beamSeconds))
            elif settings['Model'] == 'Exponential + Linear':
                dc = rawRatio / (1 + (params[1] / params[0]) * beamSeconds + (params[2] / params[0]) * np.exp(-params[3] * beamSeconds))

            data.createTimeSeries('DC '+ratio['name'], data.Intermediate, indexChannel.time(), dc, commonProps)

            plot = Plot(settings['FitsWidget'])
            g = plot.addGraph()
            g.setData(DHFt, DHFr)
            g2 = plot.addGraph()
            g2.setColor(QColor(255, 0, 0))

            if settings['Model'] == 'Linear':
                g2.setData(DHFt, downholeFunc_Lin(DHFt, params[0], params[1]))
            elif settings['Model'] == 'Exponential':
                g2.setData(DHFt, downholeFunc_Exp(DHFt, params[0], params[1], params[2]))
            elif settings['Model'] == 'Exponential + Linear':
                g2.setData(DHFt, downholeFunc_Exp_Lin(DHFt, params[0], params[1], params[2], params[3]))
            
            # TODO: Add uncertainties to plot as text
            
            plot.left().label = ratio['name']
            plot.bottom().label = 'Time (s)'
            plot.setToolsVisible(False)
            plot.rescaleAxes()
            plot.replot()
            settings['FitsWidget'].addTab(plot, ratio['name'])

            ratioToCalibrate = dc
            ratioToCalibrateName = 'DC ' + ratio['name']

        rm = data.referenceMaterialData(settings["ReferenceMaterial"])
        rmValue = rm[ratio['rmName']].value()
        rmSpline = data.spline(settings['ReferenceMaterial'], ratioToCalibrateName).data()
        finalRatio = (rmValue/rmSpline)*ratioToCalibrate
        data.createTimeSeries('Final ' + ratio['name'], data.Output, indexChannel.time(), finalRatio, commonProps)

        finalAge = ratio['age'](finalRatio)/1e6
        data.createTimeSeries('Final ' + ratio['name'] + ' age', data.Output, indexChannel.time(), finalAge, commonProps)
   
    for i, ratio in enumerate(ratios):
        drs.message('Working on ' + ratio['name'])
        drs.progress(50 + 50*float(i)/float(len(ratios)))
        try:
            processRatio(ratio)
        except RuntimeError as err:
            IoLog.warning('Could not process ratio %s: %s'%(ratio['name'], err)) 

    drs.message("Finished!")
    drs.progress(100)
    drs.finished()
    

def settingsWidget():
    widget = QtGui.QWidget()
    formLayout = QtGui.QFormLayout()
    formLayout.setFieldGrowthPolicy(QtGui.QFormLayout.AllNonFixedFieldsGrow )
    formLayout.setFormAlignment(Qt.AlignHCenter | Qt.AlignTop)
    widget.setLayout(formLayout)

    # Some default settings
    drs.setSetting("ReferenceMaterial", "Z_91500")
    drs.setSetting("Model", "Linear")
    drs.setSetting("StartTrim", 0)
    drs.setSetting("EndTrim", 0)

    rmComboBox = QtGui.QComboBox(widget)
    rmComboBox.setFixedWidth(150)
    rmComboBox.addItems(data.referenceMaterialNames())
    rmComboBox.setCurrentText(drs.setting("ReferenceMaterial"))
    rmComboBox.currentTextChanged.connect(lambda s: drs.setSetting("ReferenceMaterial", str(s)))
    formLayout.addRow("Reference material", rmComboBox)

    modelComboBox = QtGui.QComboBox(widget)
    modelComboBox.setFixedWidth(150)
    modelComboBox.addItems([ 'Linear', 'Exponential', 'Exponential + Linear'])
    modelComboBox.setCurrentText(drs.setting("Model"))
    modelComboBox.currentTextChanged.connect(lambda s: drs.setSetting("Model", str(s)))
    formLayout.addRow("Fit Model", modelComboBox)

    startLineEdit = QtGui.QLineEdit(widget)
    startLineEdit.setFixedWidth(150)
    startLineEdit.setText(str(drs.setting("StartTrim")))
    startLineEdit.textEdited.connect(lambda s: drs.setSetting("StartTrim", float(s)))
    formLayout.addRow("Start trim (s)", startLineEdit)

    endLineEdit = QtGui.QLineEdit(widget)
    endLineEdit.setFixedWidth(150)
    endLineEdit.setText(str(drs.setting("EndTrim")))
    endLineEdit.textEdited.connect(lambda s: drs.setSetting("EndTrim", float(s)))
    formLayout.addRow("End trim (s)", endLineEdit)
    
    tabWidget = QtGui.QTabWidget(widget)
    tabWidget.setFixedSize(500,400)
    formLayout.addRow('Down-hole fits', tabWidget)
    drs.setSetting('FitsWidget', tabWidget)

    # Restore settings
    try:
        settings = drs.settings()
        print('Restoring settings...')
        print(settings)
        rmComboBox.setCurrentText(settings["ReferenceMaterial"])
        startLineEdit.setText(str(settings['StartTrim']))
        endLineEdit.setText(str(settings['EndTrim']))
    except KeyError:
        pass

    drs.setSettingsWidget(widget)
