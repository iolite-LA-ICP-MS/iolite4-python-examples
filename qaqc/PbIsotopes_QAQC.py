#/ Type: QAQC
#/ Name: Pb Isotopes QAQC
#/ Authors: Bence Paul and Joe Petrus
#/ Description: This module compares measured Pb isotopes with accepted values
#/ References: None
#/ Version: 1.0
#/ Contact: support@iolite-software.com

from iolite.QtGui import QWidget, QFormLayout, QComboBox, QDoubleSpinBox
from iolite.QtCore import Qt
from iolite.ui import IolitePlotPyInterface as Plot
from iolite.ui import QCPErrorBars, QCPRange

import numpy as np


plot = Plot()

def update():
    print('Running Pb Isotopes QAQC update...')
    qaqc.clearReport()
    plot.clearGraphs()

    if "Final 208Pb/206Pb" not in data.timeSeriesNames():
        qaqc.pushHtml('This module must be used after running the Pb Isotopes DRS. Aborting.')
        qaqc.finished(qaqc.Error)
        return

    rmName = qaqc.settings()["RefMat"]
    sg = data.selectionGroup(qaqc.settings()["RefMat"])
    if sg is None:
        qaqc.pushHtml('Reference material not found. Aborting.')
        qaqc.finished(qaqc.Error)
        return

    # These lists will be used to create a plot
    x = []
    y = []
    yerr = []
    std_y = []
    std_yerr = []
    ticks = [] # A list of strings for the x-axis labels

    qaqc.pushHtml(f'<h3>Reference material: {rmName}</h3>')
    qaqc.pushHtml(f'<h3>Allowable difference: {float(qaqc.settings()["AllowableDiff_pct"]):.2f} %</h3>')
    qaqc.pushHtml('')
    qaqc.pushHtml('<table cellpadding="10" border="1">')
    qaqc.pushHtml('<tr><th>Ratio</th><th>Measured</th><th>Accepted</th><th>% Difference</th></tr>')

    if "206Pb/204Pb" in data.referenceMaterialData(rmName).keys():
        std_64 = data.referenceMaterialData(rmName)["206Pb/204Pb"].value()
        std_64_uncert = data.referenceMaterialData(rmName)["206Pb/204Pb"].uncertainty()
        tsd_64 = data.timeSeries("Final 206Pb/204Pb")
        measured_64 = data.groupResult(sg, tsd_64).value()
        measured_64_uncert = data.groupResult(sg, tsd_64).uncertainty()
        qaqc.pushHtml(f'<tr><td>206Pb/204Pb</td><td>{measured_64:.3f}±{measured_64_uncert:.3f}</td><td>{std_64}±{std_64_uncert}</td><td>{100*(measured_64-std_64)/std_64:.2f}%</td></tr>')
        x.append(1)
        y.append(measured_64)
        yerr.append(measured_64_uncert)
        std_y.append(std_64)
        std_yerr.append(std_64_uncert)
        ticks.append('206Pb/204Pb')

    if "207Pb/204Pb" in data.referenceMaterialData(rmName).keys():
        std_74 = data.referenceMaterialData(rmName)["207Pb/204Pb"].value()
        std_74_uncert = data.referenceMaterialData(rmName)["207Pb/204Pb"].uncertainty()
        tsd_74 = data.timeSeries("Final 207Pb/204Pb")
        measured_74 = data.groupResult(sg, tsd_74).value()
        measured_74_uncert = data.groupResult(sg, tsd_74).uncertainty()
        qaqc.pushHtml(f'<tr><td>207Pb/204Pb</td><td>{measured_74:.3f}±{measured_74_uncert:.3f}</td><td>{std_74}±{std_74_uncert}</td><td>{100*(measured_74-std_74)/std_74:.2f}%</td></tr>')
        x.append(2)
        y.append(measured_74)
        yerr.append(measured_74_uncert)
        std_y.append(std_74)
        std_yerr.append(std_74_uncert)
        ticks.append('207Pb/204Pb')

    if "208Pb/204Pb" in data.referenceMaterialData(rmName).keys():
        std_84 = data.referenceMaterialData(rmName)["208Pb/204Pb"].value()
        std_84_uncert = data.referenceMaterialData(rmName)["208Pb/204Pb"].uncertainty()
        tsd_84 = data.timeSeries("Final 208Pb/204Pb")
        measured_84 = data.groupResult(sg, tsd_84).value()
        measured_84_uncert = data.groupResult(sg, tsd_84).uncertainty()
        qaqc.pushHtml(f'<tr><td>208Pb/204Pb</td><td>{measured_84:.3f}±{measured_84_uncert:.3f}</td><td>{std_84}±{std_84_uncert}</td><td>{100*(measured_84-std_84)/std_84:.2f}%</td></tr>')
        x.append(3)
        y.append(measured_84)
        yerr.append(measured_84_uncert)
        std_y.append(std_84)
        std_yerr.append(std_84_uncert)
        ticks.append('208Pb/204Pb')

    if "208Pb/206Pb" in data.referenceMaterialData(rmName).keys():
        std_86 = data.referenceMaterialData(rmName)["208Pb/206Pb"].value()
        std_86_uncert = data.referenceMaterialData(rmName)["208Pb/206Pb"].uncertainty()
        tsd_86 = data.timeSeries("Final 208Pb/206Pb")
        measured_86 = data.groupResult(sg, tsd_86).value()
        measured_86_uncert = data.groupResult(sg, tsd_86).uncertainty()
        qaqc.pushHtml(f'<tr><td>208Pb/206Pb</td><td>{measured_86:.3f}±{measured_86_uncert:.3f}</td><td>{std_86}±{std_86_uncert}</td><td>{100*(measured_86-std_86)/std_86:.2f}%</td></tr>')
        x.append(4)
        y.append(measured_86)
        yerr.append(measured_86_uncert)
        std_y.append(std_86)
        std_yerr.append(std_86_uncert)
        ticks.append('208Pb/206Pb')
    
    if "207Pb/206Pb" in data.referenceMaterialData(rmName).keys():
        std_76 = data.referenceMaterialData(rmName)["207Pb/206Pb"].value()
        std_76_uncert = data.referenceMaterialData(rmName)["207Pb/206Pb"].uncertainty()
        tsd_76 = data.timeSeries("Final 207Pb/206Pb")
        measured_76 = data.groupResult(sg, tsd_76).value()
        measured_76_uncert = data.groupResult(sg, tsd_76).uncertainty()
        qaqc.pushHtml(f'<tr><td>207Pb/206Pb</td><td>{measured_76:.3f}±{measured_76_uncert:.3f}</td><td>{std_76}±{std_76_uncert}</td><td>{100*(measured_76-std_76)/std_76:.2f}%</td></tr>')
        x.append(5)
        y.append(measured_76)
        yerr.append(measured_76_uncert)
        std_y.append(std_76)
        std_yerr.append(std_76_uncert)
        ticks.append('207Pb/206Pb')
    

    tsd_208 = data.timeSeriesByMass(data.Intermediate, 208., 0.1)
    qaqc.pushHtml(f'<tr><td>208Pb Intensity</td><td>{data.groupResult(sg, tsd_208).value():.3f}</td><td></td><td></td></tr>')
    tsd_207 = data.timeSeriesByMass(data.Intermediate, 207., 0.1)
    qaqc.pushHtml(f'<tr><td>207Pb Intensity</td><td>{data.groupResult(sg, tsd_207).value():.3f}</td><td></td><td></td></tr>')
    tsd_206 = data.timeSeriesByMass(data.Intermediate, 206., 0.1)
    qaqc.pushHtml(f'<tr><td>206Pb Intensity</td><td>{data.groupResult(sg, tsd_206).value():.3f}</td><td></td><td></td></tr>')
    tsd_total204 = data.timeSeriesByMass(data.Intermediate, 204., 0.1)
    qaqc.pushHtml(f'<tr><td>Raw m204 Intensity</td><td>{data.groupResult(sg, tsd_total204).value():.4f}</td><td></td><td></td></tr>')
    tsd_202 = data.timeSeriesByMass(data.Intermediate, 202., 0.1)
    qaqc.pushHtml(f'<tr><td>202Hg Intensity</td><td>{data.groupResult(sg, tsd_202).value():.4f}</td><td></td><td></td></tr>')
    tsd_Hg204 = data.timeSeries("Hg204")
    qaqc.pushHtml(f'<tr><td>Calculated 204Hg Intensity</td><td>{data.groupResult(sg, tsd_Hg204).value():.4f}</td><td></td><td></td></tr>')
    corr204_tsd = data.timeSeries("Pb204_corr")
    qaqc.pushHtml(f'<tr><td>Corrected 204Pb Intensity</td><td>{data.groupResult(sg, corr204_tsd).value():.4f}</td><td></td><td></td></tr>')

    qaqc.pushHtml('</table>')

    graph = plot.addGraph()
    graph.setData(x,y)
    graph.setScatterStyle('ssCircle', 8, Qt.blue, Qt.blue)
    graph.setLineStyle('lsNone')

    yErrBars = QCPErrorBars(plot.bottom(), plot.left())
    yErrBars.setDataPlottable(graph)
    yErrBars.errorType = QCPErrorBars.etValueError
    yErrBars.setData(yerr)

    stdGraph = plot.addGraph()
    stdGraph.setData(x,std_y)
    stdGraph.setScatterStyle('ssSquare', 8, Qt.red, Qt.red)
    stdGraph.setLineStyle('lsNone')

    stdYErrBars = QCPErrorBars(plot.bottom(), plot.left())
    stdYErrBars.setDataPlottable(stdGraph)
    stdYErrBars.errorType = QCPErrorBars.etValueError
    stdYErrBars.setData(std_yerr)

    plot.bottom().setRange(QCPRange(0, 6))
    ratio_ticks = [i+1 for i in range(len(ticks))]
    plot.bottom().setTicks(ratio_ticks, ticks)
    plot.bottom().setTickLabelRotation(90)
    plot.bottom().label = 'Ratio'

    plot.left().rescale()
    plot.left().scaleRange(1.1)
    plot.left().label = 'Value'
    plot.replot()

    qaqc.pushFigure(plot.toPixmap(600, 400))

    if ("206Pb/204Pb" in data.referenceMaterialData(rmName).keys() and
        np.abs(100*(measured_64-std_64)/std_64) > float(qaqc.settings()["AllowableDiff_pct"]) ):
        qaqc.pushHtml('<h3>206Pb/204Pb difference is greater than the allowable % difference.</h3>')
        qaqc.finished(qaqc.Error)
        return

    if ("207Pb/204Pb" in data.referenceMaterialData(rmName).keys() and
        np.abs(100*(measured_74-std_74)/std_74) > float(qaqc.settings()["AllowableDiff_pct"]) ):
        qaqc.pushHtml('<h3>207Pb/204Pb difference is greater than the allowable % difference.</h3>')
        qaqc.finished(qaqc.Error)
        return

    if ("208Pb/204Pb" in data.referenceMaterialData(rmName).keys() and
        np.abs(100*(measured_84-std_84)/std_84) > float(qaqc.settings()["AllowableDiff_pct"]) ):
        qaqc.pushHtml('<h3>208Pb/204Pb difference is greater than the allowable % difference.</h3>')
        qaqc.finished(qaqc.Error)
        return
    
    if ("208Pb/206Pb" in data.referenceMaterialData(rmName).keys() and
        np.abs(100*(measured_86-std_86)/std_86) > float(qaqc.settings()["AllowableDiff_pct"]) ):
        qaqc.pushHtml('<h3>208Pb/206Pb difference is greater than the allowable % difference.</h3>')
        qaqc.finished(qaqc.Error)
        return
    
    if ("207Pb/206Pb" in data.referenceMaterialData(rmName).keys() and
        np.abs(100*(measured_76-std_76)/std_76) > float(qaqc.settings()["AllowableDiff_pct"]) ):
        qaqc.pushHtml('<h3>207Pb/206Pb difference is greater than the allowable % difference.</h3>')
        qaqc.finished(qaqc.Error)
        return
    
    qaqc.finished(qaqc.Success)


def settingsWidget():
    qaqc.setDefaultSetting("RefMat", "G_BCR2G")
    qaqc.setDefaultSetting("AllowableDiff_pct", 10)

    widget = QWidget()
    formLayout = QFormLayout()
    widget.setLayout(formLayout)

    refMatComboBox = QComboBox()
    refMatComboBox.addItems(data.selectionGroupNames(data.ReferenceMaterial))
    if qaqc.settings()["RefMat"] in data.selectionGroupNames(data.ReferenceMaterial):
        refMatComboBox.setCurrentText(qaqc.settings()["RefMat"])
    else:
        refMatComboBox.setCurrentText(data.selectionGroupNames(data.ReferenceMaterial)[0])
    
    qaqc.setSetting('RefMat', refMatComboBox.currentText)
    refMatComboBox.currentTextChanged.connect(lambda x: qaqc.setSetting('RefMat', x))
    formLayout.addRow('Reference material:', refMatComboBox)

    allowableDiffSpinBox = QDoubleSpinBox()
    allowableDiffSpinBox.setRange(0, 1000)
    allowableDiffSpinBox.setDecimals(2)
    allowableDiffSpinBox.setValue(qaqc.settings()["AllowableDiff_pct"])

    def updateDiffValue(x):
        qaqc.setSetting('AllowableDiff_pct', x)
        update()

    allowableDiffSpinBox.valueChanged.connect(updateDiffValue)
    formLayout.addRow('Allowable % difference:', allowableDiffSpinBox)

    qaqc.setSettingsWidget(widget)
