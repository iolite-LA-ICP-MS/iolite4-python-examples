# A python-based QA/QC module for iolite 4 starts with some metadata
#/ Type: QAQC
#/ Name: QA/QC Introduction
#/ Authors: Joe Petrus and Bence Paul
#/ Description: This is an example python-based plugin for QA/QC in iolite 4
#/ References: None
#/ Version: 1.0
#/ Contact: support@iolite-software.com

"""
Qt imports can be done through 'iolite', e.g.
from iolite.QtGui import QLabel

Before functions are called a few additional objects are
added to the module:

data	an interface to iolite's C++ data. E.g. you can get
        existing time series data or selection groups with it
        as well as make new ones.

IoLog	an interface to iolite's logging facility. You can add
        messages with, e.g., IoLog.debug('My message')

qaqc	an interface to the PythonQAQC C++ class in iolite from 
		which some built-in features can be accessed, e.g.,
		pushHtml(html)
"""

from iolite.QtGui import QWidget, QLabel, QFormLayout, QLineEdit, QComboBox, QImage
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

def update():
	"""
	This method will be called by iolite when the user triggers
	the module to run in the QA/QC UI or as part of a processing
	template. It should analyze current results to evaluate whether
	the instruments were performing as expected.

	When finished, the 'finished' signal should be emitted with a status code,
	e.g. qaqc.Success, qaqc.Error, qaqc.SuccessWithWarnings 

	As an example, we will compare the group stats of the specified
	selection group and channel to a target value. If they're within
	uncertainty, finish with success.

	"""	

	# Add the QA/QC settings to the report:
	settings = qaqc.settings()
	
	qaqc.clearReport()
	qaqc.pushHtml('<h2>QA/QC Introduction</h2>')
	qaqc.pushHtml('<p>These are the settings you supplied:<br>')
	for key in settings:
		qaqc.pushHtml('<b>%s</b>: %s<br>'%(key, settings[key]))
	qaqc.pushHtml('</p>')

	# Get the specified channel and group:
	group = data.selectionGroup(settings["GroupName"])
	channel = data.timeSeries(settings["ChannelName"])
	target = float(settings["Target"])
	allowable_diff = float(settings["AllowableDiff"])
	target_diff = target * (allowable_diff/100)

	plt.clf()
	fig = plt.gcf()	
	fig.set_size_inches(6, 4)
	fig.set_dpi(120)
	ax = fig.add_subplot(1, 1, 1)

	# Plot the individual selection data:
	x = range(0, len(group.selections()))
	y = [data.result(s, channel).value() for s in group.selections()]
	w = [data.result(s, channel).uncertaintyAs2SE() for s in group.selections()]

	plt.errorbar(x, y, yerr=w, fmt='o', capsize=5, color='b')

	# Plot group result:
	gr = data.groupResult(group, channel)
	group_patch = Rectangle( (-0.5, gr.value() - gr.uncertainty()), len(group.selections()), 2*gr.uncertainty(), color='b', alpha=0.3)
	ax.add_patch(group_patch)

	# Plot target:
	target_patch = Rectangle( (-0.5, target-target_diff), len(group.selections()), 2*target_diff, color='g', alpha=0.3)
	ax.add_patch(target_patch)

	fig.canvas.draw()
	w, h = fig.canvas.get_width_height()
	im = QImage(fig.canvas.buffer_rgba(), w, h, QImage.Format_ARGB32)	
	qaqc.pushImage(im)

	grmin, grmax = gr.value() - gr.uncertainty(), gr.value() + gr.uncertainty()
	trmin, trmax = target-target_diff, target+target_diff
	
	# Check if the group result is within the target:
	if grmin < trmax and grmax > trmin:
		qaqc.finished(qaqc.Success)
	else:
		qaqc.finished(qaqc.Error)


def settingsWidget():
	"""
	This function puts together a user interface to configure the QA/QC module.
	The widget should be kept as small as possible and will be shown as a popup
	when the user clicks the settings button.

	It is important to have the last line of this function call:
	qaqc.setSettingsWidget(widget)
	"""	

	qaqc.setDefaultSetting("GroupName", data.selectionGroupNames()[0])
	qaqc.setDefaultSetting("ChannelName", data.timeSeriesNames()[0])
	qaqc.setDefaultSetting("Target", 100)
	qaqc.setDefaultSetting("AllowableDiff", 10)

	widget = QWidget()
	formLayout = QFormLayout()
	widget.setLayout(formLayout)

	group_lineedit = QLineEdit()
	group_lineedit.setText(qaqc.settings()["GroupName"])
	group_lineedit.textChanged.connect(lambda x: qaqc.setSetting('GroupName', x))
	formLayout.addRow('Group name:', group_lineedit)

	channel_lineedit = QLineEdit()
	channel_lineedit.setText(qaqc.settings()["ChannelName"])
	channel_lineedit.textChanged.connect(lambda x: qaqc.setSetting('ChannelName', x))
	formLayout.addRow('Channel name:', channel_lineedit)

	target_lineedit = QLineEdit()
	target_lineedit.setText(qaqc.settings()["Target"])
	target_lineedit.textChanged.connect(lambda x: qaqc.setSetting('Target', float(x)))
	formLayout.addRow('Target value:', target_lineedit)

	allowable_diff_lineedit = QLineEdit()
	allowable_diff_lineedit.setText(qaqc.settings()["AllowableDiff"])
	allowable_diff_lineedit.textChanged.connect(lambda x: qaqc.setSetting('AllowableDiff', float(x)))
	formLayout.addRow('Allowable % diff:', allowable_diff_lineedit)

	qaqc.setSettingsWidget(widget)
