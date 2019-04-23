# A python-based QA/QC module for iolite 4 starts with some metadata
#/ Name: Fancy QA/QC
#/ Authors: Albert Einstein
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

from iolite.TimeSeriesData import TimeSeriesData
from iolite.QtGui import QWidget, QLabel, QFormLayout, QLineEdit, QComboBox, QImage
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage  

# Could specify some default settings like this:
#qaqc.setSetting('text', 'default')

def update():
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
	settings = qaqc.settings()
	
	qaqc.clearReport()
	qaqc.pushHtml('<h2>Example1</h2>')
	qaqc.pushHtml('<p>This is what your specified: <b>%s</b>.</p><br>'%(settings['text']))

	channels = [i.name for i in data.timeSeriesList(TimeSeriesData.tsInput)]
	d = [ [data.groupResult(g, data.timeSeries(c)).value() for c in channels] for g in data.selectionGroupList()]
	datadict = {g.name: v for g, v in zip(data.selectionGroupList(),d)}

	#table = pd.DataFrame(datadict, index=channels)
	#pd.set_option('precision',3)
	#qaqc.pushHtml(table.to_html())
	plt.clf()
	fig = plt.gcf()	
	fig.set_size_inches(8, 6)
	fig.set_dpi(120)
	#plt.plot(range(10))

	#fig.canvas.draw()
	#w, h = fig.canvas.get_width_height()
	#im = QImage(fig.canvas.buffer_rgba(), w, h, QImage.Format_ARGB32)	
	#qaqc.pushImage(im)


	# Try a dendogram
	X = np.array([])
	ls = []
	for c in data.timeSeriesList(TimeSeriesData.tsInput):
		print(c.name)
		
		if c.name == 'TotalBeam':
			continue
		ls += [c.name]
		d = np.array([])
		for s in data.selectionGroup(settings['text']).selections():
			sd = c.dataForSelection(s)
			#sd = data.result(s, c).value()
			d = np.insert(d, 0, sd)

		d = d/np.max(d)

		if len(X) == 0:
			X = d
		else: 
			print(X.shape)
			print(d.shape)
			X = np.column_stack((X, d))

	linked = linkage(np.transpose(X), 'ward')
	dendrogram(linked, labels=ls)
	fig.canvas.draw()
	w, h = fig.canvas.get_width_height()
	im = QImage(fig.canvas.buffer_rgba(), w, h, QImage.Format_ARGB32)	
	qaqc.pushImage(im)


	qaqc.finished(qaqc.Success)
	

def run():
	"""
	This method will be called by iolite when the user triggers the
	module from the QA/QC interface or as part of a processing template.
	It should 
	"""
	print('run called')


def settingsWidget():
	"""
	This function puts together a user interface to configure the DRS.
	
	It is important to have the last line of this function call:
	drs.setSettingsWidget(widget)
	"""
	print('settingsWidget called')
	
	widget = QWidget()
	formLayout = QFormLayout()
	widget.setLayout(formLayout)

	le = QLineEdit()
	le.connect('textChanged(QString)', lambda x: qaqc.setSetting('text', x))
	formLayout.addRow('Test', le)

	qaqc.setSettingsWidget(widget)
