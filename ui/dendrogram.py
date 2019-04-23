#/ Name: Dendrogram
#/ Authors: Joe Petrus and Bence Paul
#/ Description: Creates a dendrogram for the specified group using matplotlib/scipy
#/ References: None
#/ Version: 1.0
#/ Contact: support@iolite-software.com

from iolite.QtGui import QAction, QWidget, QComboBox, QLabel, QHBoxLayout, QVBoxLayout, QImage, QPixmap
from iolite.TimeSeriesData import TimeSeriesData
from iolite.SelectionGroup import SelectionGroup
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage  

widget = None

def installUIHooks(window):   
    a = QAction('Dendrogram', window)
    a.triggered.connect(create_widget)
    plugin.appendActionToMenu(["Tools", "Examples"], a)

class DendrogramWidget(QWidget):

	channel_type = TimeSeriesData.tsInput
	image_widget = None
	group_combobox = None

	def __init__(self, parent = None):
		super(QWidget, self).__init__(parent)

		self.setLayout(QVBoxLayout())
	
		header_layout = QHBoxLayout()	
		header_layout.addWidget(QLabel("Selection group:"))
	
		self.group_combobox = QComboBox()
		self.group_combobox.addItems(data.selectionGroupNames())
		self.group_combobox.currentTextChanged.connect(self.update_dendrogram)
		header_layout.addWidget(self.group_combobox)
	
		self.layout().addLayout(header_layout)
		self.image_widget = QLabel()
		self.layout().addWidget(self.image_widget)
		self.resize(600, 800)
		self.setWindowTitle("Dendrogram Example")
		self.show()

	def update_dendrogram(self):
		group_name = self.group_combobox.currentText
		print("dendrogram: update_dendrogram %s"%(group_name))
		
		plt.clf()
		fig = plt.gcf()	
		fig.set_size_inches(8, 6)
		fig.set_dpi(120)

		X = np.array([])
		labels = []
		
		for channel in data.timeSeriesList(self.channel_type):
			if channel.name == "TotalBeam":
				continue

			labels += [channel.name]
			d = np.array([])
			for s in data.selectionGroup(group_name).selections():
				selection_data = channel.dataForSelection(s)		
				d = np.insert(d, 0, selection_data)

			d = d/np.max(d)

			if len(X) == 0:
				X = d
			else: 
				X = np.column_stack((X, d))

		linked = linkage(np.transpose(X), 'ward')
		dendrogram(linked, labels=labels)
		fig.canvas.draw()
		w, h = fig.canvas.get_width_height()
		im = QImage(fig.canvas.buffer_rgba(), w, h, QImage.Format_ARGB32)

		self.image_widget.setPixmap(QPixmap.fromImage(im))


def create_widget():
	print("dendrogram: create_widget")
	global widget
	widget = DendrogramWidget()
	widget.show()
	