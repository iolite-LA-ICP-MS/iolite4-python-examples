#/ Type: UI
#/ Name: Dendrogram
#/ Authors: Joe Petrus and Bence Paul
#/ Description: Creates a dendrogram for the specified group using matplotlib/scipy
#/ References: None
#/ Version: 1.0
#/ Contact: support@iolite-software.com

from iolite.QtGui import QAction, QWidget, QComboBox, QLabel, QHBoxLayout, QVBoxLayout, QImage, QPixmap, QSizePolicy
from iolite.QtCore import Qt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvas
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage  

widget = None

def createUIElements():   
	a = QAction('Dendrogram', ui)
	a.triggered.connect(create_widget)
	ui.setAction(a)
	ui.setMenuName(['Examples'])


class DendrogramWidget(QWidget):

	channel_type = data.Input
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
		self.plot = QLabel("Select a group above...")
		self.layout().addWidget(self.plot)
		self.resize(600, 800)
		self.setWindowTitle("Dendrogram Example")
		self.show()

	def update_dendrogram(self):
		self.layout().removeWidget(self.plot)
		self.plot.deleteLater()

		group_name = self.group_combobox.currentText
		print("dendrogram: update_dendrogram %s"%(group_name))
		
		plt.clf()
		fig = plt.gcf()	
		fig.set_size_inches(8, 8)
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

			try:
				d = d/np.max(d)
			except:
				labels.pop()
				print('Error processing channel %s'%(channel.name))
				continue

			if len(X) == 0:
				X = d
			else: 
				X = np.column_stack((X, d))

		try:
			linked = linkage(np.transpose(X), 'ward')
			dendrogram(linked, labels=labels)
		except RuntimeError as r:
			print(r)
			self.plot = QLabel('There was a problem processing the selected group.')
			self.layout().addWidget(self.plot)
			return

		fig.canvas.draw()
		self.plot = FigureCanvas(fig)
		self.plot.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
		self.layout().addWidget(self.plot)

def create_widget():
	print("dendrogram: create_widget")
	global widget
	try:
		widget.show()
	except:	
		widget = DendrogramWidget()
		widget.setAttribute(Qt.WA_DeleteOnClose)
		widget.show()
	