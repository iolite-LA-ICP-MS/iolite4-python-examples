#/ Name: Density Plot
#/ Authors: Joe Petrus and Bence Paul
#/ Description: This is an example that creates a density plot
#/ References: None
#/ Version: 1.0
#/ Contact: support@iolite-software.com

from iolite.QtGui import QAction, QLabel, QSizePolicy, QWidget, QComboBox, QLabel, QToolButton
from iolite.QtGui import QHBoxLayout, QVBoxLayout, QCheckBox, QImage, QPixmap, QFileDialog
from iolite.QtGui import QDialog, QFormLayout, QDialogButtonBox, QLineEdit, QPushButton
from iolite.QtCore import QDir, Qt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvas
import numpy as np
from scipy import stats

widget = None

def kde_bandwidth(obj, fac=1./5):
    return np.power(obj.n, -1./(obj.d+4))*fac

class DensityWidget(QWidget):
    
    xrange_plot = None
    yrange_plot = None

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        super(QWidget, self).__init__(parent)
        self.setLayout(QVBoxLayout())

        self.header_layout = QHBoxLayout()
        self.layout().addLayout(self.header_layout)
        
        self.header_layout.addWidget(QLabel("Selection group"))
        self.group_combobox = QComboBox(self)        
        self.header_layout.addWidget(self.group_combobox)
        
        self.header_layout.addWidget(QLabel("Channel"))        
        self.channel_combobox = QComboBox(self)
        self.header_layout.addWidget(self.channel_combobox)
        
        axes_button = QToolButton(self)
        axes_button.setText("Setup axes")
        axes_button.clicked.connect(self.setup_axes)
        self.header_layout.addWidget(axes_button)


        spacer = QWidget(self)
        spacer.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Minimum)
        self.header_layout.addWidget(spacer)
        
        save_button = QToolButton(self)
        save_button.setText("Save")
        save_button.clicked.connect(self.save_figure)
        self.header_layout.addWidget(save_button)

        self.plot = QLabel("Select a group above...")
        self.layout().addWidget(self.plot)

        self.group_combobox.addItems(data.selectionGroupNames())
        self.group_combobox.currentTextChanged.connect(lambda: self.update_plot())
        self.channel_combobox.addItems(data.timeSeriesNames())
        self.channel_combobox.currentTextChanged.connect(lambda: self.update_plot())
        self.setWindowTitle("Density Plot")


    def update_plot(self, xrange=None, yrange=None):
        self.layout().removeWidget(self.plot)
        self.plot.deleteLater()
        group_name = self.group_combobox.currentText
        channel_name = self.channel_combobox.currentText
        print("density: update_plot %s %s"%(group_name, channel_name))

        group = data.selectionGroup(group_name)
        channel = data.timeSeries(channel_name)        
        
        x = np.empty(group.count)*np.nan
        err = np.empty(group.count)*np.nan

        for i, s in enumerate(group.selections()):
            r = data.result(s, channel)
            x[i] = r.value()
            err[i] = r.uncertaintyAs2SE()
            

        self.fig = plt.figure()
        self.fig.set_size_inches(8, 8)
        self.fig.set_dpi(120)
        ax = self.fig.add_subplot(111)

        kde = stats.gaussian_kde(x, bw_method=kde_bandwidth)

        ax.plot(x, np.zeros(x.shape), 'b|', ms=5)

        x_eval = np.linspace(0.8*x.min(), 1.2*x.max(), num=2000)
        ax.plot(x_eval, kde(x_eval), 'k-', label='kde')

        if xrange is None:
            self.xrange_plot = ax.get_xlim()
        else:
            print('xrange=%s'%(xrange))
            ax.set_xlim(xrange)

        if yrange is None:
            self.yrange_plot = ax.get_ylim()
        else:
            print('yrange=%s'%(yrange))
            ax.set_ylim(yrange)

        ax.set_xlabel(channel.name)

        self.fig.canvas.draw()
        self.plot = FigureCanvas(self.fig)
        self.plot.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.layout().addWidget(self.plot)

    def save_figure(self):
        fname = QFileDialog.getSaveFileName(self, "Save figure", QDir.homePath(), "Images (*.png, *.jpg, *.tif, *.pdf, *.svg")

        print("Trying to save to file %s"%(fname))
        if fname:
            self.fig.savefig(fname)

    def setup_axes(self):
        d = QDialog(self)
        l = QVBoxLayout()
        fl = QFormLayout()
        d.setLayout(l)

        params = {'X plot range': (self.xrange_plot, QLineEdit(), QLineEdit()),
                  'Y plot range': (self.yrange_plot, QLineEdit(), QLineEdit())}

        for key in params:
            ll = QHBoxLayout()            
            params[key][1].setText(params[key][0][0])
            params[key][2].setText(params[key][0][1])
            ll.addWidget(params[key][1])
            ll.addWidget(params[key][2])
            fl.addRow(key, ll)

        l.addLayout(fl)

        ok_button = QPushButton('Ok')
        cancel_button = QPushButton('Cancel')
        button_layout = QHBoxLayout()
        button_layout.addWidget(cancel_button)
        button_layout.addWidget(ok_button)        
        l.addLayout(button_layout)

        ok_button.clicked.connect(lambda: d.accept())
        cancel_button.clicked.connect(lambda: d.reject())

        if (d.exec() == QDialog.Accepted):
            self.xrange_plot = (float(params['X plot range'][1].text), float(params['X plot range'][2].text))
            self.yrange_plot = (float(params['Y plot range'][1].text), float(params['Y plot range'][2].text))

            self.update_plot(xrange=self.xrange_plot, yrange=self.yrange_plot)

def create_widget():
    global widget
    if widget is None:
        widget = DensityWidget()    
        widget.setAttribute(Qt.WA_DeleteOnClose)

    widget.show()  
    

def installUIHooks(window):   
    a = QAction('Density Plot', window)
    a.triggered.connect(create_widget)
    ui.appendActionToMenu(["Tools", "Examples"], a)