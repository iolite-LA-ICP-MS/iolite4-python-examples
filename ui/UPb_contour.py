#/ Name: U-Pb Contour
#/ Authors: Joe Petrus and Bence Paul
#/ Description: This is an example that creates a contour plot of U-Pb data
#/ References: None
#/ Version: 1.0
#/ Contact: support@iolite-software.com

from iolite.QtGui import QAction, QLabel, QSizePolicy, QWidget, QComboBox, QLabel, QToolButton
from iolite.QtGui import QHBoxLayout, QVBoxLayout, QCheckBox, QImage, QPixmap, QFileDialog
from iolite.QtGui import QDialog, QFormLayout, QDialogButtonBox, QLineEdit, QPushButton
from iolite.QtCore import QDir, Qt

from matplotlib.backends.backend_qt5agg import FigureCanvas

import pandas as pd
import seaborn as sb
import numpy as np
import math

widget = None

# Decay constants
l238U = 1.55125*10**(-10)
l235U = 9.8485*10**(-10)
l232Th = 4.9475*10**(-11)
U85r = 137.818

class UPbContourWidget(QWidget):
    
    xrange_mask = (0, 30)
    yrange_mask = (0, 3)
    zrange_mask = (5000, math.inf)
    xrange_plot = (0, 30)
    yrange_plot = (0, 3)
    marker_sep = 100*10**6

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        super(QWidget, self).__init__(parent)
        self.setLayout(QVBoxLayout())

        self.header_layout = QHBoxLayout()
        self.layout().addLayout(self.header_layout)
        self.header_layout.addWidget(QLabel("Selection group"))
        self.group_combobox = QComboBox(self)
        self.header_layout.addWidget(self.group_combobox)
        self.tw_checkbox = QCheckBox(self)
        self.tw_checkbox.setText("Tera-Wasserburg?")
        self.header_layout.addWidget(self.tw_checkbox)
        
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
        self.group_combobox.currentTextChanged.connect(self.update_plot)
        self.tw_checkbox.toggled.connect(self.update_plot)
        self.setWindowTitle("U-Pb contour")


    def update_plot(self):
        self.layout().removeWidget(self.plot)
        self.plot.deleteLater()

        TW = self.tw_checkbox.checked
        group_name = self.group_combobox.currentText        
        print("UPb_contour: update_plot %s %s"%(group_name, TW))

        channel_data = {}
        channel_names = ()

        if TW:
            channel_names = ('Final U238/Pb206', 'Final Pb207/Pb206', 'Zr96_CPS')
        else:
            channel_names = ('Final Pb207/U235', 'Final Pb206/U238', 'Zr96_CPS')

        for channel_name in channel_names:
            channel = data.timeSeries(channel_name)
            channel_data[channel_name] = np.array([])

            for s in data.selectionGroup(group_name).selections():
                selection_data = channel.dataForSelection(s)
                channel_data[channel_name] = np.insert(channel_data[channel_name], 0, selection_data)

        masks = [self.xrange_mask, self.yrange_mask, self.zrange_mask]
        for i, mask in enumerate(masks):
            if mask is None:
                masks[i] = (-math.inf, math.inf)

        df = pd.DataFrame( {'x': channel_data[channel_names[0]], 'y': channel_data[channel_names[1]], 'z': channel_data[channel_names[2]]} )
        
        df = df.loc[(df['x'] > masks[0][0]) &
                    (df['x'] < masks[0][1]) &
                    (df['y'] > masks[1][0]) &
                    (df['y'] < masks[1][1]) &
                    (df['z'] > masks[2][0]) &
                    (df['z'] < masks[2][1])]

        # Calculate concordia and markers
        t = np.array(range(1, 4001*10**6, 1*10**6))
        tmarkers = np.array(range(100*10**6, 4001*10**6, self.marker_sep))
    
        Xcon, Ycon = np.array([]), np.array([])
        Xmarkers, Ymarkers = np.array([]), np.array([])
    
        if TW:        
            (Xcon, Ycon) = (1/(np.exp(l238U*t)-1), (1/U85r)*(np.exp(l235U*t)-1)/(np.exp(l238U*t)-1))
            (Xmarkers, Ymarkers) = (1/(np.exp(l238U*tmarkers)-1), (1/U85r)*(np.exp(l235U*tmarkers)-1)/(np.exp(l238U*tmarkers)-1))
        else:
            (Xcon, Ycon) = (np.exp(l235U*t)-1, np.exp(l238U*t)-1)
            (Xmarkers, Ymarkers) = (np.exp(l235U*tmarkers)-1, np.exp(l238U*tmarkers)-1)

        # Do joint plot
        sb.set_style('ticks')
        self.jointgrid = sb.jointplot(df['x'], df['y'], kind="kde", height=8, space=0,  n_levels=50, shade_lowest=False, xlim=self.xrange_plot, ylim=self.yrange_plot)    
        self.jointgrid.ax_joint.plot(Xcon, Ycon, zorder=1)
        self.jointgrid.ax_joint.scatter(Xmarkers, Ymarkers, s=50, zorder=2)
        
        if TW:
            self.jointgrid.set_axis_labels(r'${^{238}\rm{U}}/{^{206}\rm{Pb}}$', r'${^{207}\rm{Pb}}/{^{206}\rm{Pb}}$', fontsize=24)
        else:
            self.jointgrid.set_axis_labels(r'${^{207}\rm{Pb}}/{^{235}\rm{U}}$', r'${^{206}\rm{Pb}}/{^{238}\rm{U}}$', fontsize=24)
            
        x_scale = 9 * 0.78 / (self.xrange_plot[1] - self.xrange_plot[0])
        y_scale = 7 * 0.80 / (self.yrange_plot[1] - self.yrange_plot[0])

        for i, a in enumerate(tmarkers):
  
            if (a==0 or i < 2):
                continue      
            
            rotation = np.degrees(np.arctan2( (Ymarkers[i] - Ymarkers[i-1]) * y_scale, (Xmarkers[i] - Xmarkers[i-1])*x_scale))
            offset = (0, 0)
            self.jointgrid.ax_joint.annotate("%i"%(a/1e6) + " Ma", xy=(Xmarkers[i], Ymarkers[i]), xytext=offset, textcoords='offset points', rotation=rotation-90, ha='right', va='center', rotation_mode='anchor', fontsize=14)



        self.plot = FigureCanvas(self.jointgrid.fig)
        self.plot.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.layout().addWidget(self.plot)

    def save_figure(self):
        fname = QFileDialog.getSaveFileName(self, "Save figure", QDir.homePath(), "Images (*.png, *.jpg, *.tif, *.pdf, *.svg")

        print("Trying to save to file %s"%(fname))
        if fname:
            self.jointgrid.savefig(fname)

    def setup_axes(self):
        d = QDialog(self)
        l = QVBoxLayout()
        fl = QFormLayout()
        d.setLayout(l)

        params = {'X mask range': (self.xrange_mask, QLineEdit(), QLineEdit()),
                  'Y mask range': (self.yrange_mask, QLineEdit(), QLineEdit()),
                  'Z mask range': (self.zrange_mask, QLineEdit(), QLineEdit()),
                  'X plot range': (self.xrange_plot, QLineEdit(), QLineEdit()),
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
            self.xrange_mask = (float(params['X mask range'][1].text), float(params['X mask range'][2].text))
            self.yrange_mask = (float(params['Y mask range'][1].text), float(params['Y mask range'][2].text))
            self.zrange_mask = (float(params['Z mask range'][1].text), float(params['Z mask range'][2].text))
            self.xrange_plot = (float(params['X plot range'][1].text), float(params['X plot range'][2].text))
            self.yrange_plot = (float(params['Y plot range'][1].text), float(params['Y plot range'][2].text))

            self.update_plot()

def create_widget():
    global widget
    print(widget)
    if widget is None:
        widget = UPbContourWidget()    
        widget.setAttribute(Qt.WA_DeleteOnClose)    
    
    widget.show()

def installUIHooks(window):   
    a = QAction('U-Pb Contour', window)
    a.triggered.connect(create_widget)
    ui.appendActionToMenu(["Tools", "Examples"], a)