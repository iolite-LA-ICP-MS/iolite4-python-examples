# A python-based UI module for iolite 4 starts with some metadata
#/ Name: U-Pb Contour
#/ Authors: Albert Einstein
#/ Description: This is an example python-based plugin for UI in iolite 4
#/ References: None
#/ Version: 1.0
#/ Contact: support@iolite-software.com

from iolite.QtGui import QAction, QLabel, QSizePolicy
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

class MyMplCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        # We want the axes cleared every time plot() is called
        self.axes.hold(False)

        #
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                QSizePolicy.Expanding,
                QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

        self.axes.plot([1, 2, 3, 4])

def new_plot():
	a = MyMplCanvas()
	a.show()





def installUIHooks(window):
    print('installUIHooks for U-Pb contour')
    print(window)
    print(dir(plugin))
    IoLog.warning('installUIHooks called')
    
    a = QAction('U-Pb Contour', window)
    a.triggered.connect(new_plot)

    plugin.appendActionToMenu(["Tools"], a)