# A python-based UI module for iolite 4 starts with some metadata
#/ Type: UI
#/ Name: Excess Uncertainty Explorer
#/ Description: A tool to help explore how excess uncertainty is calculated by iolite
#/ Authors: B Paul and J Petrus
#/ References:
#/ Version: 0.1
#/ Contact:


from iolite.QtGui import QAction, QWidget, QVBoxLayout, QToolButton, QFrame
from iolite.QtGui import QTableView, QSpinBox, QDoubleSpinBox, QSlider
from iolite.QtGui import QColor, QAbstractTableModel, QAbstractItemView
from iolite.QtGui import QHeaderView, QItemSelectionModel, QPen, QBrush, QMargins
from iolite.QtCore import Qt, QSettings, QFile, QModelIndex
from iolite.QtUiTools import QUiLoader
from iolite.ui import IolitePlotPyInterface as Plot
from iolite.ui import CommonUIPyInterface as CUI
from iolite.ui import QCPErrorBars, QCPBarsGroup, QCPMarginGroup, QCP, QCPRange

import numpy as np

explorer = None


def calculate():
    """
    Calculates the MSWD and excess uncertainty for
    the current experiment

    Returns:
    ________

    splines [list of numpy arrays] - The splines that are recaclulated after
                                     each RM value is removed. A copy of
                                     index_time is stored at position 0.
                                     Splines and the copy of index_time
                                     may be downsampled as they are just for
                                     display.

    jvals [numpy array] - The jackknifed values as calculated

    spline_vals [numpy array] - The value of the spline at each jackknifed
                                value.

    rm_avgs [numpy array] - The mean value for the non-drift corrected channel
                            for each selection

    rm_uncs [numpy array] - The uncertainties of the means, as absolute 2 SE

    rm_1rse [numpy array] - The uncertainties of the means, as relative 1 SE

    rm_times [numpy array] - The selection mid-time (in seconds) for each RM
                            selection

    jval_mean [float] - The mean jackknifed value before mean-subtraction

    mswds [list of np.float64s] - The MSWD as calculated in each loop

    excesses [list of np.float64s] - The excess uncertainty as calculated in
                                    each loop

    NOTES: I have not made any effort to remove linked selections from
    the calculation. If you have linked selections in your primary
    reference material (RM) group, you may get different result to
    those calculated automatically by iolite.

    Also, this script does not stop if you have too few selections.
    By default, iolite will not calculate excess uncertainties with
    less than eight selections for the primary RM

    Also, assumes U238 is the index channel. If not, change the
    index_time line below
    """
    # Non-drift corrected channel
    ndc_name = 'DC Pb206/U238'

    # Primary reference Material
    rm_name = 'Z_91500'

    index_time = np.copy(data.timeSeries('U238').time())

    rm_group = data.selectionGroup(rm_name)
    rm_sels = rm_group.selections()
    rm_spline_type = rm_group.splineType

    ndc = data.timeSeries(ndc_name)

    rm_avgs = np.array([data.result(sel, ndc).value() for sel in rm_sels], dtype=np.float64)
    rm_1rse = np.array([data.result(sel, ndc).uncertaintyAs2SE() / (2 * data.result(sel, ndc).value()) for sel in rm_sels], dtype=np.float64)
    rm_uncs = np.array([data.result(sel, ndc).uncertaintyAs2SE() for sel in rm_sels], dtype=np.float64)
    rm_mtimes = np.array([sel.midTimestamp for sel in rm_sels], dtype=np.float64)

    print(f'rm_mtimes: {rm_mtimes}')




    # Work out JackKnifed values:
    jvals = []

    spline_vals = []
    splines = []

    for i in range(0, len(rm_avgs)):
        # create copies to remove values from
        avgs = np.copy(rm_avgs)
        uncs = np.copy(rm_uncs)
        mtimes = np.copy(rm_mtimes)

        # Record the values of this selection
        this_avg = avgs[i]
        this_mtime = mtimes[i]

        # And now remove this selection from the group of results to be splined
        avgs = np.delete(avgs, i)
        mtimes = np.delete(mtimes, i)
        uncs = np.delete(uncs, i)

        # Recalculte spline
        new_spline = data.spline(mtimes, avgs, uncs, rm_spline_type)
        data.createTimeSeries("Spline_"+str(i), data.Intermediate, index_time, new_spline)

        # And store it for later plotting
        splines.append(new_spline)

        index_val = np.where(index_time > this_mtime)[0][0]
        spline_at_index = new_spline[index_val]

        spline_vals.append(spline_at_index)

        jvals.append((this_avg - spline_at_index) / this_avg)

    # Convert JackKnifed values to numpy array:
    jvals = np.array(jvals, dtype=np.float64)
    spline_vals = np.array(spline_vals, dtype=np.float64)

    # Now kick out the last and first value because they can't be trusted
    jvals = np.delete(jvals, 0)
    jvals = np.delete(jvals, jvals.size - 1)
    rm_1rse = np.delete(rm_1rse, 0)
    rm_1rse = np.delete(rm_1rse, rm_1rse.size - 1)

    # Now mean adjust them
    mean_jval = jvals.mean()
    jvals = jvals - mean_jval

    # Place index_time copy in position 0 of splines list.
    splines.insert(0, index_time)

    '''
    Now loop through, adding some excess error (starting at 0 excess error)
    and calculate MSWD as we go. If the MSWD is less than 1.002, end loop
    '''

    MSWD = np.inf
    excessError1SE = 0.
    smallest1SE = np.min(rm_1rse)
    maxItr = 20
    mswds = []
    excesses = []

    for counter in range(maxItr):
        excesses.append(excessError1SE)

        allSels_1SEplusExcess1SE = np.sqrt(rm_1rse**2 + excessError1SE**2)
        allSels_SqWtDev = (jvals/allSels_1SEplusExcess1SE)**2
        MSWD = np.sum(allSels_SqWtDev)/(len(rm_1rse)-1)

        mswds.append(MSWD)
        # This test_mswd below is the same calc as we do for the group stats including the excess error
        # determined by the propagate errors code.
        # You can see that the group MSWD is the same as what we see in iolite, i.e. not 1 (and it matches iolite)
        rm_pe = np.sqrt( (0.5*rm_uncs)**2 + (excessError1SE*rm_avgs)**2)
        test_mswd = np.sum( (rm_avgs - np.mean(rm_avgs))**2 / (rm_pe)**2)/(len(rm_avgs)-1.)
        print(f'{counter}: MSWD = {MSWD}\t EE = {excessError1SE}, RM MSWD = {test_mswd}')

        if MSWD < 1.002:
            break

        excessError1SE += smallest1SE*(np.sqrt(MSWD)-1)

    # Record final values
    mswds.append(MSWD)
    excesses.append(excessError1SE)

    return splines, jvals, spline_vals, rm_avgs, rm_uncs, rm_1rse, rm_mtimes, mean_jval, mswds, excesses


class JackKnifeValsTableModel(QAbstractTableModel):
    '''
    Table model to handle data shown in JackKnifed Data Table
    '''

    def __init__(self, parent, rm_avgs, rm_uncs, rm_1rse, spline_vals, mean_jval, jvals):
        super().__init__(parent)
        self.rm_avgs = rm_avgs
        self.rm_uncs = rm_uncs
        self.rm_1rse = rm_1rse
        self.spline_vals = spline_vals
        self.mean_jval = mean_jval
        self.jvals = jvals

    def rowCount(self, index=QModelIndex()):
        return len(self.rm_avgs)

    def columnCount(self, index):
        return 5

    def headerData(self, section, orientation, role):

        if role == Qt.ToolTipRole:
            if section == 0:
                return 'The mean value for this selection'
            elif section == 1:
                return 'The uncertainty of this mean'
            elif section == 2:
                return 'The uncertainty expressed as 1 RSE (%)'
            elif section == 3:
                return 'The value of the spline at this point'
            elif section == 4:
                tt_str = '''The Jack-Knifed value for this selection\r
                It is calculated as (avg_val - spline_val) / avg_val
                It is then mean subtracted (the mean of the jackknifed values is subtracted from each value)\n
                '''
                tt_mean_val = f"The mean jack-knifed value for this dataset (before mean normalisation), was {self.mean_jval:.5f}"
                return tt_str + tt_mean_val

        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            if section == 0:
                return 'Avg'
            elif section == 1:
                return '2SE'
            elif section == 2:
                return '1 RSE (%)'
            elif section == 3:
                return 'Spline Value'
            elif section == 4:
                return 'Jack Knifed Value'

        elif role == Qt.DisplayRole and orientation == Qt.Vertical:
            return section + 1

        return None

    def data(self, index, role):
        if not index.isValid():
            return None

        if role != Qt.DisplayRole:
            return None

        if index.column() == 0:
            return self.rm_avgs[index.row()]

        elif index.column() == 1:
            return self.rm_uncs[index.row()]

        elif index.column() == 2:
            if index.row() == 0:
                return None
            elif index.row() > len(self.rm_1rse):
                return None
            else:
                return self.rm_1rse[index.row()-1]

        elif index.column() == 3:
            return self.spline_vals[index.row()]

        elif index.column() == 4:
            if np.isnan(self.jvals[index.row()]):
                return None
            else:
                return self.jvals[index.row()]


class ExcessUncertExplorer(QWidget):
    '''
    This widget is used to display the JackKnifed values, spline
    and other parameters used to calculate the excess uncertainty
    in iolite
    '''

    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        self.setLayout(QVBoxLayout())
        self.layout().setContentsMargins(0, 0, 0, 0)
        settings = QSettings()
        uiPath = settings.value('Paths/UIPluginsPath')
        uiFile = QFile(uiPath + '/ExcErrorExplorer.ui')
        loader = QUiLoader()
        self.ui = loader.load(uiFile)
        self.layout().addWidget(self.ui)
        self.setWindowTitle('Excess Uncertainty Explorer')

        self.jackFrame = self.ui.findChild(QFrame, 'jackFrame')
        self.jackTable = self.ui.findChild(QTableView, 'jackTableView')
        self.loopSlider = self.ui.findChild(QSlider, 'loopSlider')
        self.loopSpinBox = self.ui.findChild(QSpinBox, 'loopSpinBox')
        self.mswdSpinBox = self.ui.findChild(QDoubleSpinBox, 'mswdSpinBox')
        self.excessSpinBox = self.ui.findChild(QDoubleSpinBox, 'excessSpinBox')
        self.mswdFrame = self.ui.findChild(QFrame, 'MSWDFrame')
        self.ratioFrame = self.ui.findChild(QFrame, 'ratioFrame')
        self.barFrame = self.ui.findChild(QFrame, 'barFrame')

        splines, jvals, spline_vals, rm_avgs, rm_uncs, rm_1rse, rm_times, mean_jval, mswds, excesses = calculate()
        self.index_time = splines.pop(0)
        self.splines = splines
        self.jack_vals = jvals
        self.spline_vals = spline_vals  # This is the value of the spline at the time of the current selection. Not to be confused with the spline array
        self.rm_avgs = rm_avgs
        self.rm_uncs = rm_uncs
        self.rm_1rse = rm_1rse
        self.rm_times = rm_times
        self.mean_jval = mean_jval
        self.mswds = mswds
        self.excesses = excesses

        # Pad out jack_vals because the first and last values
        # were removed. Add blanks here to give same no of points:
        # Can't do the same with rm_1rse because it's used in the
        # calculation with the excessUncerts when updating the MSWD
        # plot
        blank_val = np.array([np.nan], dtype=np.float64)
        self.jack_vals = np.hstack([blank_val, self.jack_vals])
        self.jack_vals = np.hstack([self.jack_vals, blank_val])

        # These are the currently selected selection and loop
        # These are displayed differently
        self.currentSel = 0
        self.currentLoop = 0

        self.initUi()

    def initUi(self):
        # Set up plots
        self.jackFrame.setLayout(QVBoxLayout())
        self.jackPlot = Plot(self.jackFrame)
        self.jackFrame.layout().addWidget(self.jackPlot)
        self.jackPlot.left().label = "DC Pb206/U238"
        self.jackPlot.bottom().setDateTime(True)
        self.jackPlot.setToolsVisible(False)
        self.jackPlot.setBackground(CUI().tabBackgroundColor())

        self.mswdBoxLayout = QVBoxLayout()
        self.mswdBoxLayout.setContentsMargins(QMargins(0,0,0,0))
        self.mswdFrame.setLayout(self.mswdBoxLayout)
        self.mswdPlot = Plot(self.mswdFrame)
        self.mswdFrame.layout().addWidget(self.mswdPlot)
        self.mswdPlot.bottom().visible = False
        self.mswdPlot.left().label = "DC Pb206/U238"
        self.mswdPlot.setToolsVisible(False)
        self.mswdPlot.setBackground(CUI().tabBackgroundColor())
        self.mswdPlot.axisRect().setMargins(QMargins(0,0,0,0))

        # All left and right margins in the MSWD tab share a group to keep their margin widths the same
        self.mswdMarginGrp = QCPMarginGroup(self.mswdPlot)

        self.ratioBoxLayout = QVBoxLayout()
        self.ratioBoxLayout.setContentsMargins(QMargins(0,0,0,0))
        self.ratioFrame.setLayout(self.ratioBoxLayout)
        self.ratioPlot = Plot(self.ratioFrame)
        self.ratioFrame.layout().addWidget(self.ratioPlot)
        self.jackPlot.bottom().setDateTime(True)
        self.ratioPlot.bottom().visible = False
        self.ratioPlot.left().label = 'Jack-Knife to\ntotal uncert'
        self.ratioPlot.setToolsVisible(False)
        self.ratioPlot.setBackground(CUI().tabBackgroundColor())
        self.ratioPlotRange = None # This will be set the first time the MSWD plot is updated. It will then stay the same thereafter.
        #self.ratioPlot.axisRect().setNoAutoMargins()
        self.ratioPlot.axisRect().setMargins(QMargins(0,0,0,0))


        self.barBoxLayout = QVBoxLayout()
        self.barBoxLayout.setContentsMargins(QMargins(0,0,0,0))
        self.mswdFrame.setLayout(self.barBoxLayout)
        self.barFrame.setLayout(QVBoxLayout())
        self.barPlot = Plot(self.barFrame)
        self.barFrame.layout().addWidget(self.barPlot)
        self.barPlot.bottom().setDateTime(True)
        self.barPlot.setToolsVisible(False)
        self.barPlot.setLegendVisible(True)
        self.barPlot.setBackground(CUI().tabBackgroundColor())

        # TODO: Margin groups aren't meant to be used between different plots,
        #       they're meant to be used on different layout elements within
        #       one plot. If we want these aligned it will required changing
        #       the ui file and getting rid of the three frames and instead
        #       having 1 frame with 1 plot and having 3 QCPAxisRects in the
        #       single plot.
        self.mswdPlot.axisRect().setMarginGroup(QCP.msLeft | QCP.msRight, self.mswdMarginGrp)
        self.ratioPlot.axisRect().setMarginGroup(QCP.msLeft | QCP.msRight, self.mswdMarginGrp)
        self.barPlot.axisRect().setMarginGroup(QCP.msLeft | QCP.msRight, self.mswdMarginGrp)

        # Set up table of JackKnifed Values
        self.jModel = JackKnifeValsTableModel(self, self.rm_avgs, self.rm_uncs, self.rm_1rse, self.spline_vals, self.mean_jval, self.jack_vals)
        self.jackTable.setModel(self.jModel)
        self.jackTable.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.jackTable.selectRow(0)
        self.jackTable.horizontalHeader().setSectionResizeMode(
            QHeaderView.Stretch)
        self.jackTable.selectionModel().selectionChanged.connect(
            self.process_table_click)

        # Set up loop slider
        self.loopSlider.setMinimum(0)
        self.loopSlider.setMaximum(len(self.mswds) - 1)
        self.loopSlider.valueChanged.connect(self.updateLoopNumber)
        self.loopSpinBox.setMinimum(1)
        self.loopSpinBox.setMaximum(len(self.mswds))
        # This is the little clunky, but the spin box shows the loop number + 1
        # Because the updateLoopNumber function works for both the slider and the
        # spinBox, we just have to subtract one from the value for the spin box
        self.loopSpinBox.valueChanged.connect(lambda v: self.updateLoopNumber(v-1))

        self.resize(900, 850)

        self.updateJackPlot()
        self.updateMSWDPlot()

    def updateLoopNumber(self, l):
        self.currentLoop = l
        self.loopSpinBox.blockSignals(True)
        self.loopSlider.blockSignals(True)
        self.updateMSWDPlot()
        self.loopSpinBox.blockSignals(False)
        self.loopSlider.blockSignals(False)

    def process_table_click(self, selectedItems, deselectedItems):
        if len(selectedItems.indexes()) < 1:
            return

        sel_row = selectedItems.indexes()[0].row()
        self.update_current_selection(sel_row)

    def update_current_selection(self, sel_num):
        if sel_num < 0:
            return
        if sel_num > len(self.rm_avgs) - 1:
            return

        self.currentSel = sel_num
        self.updateJackPlot()

    def updateJackPlot(self):

        self.jackPlot.clearGraphs()

        for si, spline in enumerate(self.splines):
            g = self.jackPlot.addGraph()
            g.setData(self.index_time, spline)
            if si == self.currentSel:
                continue  # Skip this: it's plotted last so it's on top
            else:
                g.setColor(QColor(Qt.gray))

        # Add spline for current selection
        g = self.jackPlot.addGraph()
        g.setData(self.index_time, self.splines[self.currentSel])
        splinePen = QPen()
        splinePen.setColor(Qt.red)
        splinePen.setWidthF(4)
        g.setPen(splinePen)

        # Now plot the offset between mean value and spine:
        self.jackPlot.clearItems()
        line = self.jackPlot.addLine([self.rm_times[self.currentSel],
                                      self.rm_avgs[self.currentSel]],
                                     [self.rm_times[self.currentSel],
                                      self.spline_vals[self.currentSel]],
                                     'esFlatArrow', 'esNone')

        line.pen = QPen(Qt.black, 3, Qt.SolidLine, Qt.RoundCap, Qt.RoundJoin)

        # Add Avgs and uncertainties
        g_avgs = self.jackPlot.addGraph()
        g_avgs.setData(self.rm_times, self.rm_avgs)
        g_avgs.setLineStyle('lsNone')
        g_avgs.setScatterStyle('ssCircle', 6, QColor(Qt.gray), QColor(Qt.gray))

        eby = QCPErrorBars(self.jackPlot.bottom(), self.jackPlot.left())
        eby.setData(self.rm_uncs)
        eby.setDataPlottable(g_avgs)
        eby.errorType = QCPErrorBars.etValueError

        # Now plot the current selection
        g_this_avg = self.jackPlot.addGraph()
        this_x = np.array([self.rm_times[self.currentSel]], dtype=np.float64)
        this_y = np.array([self.rm_avgs[self.currentSel]], dtype=np.float64)
        g_this_avg.setData(this_x, this_y)
        g_this_avg.setLineStyle('lsNone')
        g_this_avg.setScatterStyle('ssCircle', 6, QColor(Qt.red), QColor(Qt.red))

        self.jackPlot.rescaleAxes()
        self.jackPlot.bottom().scaleRange(1.1)
        self.jackPlot.left().scaleRange(1.1)
        self.jackPlot.replot()

    def updateMSWDPlot(self):
        '''
        This function updates both the mswdPlot and the barPlot
        '''

        self.loopSpinBox.setValue(self.currentLoop + 1)
        self.loopSlider.setValue(self.currentLoop)
        self.mswdSpinBox.setValue(self.mswds[self.currentLoop])
        self.excessSpinBox.setValue(self.excesses[self.currentLoop])

        self.mswdPlot.clearGraphs()

        g_avgs = self.mswdPlot.addGraph()
        g_avgs.setData(self.rm_times, self.rm_avgs)
        g_avgs.setLineStyle('lsNone')
        g_avgs.setScatterStyle('ssCircle', 6, QColor(Qt.gray), QColor(Qt.gray))

        eby = QCPErrorBars(self.mswdPlot.bottom(), self.mswdPlot.left())
        eby.setData(self.rm_uncs)
        eby.setDataPlottable(g_avgs)
        eby.errorType = QCPErrorBars.etValueError

        allSels_1SEplusExcess1SE = np.sqrt(self.rm_1rse**2 +
                                           self.excesses[self.currentLoop]**2)

        blank_val = np.array([np.nan], dtype=np.float64)
        allSels_1SEplusExcess1SE = np.hstack([blank_val, allSels_1SEplusExcess1SE])
        allSels_1SEplusExcess1SE = np.hstack([allSels_1SEplusExcess1SE, blank_val])

        totalUnc = allSels_1SEplusExcess1SE * self.rm_avgs * 2

        eby_tot = QCPErrorBars(self.mswdPlot.bottom(), self.mswdPlot.left())
        eby_tot.setData(totalUnc)
        eby_tot.setDataPlottable(g_avgs)
        eby_tot.whiskerWidth = 13
        eby_tot.pen = QPen(Qt.black, 2, Qt.SolidLine, Qt.RoundCap, Qt.RoundJoin)
        eby_tot.errorType = QCPErrorBars.etValueError

        self.mswdPlot.rescaleAxes()

        # Working out min and max y ranges here so that the y range can be fixed:
        maxUnc = np.nanmax(np.sqrt(self.rm_1rse**2 + np.nanmax(self.excesses)**2))
        yMax = np.nanmax(self.rm_avgs) + maxUnc
        yMin = np.nanmin(self.rm_avgs) - maxUnc

        self.mswdPlot.left().range = QCPRange(yMin, yMax)
        self.mswdPlot.bottom().scaleRange(1.1)
        self.mswdPlot.replot()

        # Will use the range of the bottom axis of the mswd plot a couple of times here
        # so assign it to it's own variables:
        timeMin = self.mswdPlot.bottom().range.lower()
        timeMax = self.mswdPlot.bottom().range.upper()

        # Show ratios in their own little plot:
        self.ratioPlot.clearGraphs()
        self.ratioPlot.clearPlottables()
        self.ratioPlot.clearItems()

        uncs = np.sqrt(self.rm_1rse**2 + self.excesses[self.currentLoop]**2)
        uncs = np.hstack([blank_val, uncs])
        uncs = np.hstack([uncs, blank_val])

        local_jacks = np.copy(self.jack_vals)
        ratios = (local_jacks / uncs)**2

        g_ratios = self.ratioPlot.addGraph()
        g_ratios.setData(self.rm_times, ratios)
        g_ratios.setLineStyle('lsNone')
        g_avgs.setScatterStyle('ssCircle', 6, QColor(Qt.gray), QColor(Qt.gray))
        g_ratios.setScatterStyle('ssSquare', 6, QColor(144, 219, 244),
                                 QColor(144, 219, 244))

        mean_ratio = np.sum(ratios[np.isfinite(ratios)]) / (len(ratios[np.isfinite(ratios)]) - 1)

        rmean_line = self.ratioPlot.addStraightLine([timeMin, mean_ratio],
                                                    [timeMax, mean_ratio])
        rmean_pen = QPen()
        rmean_pen.setColor(QColor(Qt.gray))
        rmean_pen.setWidthF(2)
        rmean_pen.setStyle(Qt.DashLine)
        rmean_line.pen = rmean_pen

        self.ratioPlot.rescaleAxes()
        self.ratioPlot.bottom().range = QCPRange(timeMin, timeMax)

        if self.ratioPlotRange is None:
            self.ratioPlot.left().scaleRange(1.1)
            self.ratioPlotRange = self.ratioPlot.left().range

        self.ratioPlot.left().range = self.ratioPlotRange

        self.ratioPlot.replot()

        # Bar Plot
        self.barPlot.clearPlottables()
        self.barPlot.clearItems()

        bars_group = QCPBarsGroup(self.barPlot)
        j_bars = self.barPlot.addBars(self.rm_times, np.abs(self.jack_vals))
        j_bars.name = "Jack-Knifed Values"
        j_bars.pen = QPen(QColor(163, 196, 243))
        j_bars.brush = QBrush(QColor(163, 196, 243))
        j_bars.width = 50
        j_bars.barsGroup = bars_group
        rse_bars = self.barPlot.addBars(self.rm_times, allSels_1SEplusExcess1SE)
        rse_bars.name = 'Total Uncertainty (2 SE)'
        rse_bars.pen = QPen(QColor(185, 251, 192))
        rse_bars.brush = QBrush(QColor(185, 251, 192))
        rse_bars.width = 50
        rse_bars.barsGroup = bars_group

        self.barPlot.rescaleAxes()
        self.barPlot.bottom().range = QCPRange(timeMin, timeMax)
        self.barPlot.replot()


def createUIElements():
    action = QAction('Excess Uncertainty Explorer', None)
    action.triggered.connect(show_explorer)
    ui.setAction(action)
    ui.setMenuName(['Tools'])


def show_explorer():
    global explorer
    explorer = ExcessUncertExplorer()
    explorer.setAttribute(Qt.WA_DeleteOnClose)
    explorer.show()
