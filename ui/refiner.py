#/ Type: UI
#/ Name: Selection Refiner
#/ Description: Refine selections using multiple algorithms
#/ Authors: J. Petrus and B. Paul
#/ References:
#/ Version: 1.0
#/ Contact: support@iolite-software.com

from iolite.QtGui import QWidget, QMessageBox, QAction, QVBoxLayout, QToolButton, QDialogButtonBox, QListWidget, QTableView, QGroupBox, QLineEdit, QComboBox, QFormLayout
from iolite.QtGui import QHBoxLayout, QProgressBar
from iolite.QtUiTools import QUiLoader
from iolite.QtCore import QFile, QSettings, QDateTime, Qt
from iolite.ui import CommonUIPyInterface as CUI
from iolite.ui import ASCriteriaModel, ComboBoxDelegate, ComparatorDelegate, ChannelsComboBox
from sklearn.cluster import MeanShift
from sklearn import preprocessing
from scipy.stats import sem
import numpy as np

rw = None

def createUIElements():
    action = QAction('Refine Selections', None)
    action.triggered.connect(main)
    ui.setAction(action)
    ui.setMenuName(['Selections'])

# Run length encoding from stack overflow :)
def rle(inarray):
    ia = np.asarray(inarray)
    n = len(ia)
    if n == 0:
        return (None, None, None)
    else:
        y = np.array(ia[1:] != ia[:-1])     # pairwise unequal (string safe)
        i = np.append(np.where(y), n - 1)   # must include last element posi
        z = np.diff(np.append(-1, i))       # run lengths
        p = np.cumsum(np.append(0, z))[:-1] # positions
        return(z, p, ia[i])

class AbstractStrategy(object):
    usesGlobalTime = False

    def settings(self):
        if not hasattr(self, 'settingsWidget'):
            self.settingsWidget = self.createSettingsWidget()

        return self.settingsWidget


class CriteriaStrategy(AbstractStrategy):
    name = 'Criteria'

    def __init__(self):
        print('Criteria strategy')

    def process(self, selection):
        t = data.timeSeries('TotalBeam').timeForSelection(selection)
        mask = np.ones(len(t))

        for c in self.criteriaModel.criteriaAsVariantList():
            print(c)
            d = data.timeSeries(c['channel']).dataForSelection(selection)
            if c['comparator'] == 'GreaterThan':
                mask *= d > c['firstvalue']
            elif c['comparator'] == 'LessThan':
                mask *= d < c['firstvalue']
            elif c['comparator'] == 'InRange':
                mask *= (d > c['firstvalue'])*(d < c['secondvalue'])
            elif c['comparator'] == 'NotInRange':
                mask *= (d < c['firstvalue']) + (d > c['secondvalue'])
            elif c['comparator'] == 'EqualTo':
                mask *= d == c['firstvalue']
            elif c['comparator'] == 'NotEqualTo':
                mask *= d != c['firstvalue']

        idx_pairs = np.where(np.diff(np.hstack(([False],mask==1,[False]))))[0].reshape(-1,2)
        return idx_pairs

    def createSettingsWidget(self):
        w = QWidget()
        w.setLayout(QFormLayout())
        w.layout().setContentsMargins(0, 0, 0, 0)

        self.criteriaModel = ASCriteriaModel()
        self.tableView = QTableView(w)
        self.tableView.setModel(self.criteriaModel)
        self.tableView.setItemDelegateForColumn(0, ComboBoxDelegate(data.timeSeriesNames(), w))
        self.tableView.setItemDelegateForColumn(1, ComparatorDelegate(w))
        buttonLayout = QHBoxLayout()
        self.addButton = QToolButton(w)
        self.addButton.setIcon(CUI().icon('plus'))
        self.addButton.clicked.connect(lambda: self.criteriaModel.insertRows(0, 1))
        self.removeButton = QToolButton(w)
        self.removeButton.setIcon(CUI().icon('minus'))
        self.removeButton.clicked.connect(lambda: self.criteriaModel.removeRow(self.tableView.currentIndex().row()))
        buttonLayout.addWidget(self.addButton)
        buttonLayout.addWidget(self.removeButton)
        w.layout().addRow(buttonLayout)
        w.layout().addRow(self.tableView)
        return w

class MeanShiftStrategy(AbstractStrategy):
    name = 'Mean shift clustering'

    def process(self, selection):
        import pandas as pd
        channelNames = self.channels.split(';')

        d = {}
        for cn in channelNames:
            d[cn] = data.timeSeries(cn).dataForSelection(selection)
        df = pd.DataFrame(d)

        # Scale the data? Not sure if this is required.
        x = df.values
        min_max_scalar = preprocessing.MinMaxScaler()
        x_scaled = min_max_scalar.fit_transform(x)
        df = pd.DataFrame(x_scaled)

        cc = MeanShift()
        clustering = cc.fit(df)
        n_clusters = len(np.unique(clustering.labels_))

        # Calculate longest stretches of same label
        z, p, _ = rle(clustering.labels_)
        max_rl = np.amax(z)
        start_index = p[np.argmax(z)]

        return np.array([[start_index, start_index + max_rl]])

    def createSettingsWidget(self):
        w = QWidget()
        w.setLayout(QFormLayout())
        w.layout().setContentsMargins(0, 0, 0, 0)
        self.channelsComboBox = ChannelsComboBox(w)
        self.channelsComboBox.setExclusive(False)
        self.channelsComboBox.channelsChanged.connect(lambda: setattr(self, 'channels', self.channelsComboBox.channelsString()))
        w.layout().addRow('Channels', self.channelsComboBox)
        return w


class RollingFuncStrategy(AbstractStrategy):
    # Not to be used directly...
    def __init__(self):
        self.func = None
        self.start = None

    def process(self, selection):
        channelName = self.channel
        mink = self.mink

        d = data.timeSeries(channelName).dataForSelection(selection)
        ct = self.start
        idx = None

        for k in range(mink, len(d)):
            for i in range(0, len(d)-k):
                t = self.func(d[i:i+k])
                if t <= ct:
                    ct = t
                    idx = [i, i+k]
        if idx:
            return np.array([idx])
        else:
            return []    

class RollingSDStrategy(RollingFuncStrategy):
    name = 'Rolling standard deviation'

    def __init__(self):
        super().__init__()
        self.func = np.std
        self.start = np.inf     

    def createSettingsWidget(self):
        w = QWidget()
        w.setLayout(QFormLayout())
        w.layout().setContentsMargins(0, 0, 0, 0)
        self.channelsComboBox = ChannelsComboBox(w)
        self.channelsComboBox.setExclusive(True)
        self.channelsComboBox.channelsChanged.connect(lambda: setattr(self, 'channel', self.channelsComboBox.channelsString()))
        w.layout().addRow('Channels', self.channelsComboBox)
        self.minkLineEdit = QLineEdit(w)
        self.minkLineEdit.textEdited.connect(lambda: setattr(self, 'mink', int(self.minkLineEdit.text)))
        w.layout().addRow('Minimum # of points', self.minkLineEdit)
        return w

class RollingSEStrategy(RollingSDStrategy):
    name = 'Rolling standard error'

    def __init__(self):
        super().__init__()
        self.func = sem

class ExtremeStrategy(RollingFuncStrategy):
    name = 'Extremes'

    def __init__(self):
        super().__init__()
        self.start = np.inf
        self.updateFunc('Maximize')

    def updateFunc(self, func):
        if 'Max' in func:
            self.func = lambda v: -np.mean(v) # Since RollingFuncStrategy is trying to minimize...
        else:
            self.func = np.mean
            
    def createSettingsWidget(self):
        w = QWidget()
        w.setLayout(QFormLayout())
        w.layout().setContentsMargins(0, 0, 0, 0)
        self.channelsComboBox = ChannelsComboBox(w)
        self.channelsComboBox.setExclusive(True)
        self.channelsComboBox.channelsChanged.connect(lambda: setattr(self, 'channel', self.channelsComboBox.channelsString()))
        w.layout().addRow('Channels', self.channelsComboBox)
        self.minkLineEdit = QLineEdit(w)
        self.minkLineEdit.textEdited.connect(lambda: setattr(self, 'mink', int(self.minkLineEdit.text)))
        w.layout().addRow('Minimum # of points', self.minkLineEdit)
        self.minmaxComboBox = QComboBox(w)
        self.minmaxComboBox.addItems(['Maximize', 'Minimize'])
        self.minmaxComboBox.activated.connect(lambda: self.updateFunc(self.minmaxComboBox.currentText))
        w.layout().addRow('Extreme', self.minmaxComboBox)
        return w       

class SyncStrategy(AbstractStrategy):
    name = 'Synchronize'
    usesGlobalTime = True

    def __init__(self):
        super().__init__()
        self.maxOffset = 10

    def process(self, selection):
        tb = data.timeSeries('TotalBeam')
        ts = tb.time()[1] - tb.time()[0]
        n = int(round(self.maxOffset/ts))
        si = tb.selectionIndices(selection)

        offset = 0
        try:
            dataArray = tb.data()[si[0]-n:si[-1]+n]
            selArray = np.pad(np.ones(si[-1] - si[0]), (n, n), 'constant', constant_values=(0, 0))
            offset = data.timeOffset(np.arange(len(dataArray), dtype='float'), dataArray, np.arange(len(selArray), dtype='float'), selArray)
        except Exception as e:
            print(f'.. problem for {selection.name}: {e}')
        
        return np.array([[si[0]+int(offset), si[-1]+int(offset)]])

    def createSettingsWidget(self):
        w = QWidget()
        w.setLayout(QFormLayout())
        w.layout().setContentsMargins(0, 0, 0, 0)
        self.maxOffsetLineEdit = QLineEdit(w)
        self.maxOffsetLineEdit.setText(str(self.maxOffset))
        self.maxOffsetLineEdit.textEdited.connect(lambda: setattr(self, 'maxOffset', int(self.maxOffsetLineEdit.text)))
        w.layout().addRow('Maximum offset (s)', self.maxOffsetLineEdit)
        return w

class RefinerWidget(QWidget):

    def __init__(self):
        super(RefinerWidget, self).__init__()
        self.setLayout(QVBoxLayout())
        self.layout().setContentsMargins(0, 0, 0, 0)
        settings = QSettings()
        uiPath = settings.value("Paths/UIPluginsPath")
        uiFile = QFile(uiPath + "/refiner.ui")
        self.ui = QUiLoader().load(uiFile)
        self.layout().addWidget(self.ui)
        self.resize(700,500)
        self.setWindowTitle("Selection Refiner")
        self.setAttribute(Qt.WA_DeleteOnClose)

        self.listWidget = self.findChild(QListWidget, 'listWidget')
        self.strategyGroupBox = self.findChild(QGroupBox, 'strategyGroupBox')
        self.strategyComboBox = self.findChild(QComboBox, 'strategyComboBox')
        self.newGroupSuffixLineEdit = self.findChild(QLineEdit, 'newGroupSuffixLineEdit')
        self.linkingGroupBox = self.findChild(QGroupBox, 'linkingGroupBox')
        self.buttonBox = self.findChild(QDialogButtonBox, 'buttonBox')
        self.progressBar = self.findChild(QProgressBar, 'progressBar')

        self.progressBar.setVisible(False)
        self.listWidget.addItems(data.selectionGroupNames())
        self.linkingGroupBox.setVisible(False)

        self.buttonBox.button(QDialogButtonBox.Close).clicked.connect(lambda: self.deleteLater())
        self.buttonBox.button(QDialogButtonBox.Apply).clicked.connect(self.apply)

        strategies = [CriteriaStrategy, MeanShiftStrategy, RollingSDStrategy, RollingSEStrategy, ExtremeStrategy, SyncStrategy]

        for strategy in strategies:
            self.strategyComboBox.addItem(strategy.name, strategy())

        self.strategyComboBox.activated.connect(self.updateSettings)
        self.strategySettings = None
        self.updateSettings()

    def updateSettings(self):
        if self.strategyGroupBox.layout().count() > 2:
            oldWidget = self.strategyGroupBox.layout().takeAt(2)
            oldWidget.widget().setParent(None)
        self.strategyGroupBox.layout().addRow(self.strategyComboBox.currentData.settings())

    def apply(self):
        groupNames = [si.text() for si in self.listWidget.selectedItems()]

        if not self.newGroupSuffixLineEdit.text:
            QMessageBox.error(self, 'Refine Selections', 'Please enter a group suffix before applying refinements')
            return
        
        self.progressBar.setValue(0)
        self.progressBar.setVisible(True)

        for groupName in groupNames:
            print('processing group %s'%groupName)
            self.progressBar.setFormat(f'%p% [{groupName}]')

            outputGroupName = '%s%s'%(groupName, self.newGroupSuffixLineEdit.text)
            if outputGroupName in data.selectionGroupNames():
                data.removeSelectionGroup(outputGroupName)

            outputGroup = data.createSelectionGroup(outputGroupName, data.selectionGroup(groupName).type)

            for si, s in enumerate(data.selectionGroup(groupName).selections()):
                try:
                    idx_pairs = self.strategyComboBox.currentData.process(s)
                    print(idx_pairs)
                    if len(idx_pairs) == 0:
                        continue

                    sorted_idx = np.argsort(idx_pairs[:,0] - idx_pairs[:,1])
                    idx_pairs = idx_pairs[sorted_idx]

                    if self.strategyComboBox.currentData.usesGlobalTime:
                        t = data.timeSeries('TotalBeam').time()
                    else:
                        t = data.timeSeries('TotalBeam').timeForSelection(s)

                    if len(idx_pairs) > 0:
                        # TODO: handle linking up to N index pairs
                        for r in [0]:
                            startTime = QDateTime.fromMSecsSinceEpoch(t[idx_pairs[r,0]]*1000.0)
                            endTime = QDateTime.fromMSecsSinceEpoch(t[idx_pairs[r,1]-1]*1000.0)
                            print('%s - %s'%(startTime, endTime))
                            data.createSelection(outputGroup, startTime, endTime, s.name + ' refined')
                except Exception as e:
                    if QMessageBox.question(self, 'Problem processing selection', f'There was a problem processing {s.name}.\n\n{e}\n\nWould you like to try to continue?', QMessageBox.Yes | QMessageBox.No) == QMessageBox.No:
                        break

                self.progressBar.setValue(100.*si/len(data.selectionGroup(groupName).selections()))
        
        self.progressBar.setVisible(False)


def main():
    global rw
    rw = RefinerWidget()
    rw.show()
