#/ Type: UI
#/ Name: ChannelPruner
#/ Description: 
#/ Authors: J. Petrus & B. Paul
#/ References:
#/ Version: 1.1
#/ Contact:

from iolite.QtGui import QAction, QDialog, QDialogButtonBox, QCheckBox, QFormLayout, QMessageBox
from math import log10
import numpy as np
import re
widget = None

def createUIElements():
    action = QAction('Channel Pruner', None)
    action.triggered.connect(create)
    ui.setAction(action)
    ui.setMenuName(['Tools'])

def create():
    global widget
    widget = Pruner()

    if widget.exec_() == widget.Rejected:
        return
    
    widget.process()
    widget.deleteLater()

def nullOrEmptyString(s):
    return s is None or type(s) is not str or len(s) == 0

class Pruner(QDialog):
    
    def __init__(self, parent=None):
        QDialog.__init__(self, parent)

        self.setLayout(QFormLayout())
        self.noElementCB = QCheckBox(self)
        self.noMassCB = QCheckBox(self)
        self.lowSignalCB = QCheckBox(self)
        self.oddNameCB = QCheckBox(self)
        self.setBLsubCB = QCheckBox(self)
        self.verboseCB = QCheckBox(self)

        self.layout().addRow('No element property', self.noElementCB)
        self.layout().addRow('No mass property', self.noMassCB)
        self.layout().addRow('Low intensity', self.lowSignalCB)
        self.layout().addRow('Odd name', self.oddNameCB)
        self.layout().addRow('Set BLsub property', self.setBLsubCB)
        self.layout().addRow('Verbose?', self.verboseCB)

        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, self)
        self.layout().addRow(bb)
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)


    def process(self):
        namesToRemove = []

        inputChannels = [c for c in data.timeSeriesList(data.Input) if 'TotalBeam' not in c.name]

        for channel in inputChannels:
            if self.noElementCB.isChecked() and nullOrEmptyString(channel.property('Element')):
                namesToRemove.append(channel.name)
                if self.verboseCB.isChecked():
                    print(f'Suggesting removal of {channel.name} due to no element set')

            if self.noMassCB.isChecked() and nullOrEmptyString(channel.property('Mass')):
                namesToRemove.append(channel.name)
                if self.verboseCB.isChecked():
                    print(f'Suggesting removal of {channel.name} due to no mass set')

            if self.lowSignalCB.isChecked() and log10(np.nanmax(channel.data()) - np.nanmin(channel.data())) < 0.2:
                namesToRemove.append(channel.name)
                if self.verboseCB.isChecked():
                    print(f'Suggesting removal of {channel.name} due to less than approx 1.5 CPS difference between channel max and min')

            if self.oddNameCB.isChecked() and re.match('^[A-Za-z]{1,2}\d{1,3}$', channel.name) is None:
                namesToRemove.append(channel.name)
                if self.verboseCB.isChecked():
                    print(f'Suggesting removal of {channel.name} due to a name that doesn\'t match the [Symbol][Mass] format')

        if QMessageBox.question(self, 'Remove channels', f'Are you sure you want to remove:\n{",".join(namesToRemove)}', QMessageBox.Yes | QMessageBox.No) == QMessageBox.No:
            return

        for name in namesToRemove:
            print(f'Removing {name}')
            try:
                data.removeTimeSeries(name)
            except:
                pass

        if self.setBLsubCB.isChecked():
            print('Setting BackgroundSubtracted property...')
            for ch in data.timeSeriesList(data.Input):
                ch.setProperty("BackgroundSubtracted", True)
