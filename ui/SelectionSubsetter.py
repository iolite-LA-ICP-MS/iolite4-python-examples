# Selection Subsetter
#/ Type: UI
#/ Name: Selection Subsetter
#/ Description: A tool for subsetting selections
#/ Authors: Bence Paul and Joe Petrus
#/ References:
#/ Version: 0.1
#/ Contact: support@iolite-software.com

from iolite.QtGui import QWidget, QVBoxLayout, QLabel, QToolButton, QAction, QComboBox, QSpinBox
from iolite.QtGui import QHBoxLayout, QDoubleSpinBox, QPushButton, QSpacerItem, QSizePolicy
from iolite.QtGui import QTabWidget, QFormLayout, QMessageBox, QKeySequence
from iolite.QtCore import Qt, Signal, QSize, QDateTime

from iolite.ui import CommonUIPyInterface as CUI
from iolite.ui import IolitePlotPyInterface as Plot

widget = None

class NumOptionsWidget(QWidget):
    numSubSelsChanged = Signal(int)

    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        self.formLayout = QFormLayout()
        self.setLayout(self.formLayout)
        self.layout().setRowWrapPolicy(QFormLayout.DontWrapRows)
        self.layout().setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)
        self.layout().setFormAlignment(Qt.AlignHCenter | Qt.AlignVCenter)

        self.numSubSelsSpinBox = QSpinBox(self)
        self.numSubSelsSpinBox.setMinimum(1)
        self.numSubSelsSpinBox.setMaximum(100000)
        self.numSubSelsSpinBox.setValue(2)
        self.numSubSelsSpinBox.setSingleStep(1)
        self.numSubSelsSpinBox.valueChanged.connect(
            lambda: self.numSubSelsChanged.emit(self.numSubSelsSpinBox.value)
        )

        self.formLayout.addRow('Number of sub-selections:', self.numSubSelsSpinBox)


class TimeOptionsWidget(QWidget):
    timeOptionsChanged = Signal(float, float) 

    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        self.formLayout = QFormLayout()
        self.setLayout(self.formLayout)
        self.layout().setRowWrapPolicy(QFormLayout.DontWrapRows)
        self.layout().setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)
        self.layout().setFormAlignment(Qt.AlignHCenter | Qt.AlignVCenter)

        self.subSelDurationDoubleSpinBox = QDoubleSpinBox(self)
        self.subSelDurationDoubleSpinBox.setDecimals(6)
        self.subSelDurationDoubleSpinBox.setSingleStep(0.001)
        self.subSelDurationDoubleSpinBox.setValue(0.1)
        self.subSelDurationDoubleSpinBox.valueChanged.connect(
            lambda: self.timeOptionsChanged.emit(self.subSelDurationDoubleSpinBox.value,
                                                 self.subSelGapDoubleSpinBox.value)
        )

        self.subSelGapDoubleSpinBox = QDoubleSpinBox(self)
        self.subSelGapDoubleSpinBox.setDecimals(6)
        self.subSelGapDoubleSpinBox.setSingleStep(0.001)
        self.subSelGapDoubleSpinBox.setValue(0.0)
        self.subSelGapDoubleSpinBox.valueChanged.connect(
            lambda: self.timeOptionsChanged.emit(self.subSelDurationDoubleSpinBox.value,
                                                 self.subSelGapDoubleSpinBox.value)
        )

        self.formLayout.addRow('Sub-selection duration (s):', self.subSelDurationDoubleSpinBox)
        self.formLayout.addRow('Sub-selection gap (s):', self.subSelGapDoubleSpinBox)

        self.setLayout(self.formLayout)


class SubsetWidget(QWidget):

    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        self.setLayout(QVBoxLayout())

        # Setting defaults here
        self.templateGroup = data.selectionGroupList(data.ReferenceMaterial | data.Sample)[0]
        self.destGroupName = 'Subset group'
        self.subSetMethod = 0
        self.numSubSels = 2
        self.subSelDuration = 0.1
        self.subSelGap = 0.0

        self.groupFormLayout = QFormLayout()
        self.templateGroupComboBox = QComboBox(self)
        self.templateGroupComboBox.addItems(data.selectionGroupNames())
        self.templateGroupComboBox.setCurrentText(self.templateGroup.name)
        self.templateGroupComboBox.currentTextChanged.connect(self.update_template_group)
        self.destGroupComboBox = QComboBox(self)
        self.destGroupComboBox.addItems(['Subset group'])
        self.destGroupComboBox.setEditable(True)
        self.destGroupComboBox.currentTextChanged.connect(self.update_dest_group)

        self.groupFormLayout.addRow('Template group:', self.templateGroupComboBox)
        self.groupFormLayout.addRow('Destination group:', self.destGroupComboBox)

        self.layout().addLayout(self.groupFormLayout)

        self.numOptionsWidget = NumOptionsWidget(self)
        self.numOptionsWidget.numSubSelsChanged.connect(self.update_num_sub_sels)
        self.timeOptionsWidget = TimeOptionsWidget(self)
        self.timeOptionsWidget.timeOptionsChanged.connect(self.update_time_options)
        self.optionsTabWidget = QTabWidget(self)
        self.optionsTabWidget.addTab(self.numOptionsWidget, 'By Number')
        self.optionsTabWidget.addTab(self.timeOptionsWidget, 'By Time')
        self.optionsTabWidget.currentChanged.connect(self.update_method)

        self.layout().addWidget(self.optionsTabWidget)

        self.processButtonLayout = QHBoxLayout()
        self.processButtonLayout.addSpacerItem(QSpacerItem(10, 10, QSizePolicy.Expanding))
        self.processButton = QToolButton(self)
        self.processButton.setText('Process')
        self.processButton.setFixedSize(100, 40)
        self.processButton.setIcon(CUI().icon('cog'))
        self.processButton.setIconSize(QSize(30,30))
        self.processButton.setToolTip('Process the selections. (Ctrl/Cmd+Return)')
        self.processButton.setShortcut(QKeySequence('Ctrl+Return'))
        self.processButton.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        self.processButton.clicked.connect(self.process)
        self.processButtonLayout.addWidget(self.processButton)

        self.layout().addLayout(self.processButtonLayout)
        
        self.buttonsLayout = QHBoxLayout()
        self.buttonsLayout.addSpacerItem(QSpacerItem(10, 10, QSizePolicy.Expanding))
        self.cancelButton = QPushButton(self)
        self.cancelButton.setText('Cancel')
        self.cancelButton.clicked.connect(lambda: self.update_or_cancel_clicked('Cancel'))
        self.buttonsLayout.addWidget(self.cancelButton)
        self.OKbutton = QPushButton(self)
        self.OKbutton.setText('OK')
        self.OKbutton.clicked.connect(lambda: self.update_or_cancel_clicked('OK'))
        self.buttonsLayout.addWidget(self.OKbutton)

        self.layout().addLayout(self.buttonsLayout)

    def update_template_group(self):
        self.templateGroup = data.selectionGroup(self.templateGroupComboBox.currentText)

    def update_dest_group(self):
        self.destGroupName = self.destGroupComboBox.currentText
    
    def update_method(self):
        self.subSetMethod = self.optionsTabWidget.currentIndex

    def update_num_sub_sels(self):
        self.numSubSels = self.numOptionsWidget.numSubSelsSpinBox.value

    def update_time_options(self):
        self.subSelDuration = self.timeOptionsWidget.subSelDurationDoubleSpinBox.value
        self.subSelGap = self.timeOptionsWidget.subSelGapDoubleSpinBox.value

    def process(self):
        print('Process called....')

        if(self.destGroupName == ''):
            QMessageBox.warning(self, 'Invalid Destination Group', 'Please enter a valid destination group name.')
            return

        if(self.destGroupName in data.selectionGroupNames()):
            reply = QMessageBox.question(
                self,
                'Overwrite Selections',
                'This group already exists. All selections within the group will be overwritten. Do you want to proceed?',
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No
            )
            if reply == QMessageBox.No:
                return
            
            if reply == QMessageBox.Yes:
                dest_grp = data.selectionGroup(self.destGroupName)
                dest_grp.removeSelections(dest_grp.selections())
        else:
            data.createSelectionGroup(self.destGroupName, self.templateGroup.type)
            dest_grp = data.selectionGroup(self.destGroupName)

        s_list = [] # List of start times
        e_list = [] # List of end times
        l_list = [] # List of labels

        if self.subSetMethod == 0:
            print('By Number')
            for sel in self.templateGroup.selections():
                start = sel.startTime
                total_duration = sel.startTime.msecsTo(sel.endTime)
                sub_duration = total_duration / self.numSubSels
                for i in range(self.numSubSels):
                    this_start = start.addMSecs(i*sub_duration)
                    s_list.append(this_start)
                    e_list.append(this_start.addMSecs(sub_duration))
                    l_list.append(sel.name + '_sub' + str(i))
                    i += 1

        elif self.subSetMethod == 1:
            print('By Time')
            for sel in self.templateGroup.selections():
                start = sel.startTime
                end = sel.endTime
                time_d = float(sel.startTime.toMSecsSinceEpoch()) #NOTE: converting to float to give us better than millisecond precision
                i = 1
                while time_d < end.toMSecsSinceEpoch():
                    s_list.append(QDateTime().fromMSecsSinceEpoch(time_d))
                    e_list.append(QDateTime().fromMSecsSinceEpoch(time_d + self.subSelDuration*1000.))
                    time_d += self.subSelDuration*1000. + self.subSelGap*1000.
                    l_list.append(sel.name + '_sub' + str(i))
                    i += 1
            
        print(f's_list: {s_list}')
        print(f'e_list: {e_list}')
        print(f'l_list: {l_list}')

        data.createSelections(dest_grp, s_list, e_list, l_list)

        #data.setActiveSelectionGroup(self.destGroupName)
        
    
    def update_or_cancel_clicked(self, button: str):
        print('Update or cancel clicked')
        print(button)
        if(button == 'Cancel'):
            if(self.destGroupName in data.selectionGroupNames()):
                data.removeSelectionGroup(self.destGroupName)
            self.close()
        else:
            self.close()


def create_widget():
    global widget
    widget = SubsetWidget()
    widget.setAttribute(Qt.WA_DeleteOnClose)
    widget.setWindowFlags(Qt.WindowStaysOnTopHint)

    widget.show()


def createUIElements():
    action = QAction("Selection Subsetter", ui)
    action.triggered.connect(create_widget)
    ui.setAction(action)
    ui.setMenuName(['Selections'])
