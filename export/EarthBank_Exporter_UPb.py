#/ Type: Exporter
#/ Name: EarthBank U-Pb Exporter
#/ Authors: Bence Paul and Joe Petrus
#/ Description: An export script for creating an EarthBank formatted spreadsheet for importing into the EarthBank database.
#/ References: None
#/ Version: 0.1
#/ Contact: iolitesupport@icpms.com

import os
import subprocess
import sys
from shutil import copy2
from pprint import pprint
import numpy as np

# We'll use the weighted mean helper function from iolite_helpers for calculations.
from iolite_helpers import weightedMean

from iolite.ui import CommonUIPyInterface as CUI

# iolite may execute scripts with an unreliable __file__. The code object's
# filename is also used by tracebacks and normally contains the script path.
SCRIPT_PATH = next(
    (path for path in (sys._getframe().f_code.co_filename,
                      globals().get('__file__'))
     if path and os.path.isfile(path)),
    None,
)
if SCRIPT_PATH is None:
    raise RuntimeError(
        'Could not determine the exporter script location from iolite. '
        'Run the exporter using its saved Python file.'
    )
EXPORTER_DIR = os.path.dirname(os.path.abspath(SCRIPT_PATH))

if EXPORTER_DIR not in sys.path:
    sys.path.insert(0, EXPORTER_DIR)

try:
    from exporter_config import ConfigField, ExporterConfig
except ModuleNotFoundError as exc:
    if exc.name != 'exporter_config':
        raise
    raise ModuleNotFoundError(
        f'Could not find exporter_config.py in {EXPORTER_DIR}. '
        'Copy exporter_config.py from the iolite4-python-examples export '
        'directory into this folder, then rerun the exporter.'
    ) from exc

# UI file should be in the same directory as the exporter script.
UI_FILE_PATH = os.path.join(EXPORTER_DIR, 'EarthBank_Exporter_UPb.ui')
# Check if the UI file exists.
if not os.path.isfile(UI_FILE_PATH):
    raise FileNotFoundError(
        f'Could not find the UI file at {UI_FILE_PATH}. '
        'Ensure that EarthBank_Exporter_UPb.ui is in the same directory as the exporter script.'
    )

# Template file should be in the same directory as the exporter script.
TEMPLATE_PATH = os.path.join(EXPORTER_DIR, 'UPbDatapoint_LAICPMS.template.v2026-06-22.xlsx')
# Check if the template file exists.
if not os.path.isfile(TEMPLATE_PATH):
    raise FileNotFoundError(
        f'Could not find the template file at {TEMPLATE_PATH}. '
        'Ensure that UPbDatapoint_LAICPMS.template.v2026-06-22.xlsx is in the same directory as the exporter script.'
    )

from openpyxl import load_workbook

from iolite.Qt import Qt
from iolite.QtCore import QFile, QIODevice, QEventLoop, QUrl
from iolite.QtGui import QComboBox, QHeaderView, QLabel, QLineEdit, QMessageBox, \
    QStackedWidget, QTableWidget, QTableWidgetItem, QToolButton, QVBoxLayout, QWidget, \
    QDoubleSpinBox, QMenu, QDesktopServices
from iolite.QtGui import QCheckBox
from iolite.QtUiTools import QUiLoader

'''
Default configuration for the EarthBank U-Pb exporter.
'''
USERS = ['Bence Paul', 'Alan Grieg', 'Ashlea Wainwright', 'Brandon Mahan', 'Janet Hergt', 'Jon Woodhead', 'Roland Maas']
LAB = 'Melbourne Isotope Analytics,  The University of Melbourne'

LASER_WAVELENGTHS = ['193 nm', '213 nm', '266 nm']
LASER_PULSE_WIDTHS = ['< 5 ns', '5-10 ns', '> 10 ns']
CELL_MODEL = ['TV2', 'TV3', 'Helex']
CARRIER_GAS = ['He', 'Ar', 'N2']
MIXING_DEVICES = ['Squid', 'ESL Mixing bulb', 'None', 'Other']
BATCH_ID_OPTIONS = ['Session UUID', 'Session file path']

# This is the key used to store the config file's path in the QSettings.
# Only the path to the config file is stored, not the contents.
CONFIG_PATH_KEY = 'EarthBank_Exporter_UPb/config_path'

DEFAULT_CONFIG = {
    'users': USERS,
    'lab': LAB,
    'batch_id': BATCH_ID_OPTIONS,
    'default_batch_id': 'Session UUID',
    'laser_wavelengths': LASER_WAVELENGTHS,
    'laser_pulse_widths': LASER_PULSE_WIDTHS,
    'cell_models': CELL_MODEL,
    'carrier_gases': CARRIER_GAS,
    'mixing_devices': MIXING_DEVICES,
    'default_user': USERS[0],
    'default_laser_wavelength': '193 nm',
    'default_laser_pulse_width': '< 5 ns',
    'default_cell_model': 'TV2',
    'default_carrier_gas': 'He',
    'default_mixing_device': 'Squid',
    'laser_fluence': 3.0,
    'carrier_gas_flow_rate': 1.0,
}

CONFIG_CHOICES = (
    ('users', 'default_user', 'userComboBox'),
    ('batch_id', 'default_batch_id', 'batchIDComboBox'),
    ('laser_wavelengths', 'default_laser_wavelength', 'laserWavelengthComboBox'),
    ('laser_pulse_widths', 'default_laser_pulse_width', 'laserPulseWidthComboBox'),
    ('cell_models', 'default_cell_model', 'cellModelComboBox'),
    ('carrier_gases', 'default_carrier_gas', 'carrierGasComboBox'),
    ('mixing_devices', 'default_mixing_device', 'mixingDeviceComboBox'),
)


CONFIG_FIELDS = {
    **{key: ConfigField(type(value)) for key, value in DEFAULT_CONFIG.items()},
}

CONFIG_MANAGER = ExporterConfig(
    defaults=DEFAULT_CONFIG,
    fields=CONFIG_FIELDS,
    settings_key=CONFIG_PATH_KEY,
    title='EarthBank U-Pb Exporter',
    choices={selected: options for options, selected, _ in CONFIG_CHOICES},
)


class SettingsWidget(QWidget):

    def __init__(self, config_path, config, parent=None):
        super().__init__(parent)

        self.config_path = config_path
        self.config = config
        self.setAttribute(Qt.WA_DeleteOnClose, True)

        self.setLayout(QVBoxLayout())

        # Load the template file for the AGN export. This is an Excel file that contains the correct formatting for the AGN database.
        self.export_template = QFile(TEMPLATE_PATH)
        if not self.export_template.exists():
            raise RuntimeError(f'Could not find export template {self.export_template.fileName()}')

        if not self.export_template.open(QIODevice.ReadOnly):
            raise RuntimeError(f'Could not load export template file for reading: {self.export_template.fileName()}')

        # Load the UI file
        self.ui_file = QFile(UI_FILE_PATH)

        if not self.ui_file.exists():
            raise RuntimeError(f'Could not find settings ui {self.ui_file.fileName()}')

        if not self.ui_file.open(QIODevice.ReadOnly):
            raise RuntimeError(f'Could not load settings ui file for reading: {self.ui_file.fileName()}')

        ui = QUiLoader().load(self.ui_file, self)
        self.layout().addWidget(ui)
        self.layout().setContentsMargins(0, 0, 0, 0)

        self.stackedWidget = ui.findChild(QStackedWidget, 'stackedWidget')
        self.cancelButton = ui.findChild(QToolButton, 'cancelButton')
        self.continueButton = ui.findChild(QToolButton, 'continueButton')
        self.backButton = ui.findChild(QToolButton, 'backButton')
        self.stepLabel = ui.findChild(QLabel, 'stepLabel')
        self.configsButton = ui.findChild(QToolButton, 'configsButton')

        self.userComboBox = ui.findChild(QComboBox, 'userComboBox')
        self.labLineEdit = ui.findChild(QLineEdit, 'labLineEdit')
        self.litLineEdit = ui.findChild(QLineEdit, 'litLineEdit')
        self.fundingLineEdit = ui.findChild(QLineEdit, 'fundingLineEdit')
        self.techniqueComboBox = ui.findChild(QComboBox, 'techniqueComboBox')
        self.batchIDComboBox = ui.findChild(QComboBox, 'batchIDComboBox')
        self.uncertTypeComboBox = ui.findChild(QComboBox, 'uncertaintyTypeComboBox')
        self.trimSigFigsCheckBox = ui.findChild(QCheckBox, 'trimSigFigsCheckBox')
        self.sessionIDComboBox = ui.findChild(QComboBox, 'sessionIDComboBox')
        self.groupsTable = ui.findChild(QTableWidget, 'groupsTableWidget')
        self.duplicatesLabel = ui.findChild(QLabel, 'ratiosLabel')
        self.ratiosTable = ui.findChild(QTableWidget, 'ratiosTableWidget')
        self.laserFluenceDoubleSpinBox = ui.findChild(QDoubleSpinBox, 'laserFluenceDoubleSpinBox')
        self.laserWavelengthComboBox = ui.findChild(QComboBox, 'laserWavelengthComboBox')
        self.laserPulseWidthComboBox = ui.findChild(QComboBox, 'laserPulseWidthComboBox')
        self.cellModelComboBox = ui.findChild(QComboBox, 'cellModelComboBox')
        self.carrierGasComboBox = ui.findChild(QComboBox, 'carrierGasComboBox')
        self.carrierGasFlowRateSpinBox = ui.findChild(QDoubleSpinBox, 'carrierGasFlowRateDoubleSpinBox')
        self.mixingDeviceComboBox = ui.findChild(QComboBox, 'mixingDeviceComboBox')

        # ConfigsButton menu setup
        self.configsMenu = QMenu(self.configsButton)
        self.configsMenu.addAction('Reveal Config File', self.revealConfigFile)
        # Add option to save current configuration
        self.configsMenu.addAction('Save Config File', self.saveDefaults)
        self.configsMenu.addAction('Select Config File', self.changeConfig)
        self.configsButton.setMenu(self.configsMenu)
        self.configsButton.setIcon(CUI().icon('book'))
        self.configsButton.setToolTip('Configuration Options')
        self.configsButton.setPopupMode(QToolButton.InstantPopup)

        self.addResultButton = ui.findChild(QToolButton, 'addResultButton')
        self.removeResultButton = ui.findChild(QToolButton, 'removeResultButton')

        # Get the valid options from the template file
        template_wb = load_workbook(TEMPLATE_PATH)
        lookup_tables_ws = template_wb['Lookup Tables']
        lookup_values = {}
        for column in range(1, lookup_tables_ws.max_column + 1):
            lookup_name = lookup_tables_ws.cell(row=2, column=column).value
            if lookup_name in (None, 'Description'):
                continue

            values = []
            for row in range(3, lookup_tables_ws.max_row + 1):
                value = lookup_tables_ws.cell(row=row, column=column).value
                if value in (None, ''):
                    break
                values.append(value)
            lookup_values[lookup_name] = values

        # Some sanity checks to make sure that the template file contains the expected values
        if 'U-Pb Analytical Technique' not in lookup_values:
            raise RuntimeError('Could not find "U-Pb Analytical Technique" in the template file. Please check that the template file is correct.')

        if 'Mineral Type' not in lookup_values:
            raise RuntimeError('Could not find "Mineral Type" in the template file. Please check that the template file is correct.')

        self.batchIDComboBox.addItems(['Session UUID', 'Session file path'])

        # Also get all the column indices for the relevant sheets in the template file
        try:
            dps_ws = template_wb['UPb Datapoints']
            upb_ws = template_wb['UPbSpotData']
            icpms_ws = template_wb['ICPMS']
        except KeyError:
            print('Could not find expected sheet in template file')
            return

        # Get the column indices for the various data points in the template file
        self.dps_col_indices = {}
        for col in dps_ws.iter_cols(min_row=3, max_row=3):
            header = col[0].value
            if header:
                self.dps_col_indices[header] = col[0].column

        # Sanity check to make sure that we actually found some column indices
        if not self.dps_col_indices:
            raise RuntimeError('Could not find any column indices in the "UPb Datapoints" sheet of the template file.')

        # Repeat for UPbSpotData sheet
        self.upb_col_indices = {}
        for col in upb_ws.iter_cols(min_row=4, max_row=4):
            header = str(col[0].value)
            if header:
                # Skip the uncertainty, type and method columns
                if "uncertainty" in header.lower() or "type" in header.lower() or "method" in header.lower():
                    continue

                self.upb_col_indices[header] = col[0].column

        if not self.upb_col_indices:
            raise RuntimeError('Could not find any column indices in the "UPbSpotData" sheet of the template file.')

        # Repeat for ICP-MS sheet
        self.icpms_col_indices = {}
        for col in icpms_ws.iter_cols(min_row=4, max_row=4):
            header = col[0].value
            if header:
                self.icpms_col_indices[header] = col[0].column

        if not self.icpms_col_indices:
            raise RuntimeError('Could not find any column indices in the "ICPMS" sheet of the template file.')

        # Fill mineral type combo box
        self.mineralTypeComboBox = ui.findChild(QComboBox, 'mineralComboBox')
        self.mineralTypeComboBox.addItems(lookup_values['Mineral Type'])
        self.mineralTypeComboBox.setCurrentText('Zircon')
        
        self.techniqueComboBox.addItems(lookup_values['U-Pb Analytical Technique'])
        self.techniqueComboBox.setCurrentText('LA-ICP-MS')
        
        self.uncertTypeComboBox.addItems(['2 standard error', '1σ', '2σ', '95%'])

        self.applyConfig()

        '''
        Get primary RMs (calibrants) as these should not be included in exported results.
        TODO: Implement getting primary RM from current DRS settings NOTE: This could be 
        as simple as adding a property to the intermediate and output channels from the DRS.
        '''
        primaryRM = ''

        for ch in data.timeSeriesList(data.Intermediate):
            extStd = ch.property('Primary RM')
            if extStd is not None:
                primaryRM = extStd

        # Set up the groups table
        self.groupsTable.setColumnCount(6)
        self.groupsTable.setHorizontalHeaderLabels(['Group Name', 'Type', 'Selection Count', 'Weighted Mean age', 'Treat as Age Group?', 'Mount ID'])
        self.groupsTable.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.groupsTable.setEditTriggers(QTableWidget.DoubleClicked | QTableWidget.EditKeyPressed)

        def createGroupTableEntry(group, type):
            row = self.groupsTable.rowCount
            self.groupsTable.insertRow(row)

            groupNameItem = QTableWidgetItem(group.name)
            if group.name is primaryRM:
                groupNameItem.setCheckState(Qt.Unchecked)
            else:
                groupNameItem.setCheckState(Qt.Checked)
            self.groupsTable.setItem(row, 0, groupNameItem)

            groupTypeItem = QTableWidgetItem(type)
            self.groupsTable.setItem(row, 1, groupTypeItem)

            selCountItem = QTableWidgetItem(str(group.count))
            selCountItem.setTextAlignment(Qt.AlignHCenter)
            selCountItem.setFlags(selCountItem.flags() & ~Qt.ItemIsEditable)
            self.groupsTable.setItem(row, 2, selCountItem)

            # Get group 6/38 age and uncertainty
            grp_w_mean, grp_w_unc = self.calc_group_6_38_age(group)

            groupMeanItem = QTableWidgetItem(f'{grp_w_mean} ± {grp_w_unc}')
            groupMeanItem.setTextAlignment(Qt.AlignHCenter)
            groupMeanItem.setFlags(groupMeanItem.flags() & ~Qt.ItemIsEditable)
            self.groupsTable.setItem(row, 3, groupMeanItem)

            treatAsAgeGroupItem = QTableWidgetItem()
            treatAsAgeGroupItem.setFlags(treatAsAgeGroupItem.flags() | Qt.ItemIsUserCheckable)
            treatAsAgeGroupItem.setCheckState(Qt.Unchecked)
            self.groupsTable.setItem(row, 4, treatAsAgeGroupItem)

            mount_ID_item = QTableWidgetItem()
            mount_ID_item.setText(f'Mount ID {row + 1}')
            mount_ID_item.setTextAlignment(Qt.AlignHCenter)
            mount_ID_item.setFlags(mount_ID_item.flags() | Qt.ItemIsEditable)
            self.groupsTable.setItem(row, 5, mount_ID_item)


        for group in data.selectionGroupList(data.ReferenceMaterial):
            createGroupTableEntry(group, 'Reference Material')

        for group in data.selectionGroupList(data.Sample):
            createGroupTableEntry(group, 'Sample')

        # All this just to properly size the columns based on their content and header text
        self.groupsTable.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        self.groupsTable.resizeColumnsToContents()
        header_font_metrics = self.groupsTable.horizontalHeader().fontMetrics()
        for column in range(self.groupsTable.columnCount):
            header_text = self.groupsTable.horizontalHeaderItem(column).text()
            header_width = header_font_metrics.horizontalAdvance(header_text) + 30
            header_width = max(header_width, self.groupsTable.horizontalHeader().sectionSizeHint(column))
            self.groupsTable.setColumnWidth(column, max(self.groupsTable.columnWidth(column), header_width))

        # Set up the ratios table
        self.ratiosTable.setColumnCount(2)
        self.ratiosTable.setHorizontalHeaderLabels(['Result', 'Channel'])
        self.ratiosTable.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.ratiosTable.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        self.ratiosTable.verticalHeader().hide() 

        self.addResultButton.clicked.connect(self.addResultRow)
        self.removeResultButton.clicked.connect(self.removeResultRow)

        U_Pb_field_defaults = {
            'U concentration (µg.g-1)': 'Approx_U_PPM',
            'Th concentration (µg.g-1)': 'Approx_Th_PPM',
            '206Pb intensity (mV or cps)': 'Pb206_CPS',
            '204Pb intensity (mV or cps)': 'Pb204_CPS',
            'Total Pb concentration (µg.g-1)': 'Approx_Pb_PPM',
            '206Pb/238U Ratio': 'Final Pb206/U238',
            '238U/206Pb Ratio': 'Final U238/Pb206',
            '207Pb/206Pb Ratio': 'Final Pb207/Pb206',
            '207Pb/235U Ratio': 'Final Pb207/U235',
            '208Pb/232Th Ratio': 'Final Pb208/Th232',
            '206Pb/238U Date (Ma)': 'Final Pb206/U238 age',
            '207Pb/206Pb Date (Ma)': 'Final Pb207/Pb206 age', 
            '207Pb/235U Date (Ma)': 'Final Pb207/U235 age', 
            '208Pb/232Th Date (Ma)': 'Final Pb208/Th232 age',
        }

        channels = [ch.name for ch in data.timeSeriesList()]
        channels.insert(0, 'None')

        for option, default_channel in U_Pb_field_defaults.items():
            row = self.ratiosTable.rowCount
            self.ratiosTable.insertRow(row)
            self.ratiosTable.setItem(row, 0, QTableWidgetItem(option))
            chItem = QComboBox()
            chItem.addItems(channels)
            chItem.setCurrentText(default_channel if default_channel in channels else 'None')
            self.ratiosTable.setCellWidget(row, 1, chItem)

        self.backButton.clicked.connect(self.previousTab)
        self.continueButton.clicked.connect(self.nextTab)
        self.cancelButton.clicked.connect(lambda: self.close())

        self.stackedWidget.setCurrentIndex(0)
        self.stepLabel.setText(f'Step {self.stackedWidget.currentIndex+1} of 4')

    def applyConfig(self):
        for options_key, selected_key, widget_name in CONFIG_CHOICES:
            combo = getattr(self, widget_name)
            combo.clear()
            combo.addItems(self.config[options_key])
            combo.setCurrentText(self.config[selected_key])
        self.labLineEdit.setText(self.config['lab'])
        self.laserFluenceDoubleSpinBox.setValue(self.config['laser_fluence'])
        self.carrierGasFlowRateSpinBox.setValue(self.config['carrier_gas_flow_rate'])

    def changeConfig(self):
        result = CONFIG_MANAGER.choose(self, use_remembered=False)
        if result is not None:
            self.config_path, self.config = result
            self.applyConfig()

    def saveDefaults(self):
        updates = {'lab': self.labLineEdit.text}
        for options_key, selected_key, widget_name in CONFIG_CHOICES:
            updates[options_key] = self.config[options_key]
            updates[selected_key] = getattr(self, widget_name).currentText
        updates['laser_fluence'] = self.laserFluenceDoubleSpinBox.value
        updates['carrier_gas_flow_rate'] = self.carrierGasFlowRateSpinBox.value
        config = CONFIG_MANAGER.save_updates(self.config_path, updates, self)
        if config is not None:
            self.config = config

    def revealConfigFile(self):
        if self.config_path and os.path.exists(self.config_path):
            QDesktopServices.openUrl(QUrl.fromLocalFile(self.config_path))
        else:
            QMessageBox.warning(self, 'Config file not found', 'The configuration file could not be found.')

    def previousTab(self):
        current_index = self.stackedWidget.currentIndex
        self.stackedWidget.setCurrentIndex(current_index-1)
        self.stepLabel.setText(f'Step {self.stackedWidget.currentIndex+1} of 4')
        self.continueButton.setText("Continue")

    def nextTab(self):
        current_index = self.stackedWidget.currentIndex
        # Check that at least one group is selected in groupsTable
        if current_index == 2 and all(self.groupsTable.item(row, 0).checkState() != Qt.Checked for row in range(self.groupsTable.rowCount)):
            QMessageBox.warning(self, 'No groups selected', 'Please select at least one group before continuing.')
            return
        # Check that not all groups are selected in groupsTable (the user should leave at least one group unselected... most likely)
        if current_index == 2 and all(self.groupsTable.item(row, 0).checkState() == Qt.Checked for row in range(self.groupsTable.rowCount)):
            response = QMessageBox.question(
                self,
                'All groups selected',
                'All groups are selected. Click OK to continue, or Cancel to deselect one of the groups.',
                QMessageBox.Ok | QMessageBox.Cancel,
                QMessageBox.Cancel,
            )
            if response == QMessageBox.Cancel:
                return

        # If we're already on the last tab, export the data
        if current_index == self.stackedWidget.count - 1:
            self.exportData()
            print('Export complete.')
            return
        self.stackedWidget.setCurrentIndex(current_index+1)
        self.stepLabel.setText(f'Step {self.stackedWidget.currentIndex+1} of 4')
        if self.stackedWidget.currentIndex == self.stackedWidget.count - 1:
            self.continueButton.setText("Export Data")
        else:
            self.continueButton.setText("Continue")


    def calc_group_6_38_age(self, group):
        # Make sure to pass the group OBJECT not name string

        # Get 206/238 ages and uncertainties for this group:
        age206_238_ch = data.timeSeries('Final Pb206/U238 age')
        results = [data.result(sel, age206_238_ch) for sel in group.selections()]
        ages = np.array([result.value() for result in results])
        uncertainties = np.array([result.propagatedUncertainty() for result in results])

        mean_results = weightedMean(ages, uncertainties)
        grp_w_mean, grp_w_unc = mean_results['internal']

        # Add the mean and uncertainty to the table, but trim the numbers to a reasonable number of decimal places
        # using significant figures of the uncertainty
        # Determine the number of decimal places based on the uncertainty
        if grp_w_unc > 0:
            decimals = -int(np.floor(np.log10(grp_w_unc)))
        else:
            decimals = 2
        grp_w_mean = round(grp_w_mean, decimals)
        grp_w_unc = round(grp_w_unc, decimals)

        return grp_w_mean, grp_w_unc

    def addResultRow(self):
        row_count = self.ratiosTable.rowCount
        self.ratiosTable.insertRow(row_count)
        # Add a combobox in the first column to select the type of ratio
        ratio_combo = QComboBox()
        # Set the items to the dict keys of upb_col_indices
        upb_keys = [key for key in self.upb_col_indices.keys()]
        ratio_combo.addItems(upb_keys)
        ratio_combo.setCurrentIndex(27)
        self.ratiosTable.setCellWidget(row_count, 0, ratio_combo)

        # Create a combobox with all channels available in it
        channel_combo = QComboBox()
        channel_combo.addItem('None')
        channel_combo.addItems(data.timeSeriesNames())
        channel_combo.setCurrentIndex(0)
        self.ratiosTable.setCellWidget(row_count, 1, channel_combo)

        # Move the focus to the newly added row
        self.ratiosTable.scrollToItem(self.ratiosTable.item(row_count, 0))
        self.ratiosTable.setCurrentCell(row_count, 0)

    def removeResultRow(self):
        last_row = self.ratiosTable.rowCount - 1
        if last_row > 0:
            self.ratiosTable.removeRow(last_row)
            self.ratiosTable.setCurrentCell(last_row - 1, 0)

    def exportData(self):
        print('Exporting data...')

        try:
            # NOTE: the export_filepath variable is automatically set for the script by iolite, 
            # and is the path to the file that the user has selected to export to.
            fp = copy2(TEMPLATE_PATH, export_filepath)
        except FileNotFoundError:
            print('Could not find template file')
            return
        except PermissionError:
            print('Permission denied to write to export file')
            return
        except:
            print('Could not create copy of template file')
            return

        wb = load_workbook(fp)
        dps_ws = wb['UPb Datapoints']
        upb_ws = wb['UPbSpotData']
        icpms_ws = wb['ICPMS']
        agegrp_ws = wb['UPbAgeGroup']

        # Get dict of groups to export, with their mount IDs
        groups_to_export = {}
        for row in range(self.groupsTable.rowCount):
            grpName = self.groupsTable.item(row, 0).text()
            if self.groupsTable.item(row, 0).checkState() == Qt.Checked:
                this_dict = {}
                this_dict['Age Group'] = self.groupsTable.item(row, 4).checkState() == Qt.Checked if self.groupsTable.item(row, 4) is not None else None
                this_dict['Mount ID'] = self.groupsTable.item(row, 5).text() if self.groupsTable.item(row, 5) is not None else None
                groups_to_export[grpName] = this_dict
            
        groupCounter = 5
        rowNo = 5
        ageGroupCounter = 5

        for g in groups_to_export.keys():
            print(f'Exporting group: {g}')
            group = data.selectionGroup(g)
            datapoint_name = group.name
            dps_ws.cell(groupCounter, 1, value=datapoint_name)

            if group.type == data.Sample:
                dps_ws.cell(groupCounter, 2, value=group.name)
            elif group.type == data.ReferenceMaterial:
                dps_ws.cell(groupCounter, 3, value=group.name)
            
            # Export metadata to the template file
            dps_ws.cell(groupCounter, self.dps_col_indices['U-Pb Analytical Technique'], value = self.techniqueComboBox.currentText)

            if self.batchIDComboBox.currentText == 'Session UUID':
                dps_ws.cell(groupCounter, 4, value = data.sessionUUID())
            elif self.batchIDComboBox.currentText == 'Session file path':
                dps_ws.cell(groupCounter, 4, value = data.sessionFilePath())

            dps_ws.cell(groupCounter, 5, value = data.sessionUUID())
            dps_ws.cell(groupCounter, 6, value = data.sessionFilePath())
            dps_ws.cell(groupCounter, self.dps_col_indices['Analyst'], value = self.userComboBox.currentText)
            dps_ws.cell(groupCounter, self.dps_col_indices['Laboratory'], value = self.labLineEdit.text)
            dps_ws.cell(groupCounter, self.dps_col_indices['Mineral Type'], value = self.mineralTypeComboBox.currentText)
            # Note the typo in 'Associated Litterature' which is in the template. When this gets changed, we'll need to update this line
            dps_ws.cell(groupCounter, self.dps_col_indices['Associated litterature'], value = self.litLineEdit.text)
            dps_ws.cell(groupCounter, self.dps_col_indices['Funding'], value = self.fundingLineEdit.text)
            grp_dict = groups_to_export[g]
            
            dps_ws.cell(groupCounter, self.dps_col_indices['Mount ID'], value = grp_dict['Mount ID'])
            # If we're treating this group as an age group, add to the UPbAgeGroup sheet
            if grp_dict['Age Group']:
                agegrp_ws.cell(ageGroupCounter, 1, value=datapoint_name)
                agegrp_ws.cell(ageGroupCounter, 3, value=group.name)
                agegrp_ws.cell(ageGroupCounter, 5, value=len(group.selections()))
                # TODO: Add group age data here, depending on calculations to come...
                grp_w_mean, grp_w_unc = self.calc_group_6_38_age(group)
                agegrp_ws.cell(ageGroupCounter, 7, value=grp_w_mean)
                agegrp_ws.cell(ageGroupCounter, 8, value=grp_w_unc)
                agegrp_ws.cell(ageGroupCounter, 9, value='2 SE')

                ageGroupCounter += 1


            # Now export U-Pb data
            if rowNo == 5:
                rowNo += 1 # Skip the header row for U-Pb data

            # Keep track of repetition rates, spot heights, and spot widths for each selection
            rep_rates = set()
            spot_heights = set()
            spot_widths = set()
            
            for sel in group.selections():
                # Exclude component selections (linked selections will be exported, but components will not)
                if sel.hasLinks():
                    continue

                upb_ws.cell(row=rowNo, column=1, value=datapoint_name)
                upb_ws.cell(row=rowNo, column=3, value=sel.name)
                upb_ws.cell(row=rowNo, column=4, value=sel.startTime.toString('yyyy-MM-dd HH:mm:ss'))
                upb_ws.cell(row=rowNo, column=6, value=sel.comment)

                rep_rates.add(sel.property('Rep Rate'))
                spot_heights.add(sel.property('Spot Height'))
                spot_widths.add(sel.property('Spot Width'))

                for row in range(self.ratiosTable.rowCount):
                    field_name = self.ratiosTable.item(row, 0).text() if self.ratiosTable.item(row, 0) is not None else None
                    ch_name = self.ratiosTable.cellWidget(row, 1).currentText if self.ratiosTable.cellWidget(row, 1) is not None else None

                    # The field_name might be None if the user has added a field, so check for that here:
                    if field_name is None:
                        field_name = self.ratiosTable.cellWidget(row, 0).currentText if self.ratiosTable.cellWidget(row, 0) is not None else None

                    if field_name not in self.upb_col_indices:
                        print(f'Field name {field_name} at row {row} is not in the UPb column indices')
                        continue

                    # The user may have chosen "None" as the channel name, which just means they don't want to export this value
                    if ch_name is None or ch_name.lower() == 'none':
                        continue

                    result = data.result(sel, data.timeSeries(ch_name)).value()
                    uncert = data.result(sel, data.timeSeries(ch_name)).propagatedUncertainty()

                    upb_ws.cell(row=rowNo, column=self.upb_col_indices[field_name], value=result)
                    if isinstance(uncert, float) and (uncert != uncert):  # Check for NaN
                        pass
                        #print(f'No propagated uncertainty for field {field_name} at row {row}')
                    else:
                        upb_ws.cell(row=rowNo, column=self.upb_col_indices[field_name] + 1, value=uncert)

                rowNo += 1

            # If any of the rep rate, spot height, or spot width have a single unique value, export them to the ICP-MS sheet
            if len(rep_rates) == 1 or len(spot_heights) == 1 or len(spot_widths) == 1:
                icpms_ws.cell(row=groupCounter, column=self.icpms_col_indices['datapointName'], value=datapoint_name)
            if len(rep_rates) == 1:
                icpms_ws.cell(row=groupCounter, column=self.icpms_col_indices['laserMetadata_laserRepetitionRate'], value=rep_rates.pop())
            else:
                print(f'Not exporting rep rates for group {group.name} because there are multiple rep rates reported: {rep_rates}')
            if len(spot_heights) == 1 and len(spot_widths) == 1:
                spot_height = next(iter(spot_heights))
                spot_width = next(iter(spot_widths))

                if spot_height == spot_width:
                    icpms_ws.cell(
                        row=groupCounter,
                        column=self.icpms_col_indices['laserMetadata_laserSpotSize'],
                        value=spot_height
                    )
            else:
                print(f'Not exporting spot size for group {group.name} because there are multiple spot heights reported: {spot_heights} or multiple spot widths reported: {spot_widths}')

            # Export other laser metadata
            icpms_ws.cell(row=groupCounter, column=self.icpms_col_indices['laserMetadata_laserFluence'], value=self.laserFluenceDoubleSpinBox.value)
            # icpms_ws.cell(row=groupCounter, column=icpms_col_indices['laserMetadata_laserWavelength'], value=self.laserWavelengthComboBox.currentText)
            icpms_ws.cell(row=groupCounter, column=self.icpms_col_indices['laserMetadata_pulseWidthValue'], value=self.laserPulseWidthComboBox.currentText)
            icpms_ws.cell(row=groupCounter, column=self.icpms_col_indices['laserMetadata_cellModel'], value=self.cellModelComboBox.currentText)
            icpms_ws.cell(row=groupCounter, column=self.icpms_col_indices['laserMetadata_carrierGas'], value=self.carrierGasComboBox.currentText)
            icpms_ws.cell(row=groupCounter, column=self.icpms_col_indices['laserMetadata_mixingDevice'], value=self.mixingDeviceComboBox.currentText)
            icpms_ws.cell(row=groupCounter, column=self.icpms_col_indices['laserMetadata_carrierGasFlow'], value=self.carrierGasFlowRateSpinBox.value)

            groupCounter += 1

        wb.save(fp)

        # Now use the system open the file in Excel
        subprocess.run(['open', export_filepath])
        self.close()

configuration = CONFIG_MANAGER.choose()
if configuration is not None:
    config_path, config = configuration
    widget = SettingsWidget(config_path, config)
    widget.show()
    loop = QEventLoop()
    widget.destroyed.connect(lambda *args: loop.quit())
    loop.exec() # wait ...
    print('finished')
