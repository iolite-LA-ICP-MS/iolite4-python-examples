#/ Type: DRS
#/ Name: in situ Beta system dating v0.6
#/ Authors: Shitou Wu, Angning Yu
#/ Description: A basic Lu-Hf isotope DRS, based on Wu et al. JAAS. 2023 (Pb added)
#/ References: Wu et al. JAAS. 2023
#/ Version: 0.6
#/ Contact: shitou.wu@mail.iggcas.ac.cn

email = "shitou.wu@mail.iggcas.ac.cn"
copyright_message = f'''
Descriptions:
written by Shitou Wu, Angning Yu
contact Email: {email}
'''

SYSTEMS = [
    'LuHf',
    'RbSr',
    'KCa',
    'PbPb',
    'ReOs',
    'MoDau',
    'LaBa',
    ]
SYSTEM_DRS_DICT = {}
SYSTEM_WIDGET_DICT = {}

from iolite import QtGui
from iolite.types import Result
from iolite.QtGui import QSizePolicy, QWidget
from iolite.QtCore import QSettings
import numpy as np
import pprint
from functools import partial
import itertools
import os, sys, time

########
from iolite.QtGui import QHBoxLayout, QVBoxLayout, QWidget, QPushButton, QLabel
from iolite.QtGui import QLineEdit, QSizePolicy

def recursiveInitUI(widgets, ori='V', add_stretch=False):
    if ori=='V':
        outerbox, parent, child_ori = QVBoxLayout(), QVBoxLayout(), 'H'
    else:
        outerbox, parent, child_ori = QHBoxLayout(), QHBoxLayout(), 'V'
    parent.setSpacing(0)
    parent.setContentsMargins(12,2,12,2)
    for w in widgets:
        if isinstance(w, list):
            child_layout = recursiveInitUI(w, child_ori)
            parent.addLayout(child_layout)
        elif isinstance(w, str):
            parent.addWidget(QLabel(w))
        elif hasattr(w, 'setVisible'):
            parent.addWidget(w)
        else:
            parent.addLayout(w)
    if add_stretch:
        parent.addStretch()
    return parent



def expand_agg(data, channel_name, agg_type='median', groups=[]):
    chn = data.timeSeries(channel_name)
    signal = chn.data()
    all_sel = itertools.chain(
        *[data.selectionGroup(grp).selections() for grp in groups])
    X = signal.copy()
    basegrp = data.selectionGroupList(data.Baseline)[0]
    base = data.spline(basegrp.name, channel_name).data()
    for sel in all_sel:
        ids = chn.selectionIndices(sel)
        sel_value = signal[ids]
        sel_base = base[ids]
        if agg_type == 'mean':
            X[ids] = np.nanmean(sel_value)
        elif agg_type == 'median':
            X[ids] = np.nanmedian(sel_value)
        elif agg_type == 'max':
            X[ids] = np.nanmax(sel_value)
        elif agg_type == 'min':
            X[ids] = np.nanmin(sel_value)
        elif agg_type == 'positive':
            neg_idx = np.where(sel_value - sel_base <=0)[0]
            median = np.nanmedian(sel_value)
            sel_value[neg_idx] = median
            X[ids] = sel_value
        else: agg_value = 0
    return X
from iolite.QtCore import Signal

##ALL_AGG_TYPES = ['median', 'mean', 'max', 'min', 'positive']
##class ExpandAggWidget(QWidget):
##    channel_created = Signal(str)
##    def __init__(self, data, drs):
##        super().__init__(self)
##        self.data = data
##        self.drs = drs
##        self.drs_settings = {}
##        formLayout = QtGui.QFormLayout()
##
##
##        timeSeriesNames = self.data.timeSeriesNames(self.data.Input)
##        defaultChannelName = ""
##        if timeSeriesNames:
##            defaultChannelName = timeSeriesNames[0]
##        self.set_setting('channelName', defaultChannelName)
##        self.set_setting('aggType', ALL_AGG_TYPES[0])
##        self.set_setting('groups', [])
##
##        channelComboBox = QtGui.QComboBox()
##        timeSeriesNames = self.data.timeSeriesNames()
##        channelComboBox.addItems(timeSeriesNames)
##        channelComboBox.setCurrentText(self.drs_settings['channelName'])
##        channelComboBox.currentTextChanged.connect(
##            partial(self.set_setting, 'channelName'))
##        formLayout.addRow('channelName', channelComboBox)
##        self.channelComboBox = channelComboBox
##
##        aggTypeComboBox = QtGui.QComboBox()
##        aggTypeComboBox.addItems(ALL_AGG_TYPES)
##        aggTypeComboBox.setCurrentText(self.drs_settings['aggType'])
##        aggTypeComboBox.currentTextChanged.connect(
##            partial(self.set_setting, 'aggType'))
##        formLayout.addRow('aggType', aggTypeComboBox)
##
##        groups = self.data.selectionGroupNames(self.data.ReferenceMaterial | self.data.Sample)
##        for grp in groups:
##            checkbox = QtGui.QCheckBox(grp)
##            checkbox.toggled.connect(partial(self.toggle_group, grp))
##            formLayout.addRow(grp, checkbox)
##        createChannelBtn = QPushButton('createChannel')
##        createChannelBtn.clicked.connect(self.create_channel)
##        self.setLayout(recursiveInitUI([
##            formLayout,
##            createChannelBtn,
##            ], add_stretch=1))
##        self.channel_created.connect(self.update_channels)
##
##    def update_channels(self, new_channel):
##        old_text = self.channelComboBox.currentText
##        self.channelComboBox.clear()
##        timeSeriesNames = self.data.timeSeriesNames()
##        self.channelComboBox.addItems(timeSeriesNames)
##        self.channelComboBox.setCurrentText(old_text)
##
##    def set_setting(self, key, value):
##        self.drs_settings[key] = value
##        self.drs.setSetting(key, value)
##    def toggle_group(self, name, flag):
##        groups = self.drs_settings['groups']
##        if flag:
##            if not name in groups: groups.append(name)
##        else:
##            if name in groups: groups.remove(name)
##        self.drs.setSetting('groups', groups)
##
##    def create_channel(self):
##        name = self.drs_settings['channelName']
##        agg_type = self.drs_settings['aggType']
##        print('create_channel', 'name', name, 'agg_type', agg_type)
##        groups = self.drs_settings['groups']
##        chn = self.data.timeSeries(name)
##        X = expand_agg(self.data, name, agg_type, groups)
##        new_name = f'{name}_{agg_type}'
##        self.data.createTimeSeries(new_name, self.data.Input, chn.time(), X)
##        self.channel_created.emit(new_name)

class IsotopeWidget(QWidget):
    channel_created = Signal(str)
    def __init__(self, data, chnA, chnB, factor):
        QWidget.__init__(self)
        self.data = data
        self.chnA = chnA
        self.chnB = chnB
        self.BAfactor = factor
        Btn = QPushButton(f'create {chnA} by channel/factor')
        self.b_input = QtGui.QComboBox(self)
##        timeSeriesNames = self.data.timeSeriesNames()
        self.b_input.addItems(self.data.timeSeriesNames())
        self.factor_edit = QLineEdit()
        self.factor_edit.setText('{0:.4f}'.format(self.BAfactor))
        Btn.clicked.connect(self.doAB)
        self.setLayout(recursiveInitUI([
            Btn, self.b_input, self.factor_edit, 
            ]))
    def doAB(self):
        try:
            factor = float(self.factor_edit.text)
            chnB = str(self.b_input.currentText)
            chn = self.data.timeSeries(chnB)
            X = chn.data() / factor
            new_name = f'{self.chnA}from{chnB}'
            self.data.createTimeSeries(
                new_name, self.data.Input, chn.time(), X)
            self.channel_created.emit(new_name)
        except Exception as e:
            print('IsotopeWidget.doAB', type(e), e)


########

here_path = QSettings().value('paths/DataReductionSchemesPath')
if not here_path:
    raise ValueError("DataReductionSchemesPath not found in QSettings")
here_path = os.path.abspath(os.path.normpath(here_path))

# 验证路径是否存在
if not os.path.exists(here_path):
    raise FileNotFoundError(f"DataReductionSchemesPath does not exist: {here_path}")

# 验证必要的文件和目录
##expand_agg_path = os.path.join(here_path, 'expand_agg.py')
beta_system_dir = os.path.join(here_path, 'beta_system')

##if not os.path.exists(expand_agg_path):
##    raise FileNotFoundError(f"expand_agg.py not found at: {expand_agg_path}")
if not os.path.exists(beta_system_dir):
    raise FileNotFoundError(f"beta_system directory not found at: {beta_system_dir}")

# 临时添加路径到 sys.path 以处理相对导入
original_sys_path = sys.path.copy()
try:
    # 添加必要的路径
    sys.path.insert(0, here_path)
    sys.path.insert(0, os.path.join(here_path, 'beta_system'))
    
    # 现在可以正常导入了
##    from expand_agg import ExpandAggWidget, IsotopeWidget, recursiveInitUI
    if 'LuHf' in SYSTEMS:
        try:
            from beta_system.luhf import LuHfDRS as beta_drs_class
            beta_drs = beta_drs_class(data, drs)
            SYSTEM_DRS_DICT['LuHf'] = beta_drs
            beta_drs.set_logger(IoLog)
            SYSTEM_WIDGET_DICT['LuHf'] = beta_drs.settingsWidget()
        except Exception as e:
            print(f'importimg Luhf failed with {e}')

    if 'RbSr' in SYSTEMS:
        try:
            from beta_system.rbsr import RbSrDRS as beta_drs_class
            beta_drs = beta_drs_class(data, drs)
            SYSTEM_DRS_DICT['RbSr'] = beta_drs
            beta_drs.set_logger(IoLog)
            SYSTEM_WIDGET_DICT['RbSr'] = beta_drs.settingsWidget()
        except Exception as e:
            print(f'importimg RbSr failed with {e}')

    if 'KCa' in SYSTEMS:
        try:
            from beta_system.kca import KCaDRS as beta_drs_class
            beta_drs = beta_drs_class(data, drs)
            SYSTEM_DRS_DICT['KCa'] = beta_drs
            beta_drs.set_logger(IoLog)
            SYSTEM_WIDGET_DICT['KCa'] = beta_drs.settingsWidget()
        except Exception as e:
            print(f'importimg KCa failed with {e}')

    if 'ReOs' in SYSTEMS:
        try:
            from beta_system.reos import ReOsDRS as beta_drs_class
            beta_drs = beta_drs_class(data, drs)
            SYSTEM_DRS_DICT['ReOs'] = beta_drs
            beta_drs.set_logger(IoLog)
            SYSTEM_WIDGET_DICT['ReOs'] = beta_drs.settingsWidget()
        except Exception as e:
            print(f'importimg ReOs failed with {e}')

    if 'PbPb' in SYSTEMS:
        try:
            from beta_system.pbpb import PbPbDRS as beta_drs_class
            beta_drs = beta_drs_class(data, drs)
            SYSTEM_DRS_DICT['PbPb'] = beta_drs
            beta_drs.set_logger(IoLog)
            SYSTEM_WIDGET_DICT['PbPb'] = beta_drs.settingsWidget()
        except Exception as e:
            print(f'importimg PbPb failed with {e}')

    if 'MoDau' in SYSTEMS:
        try:
            from beta_system.modau import MoDauDRS as beta_drs_class
            beta_drs = beta_drs_class(data, drs)
            SYSTEM_DRS_DICT['MoDau'] = beta_drs
            beta_drs.set_logger(IoLog)
            SYSTEM_WIDGET_DICT['MoDau'] = beta_drs.settingsWidget()
        except Exception as e:
            print(f'importimg MoDau failed with {e}')

    if 'LaBa' in SYSTEMS:
        try:
            from beta_system.laba import LaBaDRS as beta_drs_class
            beta_drs = beta_drs_class(data, drs)
            SYSTEM_DRS_DICT['LaBa'] = beta_drs
            beta_drs.set_logger(IoLog)
            SYSTEM_WIDGET_DICT['LaBa'] = beta_drs.settingsWidget()
        except Exception as e:
            print(f'importimg LaBa failed with {e}')

finally:
    # 恢复原来的 sys.path
    sys.path = original_sys_path

def runDRS():
    drs.message("Starting beta system isotopes DRS...")
    drs.progress(0)
    drs.progress(1)

    if SYSTEM_DRS_DICT:
        beta_drs = list(SYSTEM_DRS_DICT.values())[0]
        beta_drs.baseline_sub()
    for beta_system in SYSTEMS:
        if beta_system not in SYSTEM_DRS_DICT:
            continue
        beta_drs = SYSTEM_DRS_DICT[beta_system]
        if beta_drs.should_run():
            beta_drs.run()

    drs.message("Finished!")
    drs.progress(80)
    drs.progress(100)
    drs.finished()

from iolite.QtGui import (QWidget, QHBoxLayout, QVBoxLayout, QScrollArea, 
                             QFrame, QLabel, QPushButton, QGroupBox)
from iolite.QtCore import Qt

class SettingsWidget(QWidget):
    def __init__(self, pre_w, system_widgets, parent=None):
        super().__init__(parent)
        self.pre_w = pre_w
        self.system_widgets = list(system_widgets)
        self.init_ui()
    
    def init_ui(self):
        # 创建主布局（水平布局）
        main_layout = QHBoxLayout(self)
        main_layout.setContentsMargins(10, 10, 10, 10)
        main_layout.setSpacing(15)
        
        # ========== 左侧：pre_w 区域 ==========
        vbox = QVBoxLayout()
        vbox.addWidget(self.pre_w)
        vbox.addStretch()
##        main_layout.addWidget(self.pre_w, 1)  # 设置拉伸因子为1
        main_layout.addLayout(vbox, 1)
        
        # ========== 右侧：滚动区域 ==========
        # 创建滚动区域
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)  # 允许内容自适应大小
        scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)  # 禁用水平滚动条
        scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)     # 根据需要显示垂直滚动条
        
        # 创建滚动区域的内容容器
        scroll_content = QWidget()
        scroll_content.setObjectName("scroll_content")
        
        # 创建内容布局（垂直布局）
        content_layout = QVBoxLayout(scroll_content)
        content_layout.setContentsMargins(10, 10, 10, 10)
        content_layout.setSpacing(20)
        content_layout.setAlignment(Qt.AlignTop)  # 顶部对齐

        for w in self.system_widgets:
            content_layout.addWidget(w)
        
        # 添加弹簧，将内容推到顶部（可选）
        content_layout.addStretch()

        copyright_label = QLabel()
##        copyright_label.setTextFormat(Qt.RichText)
        copyright_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        copyright_label.setOpenExternalLinks(True)
        copyright_label.setText(copyright_message.strip())
        copyright_label.setStyleSheet("""
            QLabel {
                color: #888;
                font-family: 'Courier New', monospace;
                font-size: 11px;
                background-color: #1e1e1e;
                padding: 12px;
                border: 1px solid #333;
                border-radius: 6px;
                margin-top: 20px;
            }
        """)
        copyright_label.setToolTip("Select and copy (Ctrl+C)")
        content_layout.addWidget(copyright_label)

        # 将内容容器设置到滚动区域
        scroll_area.setWidget(scroll_content)
        
        # 将滚动区域添加到主布局，设置拉伸因子为3（右侧占更多空间）
        main_layout.addWidget(scroll_area, 4)
        
        # 设置窗口样式
        self.setStyleSheet("""
            QWidget#scroll_content {
                background-color: #000;
                border-radius: 5px;
            }
            QScrollArea {
                border: 1px solid #ddd;
                border-radius: 5px;
                background-color: white;
            }
            QGroupBox {
                font-weight: bold;
                border: 1px solid #ccc;
                border-radius: 5px;
                margin-top: 10px;
                padding-top: 10px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px 0 5px;
            }
        """)
    
    
    def add_system_widget(self, widget):
        """动态添加体系设置控件到滚动区域"""
        scroll_content = self.findChild(QWidget, "scroll_content")
        if scroll_content:
            layout = scroll_content.layout()
            # 在弹簧之前插入
            layout.insertWidget(layout.count() - 1, widget)
            self.system_widgets.append(widget)
    
    def remove_system_widget(self, widget):
        """动态移除体系设置控件"""
        scroll_content = self.findChild(QWidget, "scroll_content")
        if scroll_content and widget in self.system_widgets:
            layout = scroll_content.layout()
            layout.removeWidget(widget)
            widget.deleteLater()
            self.system_widgets.remove(widget)

def settingsWidget():
    sr86_107_w = IsotopeWidget(data, 'Sr86', 'Sr107', 8.374)
    hf177_w = IsotopeWidget(data, 'Hf177', 'Hf177', 1.8865)

    for w in [sr86_107_w, hf177_w,]:
        for beta_drs in SYSTEM_DRS_DICT.values():
            w.channel_created.connect(beta_drs.update_channels)
    beta_widgets = [SYSTEM_WIDGET_DICT[beta_system] for beta_system in SYSTEMS
                    if beta_system in SYSTEM_WIDGET_DICT]
    pre_w = QWidget()
    calc_gbox = QtGui.QGroupBox('Calculated Intensity')
    calc_gbox.setLayout(recursiveInitUI([sr86_107_w, hf177_w,]))
    pre_w.setLayout(recursiveInitUI([calc_gbox]))
    pre_w.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
    widget = SettingsWidget(pre_w, beta_widgets)

    drs.setSettingsWidget(widget)
