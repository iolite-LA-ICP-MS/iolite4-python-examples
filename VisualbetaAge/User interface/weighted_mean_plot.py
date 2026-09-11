#/ Type: UI
#/ Name: Weighted Average Plot
#/ Authors: Pro Shitou Wu
#/ Description: Calculate weighted average and MSWD for selections
#/ Version: 0.2
#/ Contact: shitou.wu@mail.iggcas.ac.cn, info@chemlabcorp.com

from iolite import QtGui
from iolite.QtGui import (
    QAction, QMessageBox, QWidget, QLabel, QVBoxLayout, QHBoxLayout,
    QSizePolicy)
from iolite.QtCore import Qt, Signal

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

import numpy as np
from scipy import stats
import math
from functools import partial

def createUIElements():
    action = QAction('Weighted Mean/Plot', ui)
    action.triggered.connect(show_widget)
    ui.setMenuName(['VisualAge_β decay'])
    ui.setAction(action)

def show_widget():
    widget = WeightedAverageWidget()
    widget.show()
    widget.update_plot()

def format_uncertainty(value, uncertainty, unit=" Ma"):
    """
    格式化数值和误差，年龄用4位有效数字，误差最低位跟年龄最低位齐平
    
    示例:
        123.4 ± 5.6 Ma
        1234 ± 66 Ma
        0.1234 ± 0.0056 Ma
    """
    if np.isnan(value) or np.isnan(uncertainty):
        return f"{value:.6f} ± {uncertainty:.6f}{unit}"
    
    # 处理零值
    if abs(value) < 1e-12:
        return f"0.0 ± {uncertainty:.1f}{unit}"
    
    # 确定年龄的格式
    abs_value = abs(value)
    
    if abs_value >= 1000:
        # 四位数以上，保留整数
        value_fmt = f"{value:.0f}"
        uncertainty_fmt = f"{uncertainty:.0f}"
    elif abs_value >= 100:
        # 三位数，保留一位小数
        value_fmt = f"{value:.1f}"
        # 找到小数位数
        decimal_places = 1
        uncertainty_fmt = f"{uncertainty:.{decimal_places}f}"
    elif abs_value >= 10:
        # 两位数，保留两位小数
        value_fmt = f"{value:.2f}"
        decimal_places = 2
        uncertainty_fmt = f"{uncertainty:.{decimal_places}f}"
    elif abs_value >= 1:
        # 一位数，保留三位小数
        value_fmt = f"{value:.3f}"
        decimal_places = 3
        uncertainty_fmt = f"{uncertainty:.{decimal_places}f}"
    else:
        # 小于1的值，使用科学记数法或保留4位有效数字
        # 确定小数点后的位数
        if abs_value >= 0.1:
            decimal_places = 4
        elif abs_value >= 0.01:
            decimal_places = 5
        elif abs_value >= 0.001:
            decimal_places = 6
        else:
            decimal_places = 7
        
        value_fmt = f"{value:.{decimal_places}f}"
        
        # 对于非常小的值，可能需要调整不确定性的显示
        if uncertainty < 10 ** (-decimal_places):
            # 如果误差太小，增加显示位数
            uncertainty_fmt = f"{uncertainty:.{decimal_places+1}f}"
        else:
            uncertainty_fmt = f"{uncertainty:.{decimal_places}f}"
    
    return f"{value_fmt} ± {uncertainty_fmt}{unit}"


class WeightedAverageWidget(QWidget):
    Count = 0
    
    def __init__(self):
        QWidget.__init__(self)
        WeightedAverageWidget.Count += 1
        self.setWindowTitle(f'Weighted Mean[{WeightedAverageWidget.Count}]')
        
        # UI elements
        self.label = QLabel()
        self.canvas = WeightedAverageCanvas()
        self.canvas.got_argmin.connect(self.disp_sel_name)
        
        self.settings_widget = WeightedAverageSettingsWidget()
        self.settings_widget.channel_changed.connect(self.update_plot)

        self.color_setting_w = ColorChannelWidget()
        self.color_setting_w.color_channel_changed.connect(self.update_plot)
        self.color_setting_w.refresh_btn.clicked.connect(self.reset_color_range)

        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.settings_widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.color_setting_w.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        
        # Layout
        vbox = QVBoxLayout()
        vbox.addWidget(self.canvas)
        vbox.addWidget(self.label)

        # 设置区域布局
        hbox = QHBoxLayout()
        hbox.addWidget(self.settings_widget)
        hbox.addStretch()
        vbox.addLayout(hbox)
        
        # 颜色设置区域
        h2box = QHBoxLayout()
        h2box.addWidget(self.color_setting_w)
        h2box.addStretch()
        vbox.addLayout(h2box)

        self.setLayout(vbox)
        
        # Connect signals
        data.activeSelectionChanged.connect(self.update_plot)
        data.selectionChanged.connect(self.update_plot)
        data.activeSelectionGroupChanged.connect(self.update_plot)
        
        self.resize(800, 480)
        self.setWindowFlags(Qt.WindowStaysOnTopHint)

        # 颜色范围相关变量
        self.color_range_initialized = False
        self.color_min = 0
        self.color_max = 1
        self.cbar = None

    def update_plot(self):
        if not self.isVisible(): 
            return
            
        settings = self.settings_widget.getSettings()
        try:
            channel = data.timeSeries(settings['channel'])
        except Exception as e:
            return
            
        grp = data.activeSelectionGroup()
        if not grp:
            return
            
        selections = grp.selections()
        if not selections:
            return

        # 获取颜色通道设置
        color_settings = self.color_setting_w.getSettings()
        color_channel_name = color_settings.get('color_channel')
        color_channel = None
        if color_channel_name and color_channel_name != 'fixed_black':
            try:
                color_channel = data.timeSeries(color_channel_name)
            except:
                color_channel = None
        
        # 初始化颜色范围
        if not self.color_range_initialized and color_channel:
            self.initialize_color_range(selections, color_channel)

        # Extract values and uncertainties
        values = []
        uncertainties = []
        names = []
        colors = []
        
        for sel in selections:
            result = data.result(sel, channel)
            values.append(result.value())
            uncertainties.append(result.uncertaintyAs2SE())
            names.append(sel.name)

            # 获取颜色
            if color_channel:
                try:
                    ppm_value = data.result(sel, color_channel).value()
                    # 使用固定范围进行归一化
                    norm_value = np.clip(ppm_value, self.color_min, self.color_max)
                    normalized = (norm_value - self.color_min) / (self.color_max - self.color_min)
                    colors.append(plt.cm.viridis(normalized))
                except:
                    colors.append('black')  # 出错时使用黑色
            else:
                colors.append('black')  # 固定颜色

        values = np.array(values)
        uncertainties = np.array(uncertainties)
        
        # Calculate weighted average
        weights = 1.0 / (uncertainties ** 2)
        weighted_avg = np.sum(weights * values) / np.sum(weights)
        weighted_avg_err = 1.0 / np.sqrt(np.sum(weights))
        
        # Calculate MSWD
        if len(values) > 1:
            chi2 = np.sum(weights * (values - weighted_avg) ** 2)
            mswd = chi2 / (len(values) - 1)
        else:
            mswd = float('nan')

        unit = settings.get('unit', '')

        # 更新颜色条
        if color_channel and self.color_range_initialized:
            self.update_colorbar()

        # Update plot
        self.canvas.clear_figure()
        self.canvas.plot_data(values, uncertainties, names,
                              weighted_avg, weighted_avg_err,
                              colors=colors, unit=unit)
        self.canvas.draw()
        
        # Update label with active selection info
        sel = data.activeSelection()
        if sel is None:
            self.label.setText(f'Weighted Mean: {format_uncertainty(weighted_avg, weighted_avg_err, unit)}, MSWD: {mswd:.2f}')
        else:
            result = data.result(sel, channel)
            sel_value = result.value()
            sel_uncertainty = result.uncertaintyAs2SE()
            formatted = format_uncertainty(sel_value, sel_uncertainty)
            self.label.setText(
                f'Active: {sel.name} = {format_uncertainty(sel_value, sel_uncertainty, unit)}, '
                f'Weighted Mean: {format_uncertainty(weighted_avg, weighted_avg_err, unit)}, MSWD: {mswd:.2f}'
            )

# colors
    def initialize_color_range(self, selections, color_channel):
        """初始化颜色范围（固定范围）"""
        ppm_values = []
        for sel in selections:
            try:
                value = data.result(sel, color_channel).value()
                ppm_values.append(value)
            except:
                continue

        if len(ppm_values) > 0:
            self.color_min = min(ppm_values)
            self.color_max = max(ppm_values)
            # 避免除零错误
            if abs(self.color_max - self.color_min) < 1e-10:
                self.color_max = self.color_min + 1.0
            self.color_range_initialized = True
            
            # 更新颜色设置部件显示范围
            self.color_setting_w.set_range_label(
                f"{self.color_min:.2f} - {self.color_max:.2f}")
    
    def update_colorbar(self):
        """更新或创建colorbar"""
        if self.cbar is not None:
            try:
                norm = matplotlib.colors.Normalize(
                    vmin=self.color_min, 
                    vmax=self.color_max
                )
                
                if hasattr(self.cbar, 'mappable') and self.cbar.mappable:
                    self.cbar.mappable.set_norm(norm)
                    self.cbar.update_normal(self.cbar.mappable)
                else:
                    sm = plt.cm.ScalarMappable(
                        cmap=plt.cm.viridis, 
                        norm=norm
                    )
                    sm.set_array([])
                    self.cbar.mappable = sm
                    self.cbar.update_normal(sm)
                
                self.cbar.draw_all()
                
            except Exception as e:
                print(f"更新colorbar时出错: {e}")
                self.recreate_colorbar()
        else:
            self.recreate_colorbar()
        
        self.canvas.draw_idle()
    
    def recreate_colorbar(self):
        """重建colorbar"""
        if self.cbar is not None:
            try:
                self.cbar.remove()
            except:
                pass
        
        norm = matplotlib.colors.Normalize(
            vmin=self.color_min, 
            vmax=self.color_max
        )
        sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=norm)
        sm.set_array([])
        
        self.cbar = self.canvas.figure.colorbar(sm, ax=self.canvas.axes)
        self.cbar.set_label('Concentration (ppm)')
    
    def reset_color_range(self):
        """重置颜色范围"""
        self.color_range_initialized = False
        self.update_plot()

####

    def disp_sel_name(self, idx):
        grp = data.activeSelectionGroup()
        if not grp:
            return
            
        selections = grp.selections()
        if idx < 0 or idx >= len(selections):
            return
            
        selection = selections[idx]
        settings = self.settings_widget.getSettings()
        
        try:
            channel = data.timeSeries(settings['channel'])
            result = data.result(selection, channel)
            value = result.value()
            uncertainty = result.uncertaintyAs2SE()
            unit = settings.get('unit', '')
            self.label.setText(f'{selection.name}: {format_uncertainty(value, uncertainty, unit)}')
        except Exception:
            self.label.setText(f'{selection.name}')

class WeightedAverageCanvas(FigureCanvasQTAgg):
    got_argmin = Signal(int)
    
    def __init__(self, parent=None, width=4, height=1.6, dpi=80):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        FigureCanvasQTAgg.__init__(self, fig)
        self.setParent(parent)
        FigureCanvasQTAgg.updateGeometry(self)
        fig.canvas.mpl_connect('motion_notify_event', self.on_hover)
        
        self.x_positions = None
        self.active_idx = -1
    
    def clear_figure(self):
        self.axes.cla()
        self.x_positions = None
        self.active_idx = -1
    
    def plot_data(self, values, uncertainties, names,
                  weighted_avg, weighted_avg_err, colors=None, unit=''):
        n = len(values)
        x_pos = np.arange(n)  # Use index as x position
        
        # Store for hover detection
        self.x_positions = x_pos

        # 如果没有提供颜色，使用默认黑色
        if colors is None:
            colors = ['black'] * n

        # Plot individual points with error bars with given colors
##        self.axes.errorbar(
##            x_pos, values, yerr=uncertainties, 
##            fmt='o', capsize=3, label='Individual measurements'
##        )
        for i in range(n):
            self.axes.errorbar(
                x_pos[i], values[i], yerr=uncertainties[i],
                fmt='o', capsize=3, color=colors[i],
                ecolor=colors[i],  # 误差条颜色
                elinewidth=2  # 误差条线宽
            )

        # Plot weighted average line
        weighted_avg_formatted = format_uncertainty(weighted_avg, weighted_avg_err, unit)
        self.axes.axhline(
            y=weighted_avg, color='red', linestyle='-', 
            label=f'Weighted avg: {weighted_avg_formatted}'
        )
        
        # Plot weighted average uncertainty band
        self.axes.axhspan(
            weighted_avg - weighted_avg_err, 
            weighted_avg + weighted_avg_err, 
            color='red', alpha=0.2
        )
        
        # Highlight active selection if any
        sel = data.activeSelection()
        if sel:
            try:
                idx = names.index(sel.name)
##                self.axes.plot(x_pos[idx], values[idx], 'o', markersize=8, 
##                              markerfacecolor='none', markeredgecolor='orange', 
##                              markeredgewidth=2)
                self.axes.errorbar(
                    x_pos[idx], values[idx], yerr=uncertainties[idx],
                    fmt='o', capsize=3, color='orange', 
                    markerfacecolor='none', markersize=8,
                    markeredgewidth=2, label='Active selection')
                self.active_idx = idx
            except ValueError:
                self.active_idx = -1
        
        # Calculate MSWD for title
        if len(values) > 1:
            weights = 1.0 / (uncertainties ** 2)
            chi2 = np.sum(weights * (values - weighted_avg) ** 2)
            mswd = chi2 / (len(values) - 1)
        else:
            mswd = float('nan')
    
        # Set labels and title
        settings_widget = self.parent().settings_widget
        settings = settings_widget.getSettings()
        channel_name = settings['channel']
        
        self.axes.set_xlabel('Measurement Index')
        self.axes.set_ylabel(channel_name)
##        self.axes.set_title(f'Weighted Mean of {channel_name}')

        # Create title with weighted mean and MSWD information
        title_text = f'Weighted Mean of {channel_name}\nWeighted Mean: {weighted_avg_formatted}, MSWD: {mswd:.2f}'
        self.axes.set_title(title_text)
        
        # Set x ticks to show selection names
        self.axes.set_xticks(x_pos)
        self.axes.set_xticklabels(names, rotation=45, ha='right')
        
        # Add legend
        self.axes.legend()
        
        # Adjust layout to prevent label cutoff
        plt.tight_layout()
    
    def on_hover(self, event):
        if self.x_positions is None:
            return
            
        x = event.xdata
        if x is None:
            return
            
        # Find closest x position
        idx = np.argmin(np.abs(self.x_positions - x))
        self.got_argmin.emit(idx)

class WeightedAverageSettingsWidget(QWidget):
    channel_changed = Signal()
    
    def __init__(self):
        QWidget.__init__(self)
        self.settings = {}
        
        formLayout = QtGui.QFormLayout()
        formLayout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(formLayout)
        
        # Channel selection
        channels = data.timeSeriesNames()
        self.channel_cb = QtGui.QComboBox(self)
        self.channel_cb.addItems(channels)
        
        if channels:
            self.channel_cb.setCurrentText(channels[0])
            self.setSetting('channel', channels[0])
        
        self.channel_cb.currentTextChanged.connect(partial(self.setSetting, 'channel'))
        self.channel_cb.currentTextChanged.connect(self.emit_channel_changed)
        
        formLayout.addRow('Channel', self.channel_cb)

        # 单位输入框
        self.unit_edit = QtGui.QLineEdit(self)
        self.unit_edit.setPlaceholderText("e.g., Ma, yr, cps")
        self.unit_edit.setText(" Ma")  # 默认值，注意前面的空格
        self.setSetting('unit', " Ma")
        self.unit_edit.setMaximumWidth(100)
        self.unit_edit.returnPressed.connect(self.on_unit_changed)
        self.unit_edit.editingFinished.connect(self.on_unit_changed)

        formLayout.addRow('Unit', self.unit_edit)


    def getSettings(self):
        return self.settings
    
    def setSetting(self, key, value):
        self.settings[key] = value
    
    def emit_channel_changed(self, text):
        self.channel_changed.emit()
    
    def update_cb(self):
        channels = data.timeSeriesNames()
        self.channel_cb.blockSignals(True)
        self.channel_cb.clear()
        self.channel_cb.addItems(channels)
        self.channel_cb.blockSignals(False)

    def on_unit_changed(self, *args):
        """单位输入框内容改变时的处理"""
        unit_text = self.unit_edit.text.strip()
        # 确保单位前有一个空格（除非单位为空）
        if unit_text and not unit_text.startswith(" "):
            unit_text = " " + unit_text
        elif not unit_text:
            unit_text = ""
            
        self.unit_edit.setText(unit_text)
        self.setSetting('unit', unit_text)
        self.emit_channel_changed(None)

class ColorChannelWidget(QWidget):
    color_channel_changed = Signal()
    
    def __init__(self):
        QWidget.__init__(self)
        self.settings = {'color_channel': 'fixed_black'}
        
        thebox = QtGui.QHBoxLayout()
        thebox.setContentsMargins(0, 0, 0, 0)
        self.setLayout(thebox)
        
        # 颜色通道选择下拉框
        self.color_cb = QtGui.QComboBox(self)
        self.color_cb.currentTextChanged.connect(self.on_color_channel_changed)
        thebox.addWidget(QLabel("Color Source:"))
        thebox.addWidget(self.color_cb)
        
        # 范围显示标签
        self.range_label = QLabel("Color Range: -")
        thebox.addWidget(self.range_label)
        
        # 刷新按钮
        self.refresh_btn = QtGui.QPushButton("Refresh Color Range")
        self.refresh_btn.clicked.connect(self.emit_refresh)
        thebox.addWidget(self.refresh_btn)
        
        # 初始更新
        self.update_channels()

    def update_channels(self):
        """更新颜色通道下拉列表"""
        self.color_cb.blockSignals(True)
        current = self.color_cb.currentText
        self.color_cb.clear()
        
        channels = data.timeSeriesNames()
        
        # 添加固定颜色选项
        self.color_cb.addItem("fixed_black", "固定黑色")
        
        # 添加分隔线
        self.color_cb.insertSeparator(1)
        
        # 添加ppm相关通道（优先显示）
        ppm_n_channels = []
        ppm_channels = []
        other_channels = []
        
        for ch in channels:
            if '_ppm_n' in ch.lower():
                ppm_n_channels.append(ch)
            elif '_ppm' in ch.lower():
                ppm_channels.append(ch)
            else:
                other_channels.append(ch)

        sorted_ppm_channels = sorted(ppm_n_channels) + sorted(ppm_channels)

        # 添加ppm通道
        for ch in sorted_ppm_channels:
            self.color_cb.addItem(ch)
        
        if sorted_ppm_channels:
            self.color_cb.insertSeparator(len(sorted_ppm_channels) + 2)
        
        # 添加其他通道
        for ch in sorted(other_channels):
            self.color_cb.addItem(ch)
        
        # 恢复之前的选择
        if current and self.color_cb.findText(current) >= 0:
            self.color_cb.setCurrentText(current)
        elif sorted_ppm_channels:
            self.color_cb.setCurrentText(sorted_ppm_channels[0])
            self.settings['color_channel'] = sorted_ppm_channels[0]
        else:
            self.color_cb.setCurrentIndex(0)  # 默认固定颜色
            self.settings['color_channel'] = None
        
        self.color_cb.blockSignals(False)
        self.on_color_channel_changed(self.settings['color_channel'])

    def on_color_channel_changed(self, text):
        """颜色通道变化时的处理"""
        if text == 'fixed_black':
            self.settings['color_channel'] = None
            self.refresh_btn.setVisible(False)
            self.range_label.setText("Color Range: -")
        else:
            self.settings['color_channel'] = text
            self.refresh_btn.setVisible(True)
        
        self.color_channel_changed.emit()

    def getSettings(self):
        return self.settings.copy()

    def set_range_label(self, text):
        """设置范围显示标签"""
        self.range_label.setText(f"Color Range: {text}")

    def emit_refresh(self):
        """发射刷新信号"""
        self.color_channel_changed.emit()
