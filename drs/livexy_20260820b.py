#/ Type: UI
#/ Name: in situ Lu-Hf dating 0.51
#/ Authors: Pro Shitou Wu
#/ Description: A basic Lu-Hf isotope DRS, based on Wu et al. JAAS. 2023
#/ References: Wu et al. JAAS. 2023
#/ Version: 0.51
#/ Contact: shitou.wu@mail.iggcas.ac.cn


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
from scipy import stats, optimize
from functools import partial
import math, re
from iolite_helpers import weightedMean2D


METHODS = ['inversHf', 'isochHf', 'Pb76', 'isochSr', 'inversSr',
           'isochKCa', 'inversKCa',
           'isochBa', 'inversBa',
           'isochMoDau', 'inversMoDau',
           ]

PPM_FOR_METHODS = {
    'inversHf': ['Lu', 'Hf'],
    'isochHf': ['Lu', 'Hf'],
    'Pb76': [],
    'isochSr': ['Rb', 'Sr'],
    'inversSr': ['Rb', 'Sr'],
    'isochKCa': ['K', 'Ca'],
    'inversKCa': ['K', 'Ca'],
    'isochMoDau': [],
    'inversMoDau': [],
    'isochBa': ['La', 'Ba'],
    'inversBa': ['La', 'Ba'],
    }

SYSTEM_LAMBDAS = {
    'Lu-Hf': 1.865e-11,
    'Rb-Sr': 1.3972e-11,
    'K-Ca': 4.962e-10,
    'U-Pb': None,  # Pb-Pb体系不使用单一lambda
    'Mo-Dau': 1.865e-11,
    'La-Ba': 0.4362e-11,
}

METHOD_TO_SYSTEM = {
    'inversHf': 'Lu-Hf',
    'isochHf': 'Lu-Hf',
    'Pb76': 'U-Pb',
    'isochSr': 'Rb-Sr',
    'inversSr': 'Rb-Sr',
    'isochKCa': 'K-Ca',
    'inversKCa': 'K-Ca',
    'isochMoDau': 'Mo-Dau',
    'inversMoDau': 'Mo-Dau',
    'isochBa': 'La-Ba',
    'inversBa': 'La-Ba',
}

COLORMAP = plt.cm.viridis

def createUIElements():
    action = QAction('Live isochron v051', ui)
    action.triggered.connect(show_widget)
    ui.setMenuName(['VisualAge_β decay'])
    ui.setAction(action)

def show_widget():
    my_widget = ActiveSelWidget()
    my_widget.xyplot_setting_w.update_cb()
    my_widget.xyplot_setting_w.age_cb.setCurrentText('inversHf')
    my_widget.color_setting_w.update_for_method('inversHf')
    my_widget.show()
    my_widget.show_selection()


from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms

def confidence_ellipse(selection, chnx, chny, n_std, ax,
                       facecolor='none', **kwargs):
    x = chnx.dataForSelection(selection)
    y = chny.dataForSelection(selection)
    cov = np.cov(x, y)
##    cov_xy = cov[0, 1]
    pearson = cov[0, 1]/np.sqrt(cov[0, 0]*cov[1, 1])
    ell_radius_x = np.sqrt(1+pearson)
    ell_radius_y = np.sqrt(1-pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x*2, height=ell_radius_y*2, facecolor=facecolor, **kwargs)
    scale_x = data.result(selection, chnx).uncertaintyAs2SE()
    mean_x = data.result(selection, chnx).value()
    scale_y = data.result(selection, chny).uncertaintyAs2SE()
    mean_y = data.result(selection, chny).value()
    trans = transforms.Affine2D().rotate_deg(45).scale(scale_x, scale_y).translate(mean_x, mean_y)
    ellipse.set_transform(trans+ax.transData)
    return ax.add_patch(ellipse) #, cov_xy


def ellipse_from_errorbar(x, y, xerr, yerr, ax, facecolor='none', **kwargs):
    mean_x = x
    mean_y = y
    scale_x = xerr
    scale_y = yerr
    ellipse = Ellipse((0, 0), width=2, height=2, facecolor=facecolor, **kwargs)
    trans = transforms.Affine2D().rotate_deg(45).scale(scale_x, scale_y).translate(mean_x, mean_y)
    ellipse.set_transform(trans+ax.transData)
    return ax.add_patch(ellipse)


def format_ref_age(t_ma):
    """参考线年龄标注：<1 Ga 用 Ma 为单位，>=1 Ga（含）用 Ga 为单位。"""
    if t_ma >= 1000:
        ga = t_ma / 1000.0
        if abs(ga - round(ga)) < 1e-9:
            return f'{round(ga)} Ga'
        return f'{ga:g} Ga'
    return f'{t_ma:g} Ma'


class ActiveSelWidget(QWidget):
    Count = 0
    def __init__(self):
        QWidget.__init__(self)
        ActiveSelWidget.Count += 1
        self.setWindowTitle(f'Live XY-plot[{ActiveSelWidget.Count}]')
        self.label = QLabel()
        self.canvas = PlotCanvas()
        self.canvas.got_argmin.connect(self.disp_sel_name)
        self.xyplot_setting_w = XYPlotSettingWidget()
        self.xyplot_setting_w.channel_changed.connect(self.show_selection)
        self.fit_w = FittingOptionWidget()
        self.fit_w.setting_changed.connect(self.show_selection)

        self.color_setting_w = ColorChannelWidget()
        self.color_setting_w.color_channel_changed.connect(self.show_selection)

        self.canvas.setSizePolicy(
            QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.xyplot_setting_w.setSizePolicy(
            QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.color_setting_w.setSizePolicy(
            QSizePolicy.Expanding, QSizePolicy.Fixed)

        vbox = QVBoxLayout()
        hbox = QHBoxLayout()
        vbox.addWidget(self.canvas)
        vbox.addWidget(self.label)
        vbox.addLayout(hbox)
        hbox.addWidget(self.xyplot_setting_w)
        hbox.addWidget(self.fit_w)

        h2box = QHBoxLayout()
        vbox.addLayout(h2box)
        h2box.addWidget(self.color_setting_w)
        h2box.addStretch()
        self.rho_label = QLabel()
        vbox.addWidget(self.rho_label)
        self.setLayout(vbox)

        data.activeSelectionChanged.connect(self.show_selection)
        data.selectionChanged.connect(self.show_selection)
        data.activeSelectionGroupChanged.connect(self.show_selection)
        self.resize(800, 480)
        self.setWindowFlags(Qt.WindowStaysOnTopHint)
        self.xyplot_setting_w.got_y_interp.connect(self.fit_w.set_y_interp)
        self.xyplot_setting_w.age_cb.currentTextChanged.connect(
            self.color_setting_w.update_for_method)
        self.color_setting_w.refresh_btn.clicked.connect(self.reset_color_range)

        self.color_range_initialized = False
        self.color_min = 0
        self.color_max = 1
        self.cbar = None

##    def print_figure_info(self, step_name):
##        print(f"\n=== {step_name} ===")
##        print(f"Figure size: {self.canvas.figure.get_size_inches()}")
##        print(f"Figure DPI: {self.canvas.figure.dpi}")
##        print(f"Axes position: {self.canvas.axes.get_position()}")
##        
##        subplots_adjust_params = self.canvas.figure.subplotpars
##        print(f"Subplot params: left={subplots_adjust_params.left:.3f}, "
##              f"right={subplots_adjust_params.right:.3f}, "
##              f"bottom={subplots_adjust_params.bottom:.3f}, "
##              f"top={subplots_adjust_params.top:.3f}, "
##              f"wspace={subplots_adjust_params.wspace:.3f}, "
##              f"hspace={subplots_adjust_params.hspace:.3f}")
##        
##        print(f"Colorbar exists: {self.cbar is not None}")
##        if self.cbar:
##            print(f"Colorbar axes: {self.cbar.ax}")
##            print(f"Colorbar position: {self.cbar.ax.get_position()}")
##        
##        print(f"Number of axes in figure: {len(self.canvas.figure.axes)}")
##        for i, ax in enumerate(self.canvas.figure.axes):
##            print(f"  Axes {i}: {ax}, position: {ax.get_position()}")


    def show_selection(self):
        if not self.isVisible(): return

##        print("\n=== show_selection开始 ===")
##        self.print_figure_info("show_selection开始前")

        xy_settings = self.xyplot_setting_w.getSettings()
        try:
            chnx, chny = [data.timeSeries(xy_settings.get(f'{key}channel'))
                          for key in 'XY']
            print(f"X通道: {chnx.name}, Y通道: {chnx.name}")
        except Exception as e:
            print(f"获取通道失败: {e}")
            return

        fit_settings = self.fit_w.getSettings()
        grp = data.activeSelectionGroup()
        if not grp:
            print("没有活动的selection group")
            return

        color_settings = self.color_setting_w.getSettings()
        color_channel_name = color_settings.get('color_channel')
        color_channel = None
        if color_channel_name and color_channel_name != 'fixed_orange':
            try:
                color_channel = data.timeSeries(color_channel_name)
                print(f"颜色通道: {color_channel_name}")
            except:
                print(f"颜色通道 {color_channel_name} 不存在")
                color_channel = None
        
        print(f"颜色范围初始化状态: {self.color_range_initialized}")
        if not self.color_range_initialized and color_channel:
            self.initialize_color_range(grp.selections(), color_channel)

        self.canvas.clear_figure()
##        print("清除figure后...")
##        self.print_figure_info("清除figure后")

        self.canvas.set_ax_label(chnx.name, chny.name)
        values = self.extract_values(grp.selections(), chnx, chny)
        k, b = self.draw_linear_fit(*values[:2])

        if fit_settings['errorbar']:
            self.draw_errorbar(*values)

        if fit_settings['ellipse']:
            self.draw_colored_ellipses(
                grp.selections(), chnx, chny, color_channel, color_settings)

        if fit_settings['ellipse2']:
            self.draw_ellipse2(*values, color_channel, color_settings)

        self.canvas.setxy(values[:2])

        sel = data.activeSelection()
        if sel is None:
            self.label.setText('no active selection')
        else:
            _,_,age_msg = self.time_result_for_selection(sel)
            self.label.setText(f'active selection name: {sel.name}, {age_msg}')
            if fit_settings['errorbar']:
                self.draw_errorbar(*self.extract_values([sel], chnx, chny))
            if fit_settings['ellipse']:
                self.draw_ellipse(
                    [sel], chnx, chny, edgecolor='orange', linewidth=3)
            if fit_settings['ellipse2']:
                values_sel = self.extract_values([sel], chnx, chny)
                edgecolor = self.get_color_for_selection(sel, color_channel, color_settings)
                ellipse_from_errorbar(
                    values_sel[0][0], values_sel[1][0], 
                    values_sel[2][0], values_sel[3][0],
                    self.canvas.axes,
                    edgecolor=edgecolor, facecolor=edgecolor, alpha=0.3, linewidth=2)

        self.set_limits()
        self.draw_ref_lines(k, b, values[0])
        self.draw_err_age_band(k, b, values[0])

##        print("调整布局前...")
##        self.print_figure_info("调整布局前")
        self.canvas.adjust_layout()
##        print("调整布局后...")
##        self.print_figure_info("调整布局后")

        self.canvas.draw()
        self.update_rho_label(grp, chnx, chny)
##        print("=== show_selection结束 ===\n")


#### rho (weightedMean2D)
    def update_rho_label(self, grp, chnx, chny):
        """计算并显示weightedMean2D结果"""
        try:
            selections = grp.selections()
            if len(selections) < 2:
                self.rho_label.setText(
                    f'weightedMean2D: n={len(selections)}, '
                    f'need at least 2 selections')
                return

            x_list, y_list, sx_list, sy_list, pearson_list = \
                [], [], [], [], []
            for sel in selections:
                x_list.append(data.result(sel, chnx).value())
                y_list.append(data.result(sel, chny).value())
                sx_list.append(
                    data.result(sel, chnx).uncertaintyAs2SE() / 2)
                sy_list.append(
                    data.result(sel, chny).uncertaintyAs2SE() / 2)
                signalx = chnx.dataForSelection(sel)
                signaly = chny.dataForSelection(sel)
                mask = np.isfinite(signalx) & np.isfinite(signaly)
                if np.sum(mask) > 1:
                    pearson_list.append(
                        np.corrcoef(signalx[mask], signaly[mask])[0, 1])
                else:
                    pearson_list.append(0.0)

            twm = weightedMean2D(
                np.array(x_list), np.array(sx_list),
                np.array(y_list), np.array(sy_list),
                np.array(pearson_list))

            line1 = (
                f"x_bar={twm['x_bar']:.4g} ± {twm['sigma_x_bar']:.4g}, "
                f"y_bar={twm['y_bar']:.4g} ± {twm['sigma_y_bar']:.4g}")
            line2 = (
                f"cov_xy_bar={twm['cov_xy_bar']:.4g}, "
                f"rho_xy_bar={twm['rho_xy_bar']:.4g}, "
                f"MSWD={twm['mswd']:.4g}, "
                f"prob={twm['prob']:.4g}, "
                f"n={twm['n']}")
            self.rho_label.setText(line1 + '\n' + line2)
        except Exception as e:
            self.rho_label.setText(f'weightedMean2D error: {e}')

#### colors
    def initialize_color_range(self, selections, color_channel):
        print('starting initialize_color_range')
        ppm_values = []
        for sel in selections:
            try:
                value = data.result(sel, color_channel).value()
                ppm_values.append(value)
            except:
                continue

        print('ppm_values', ppm_values)
        if len(ppm_values) > 0:
            self.color_min = min(ppm_values)
            self.color_max = max(ppm_values)
            if abs(self.color_max - self.color_min) < 1e-10:
                self.color_max = self.color_min + 1.0
            self.color_range_initialized = True
            
            self.color_setting_w.set_range_label(
                f"{self.color_min:.2f} - {self.color_max:.2f}")

            self.update_colorbar()

##    def update_colorbar_old(self):
##        print("\n=== update_colorbar ===")
##        print(f"当前colorbar: {self.cbar}")
##        print(f"颜色范围初始化: {self.color_range_initialized}")
##
##        if self.cbar is not None:
##            print(f"移除现有colorbar: {self.cbar}")
##            try:
##                if hasattr(self.cbar, 'ax'):
##                    print(f"Colorbar axes: {self.cbar.ax}")
##                self.cbar.remove()
##                print("成功移除colorbar")
##            except:
##                print(f"移除colorbar时出错: {repr(e)}")
##
##            self.cbar = None
##
##        self.canvas.draw_idle()
##        print("移除colorbar后刷新画布")
##
##        print("创建新的colorbar...")
##        norm = matplotlib.colors.Normalize(vmin=self.color_min, vmax=self.color_max)
##        sm = plt.cm.ScalarMappable(cmap=COLORMAP, norm=norm)
##        sm.set_array([])
##        
##        self.cbar = self.canvas.figure.colorbar(sm, ax=self.canvas.axes)
##        self.cbar.set_label('Concentration (ppm)')
##
##        print(f"新colorbar创建: {self.cbar}")
##        print(f"新colorbar axes: {self.cbar.ax}")
##        print(f"新colorbar axes位置: {self.cbar.ax.get_position()}")
##
##        self.canvas.draw_idle()
##        print("创建colorbar后刷新画布")
##
##        self.canvas.figure.subplots_adjust(right=0.85)
##        print("设置了subplots_adjust(right=0.85)")

    def update_colorbar(self):
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
                        cmap=COLORMAP, 
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
        if self.cbar is not None:
            try:
                self.cbar.remove()
            except:
                pass
        
        norm = matplotlib.colors.Normalize(
            vmin=self.color_min, 
            vmax=self.color_max
        )
        sm = plt.cm.ScalarMappable(cmap=COLORMAP, norm=norm)
        sm.set_array([])
        
        self.cbar = self.canvas.figure.colorbar(sm, ax=self.canvas.axes)
        self.cbar.set_label('Concentration (ppm)')


    def get_color_for_selection(self, selection, color_channel, color_settings):
        if not color_channel:
            return 'orange'
        
        try:
            value = data.result(selection, color_channel).value()
            
            norm_value = np.clip(value, self.color_min, self.color_max)
            normalized = (norm_value - self.color_min) / (self.color_max - self.color_min)
            
            return COLORMAP(normalized)
        except:
            return 'orange'

    def draw_colored_ellipses(self, selections, chnx, chny, color_channel, color_settings):
        for sel in selections:
            edgecolor = self.get_color_for_selection(sel, color_channel, color_settings)
            confidence_ellipse(
                sel, chnx, chny, 2, self.canvas.axes,
                edgecolor=edgecolor, facecolor=edgecolor, alpha=0.3, linewidth=2)

    def draw_ellipse2(self, x, y, xerr, yerr, color_channel, color_settings):
        grp = data.activeSelectionGroup()
        if not grp:
            return
        
        selections = grp.selections()
        
        for i, sel in enumerate(selections):
            edgecolor = self.get_color_for_selection(sel, color_channel, color_settings)
            ellipse_from_errorbar(
                x[i], y[i], xerr[i], yerr[i], self.canvas.axes,
                edgecolor=edgecolor, facecolor=edgecolor, alpha=0.3, linewidth=2)

    def reset_color_range(self):
##        print("\n=== reset_color_range ===")
##        print(f"当前颜色范围: {self.color_min} - {self.color_max}")
##        print(f"当前colorbar: {self.cbar}")

        self.color_range_initialized = False
        self.show_selection()
####

    def time_result_for_selection(self, selection):
        if 'time' in data.timeSeriesNames():
            time_chn = data.timeSeries('time')
            iolite_result = data.result(selection, time_chn)
            age = iolite_result.value()
            age_err = iolite_result.uncertaintyAs2SE()
            age_msg = f'{age:.2f}±{age_err:.2f}Ma'
            print(age_msg)
            return age, age_err, age_msg
        else:
            return float('nan'), float('nan'), '-'

    def extract_values(self, selections, chnx, chny):
        x = np.array([data.result(sel, chnx).value() for sel in selections])
        y = np.array([data.result(sel, chny).value() for sel in selections])
        xerr = [data.result(sel, chnx).uncertaintyAs2SE() for sel in selections]
        yerr = [data.result(sel, chny).uncertaintyAs2SE() for sel in selections]
        return x,y,xerr,yerr

    def draw_errorbar(self, x,y,xerr,yerr):
        self.canvas.axes.errorbar(
            x,y,xerr=xerr, yerr=yerr, fmt='o', capsize=3)

    def draw_ellipse(self, selections, chnx, chny, edgecolor='red', **kwargs):
        for sel in selections:
            confidence_ellipse(
                sel, chnx, chny, 2, self.canvas.axes,
                edgecolor=edgecolor, **kwargs)

    def set_limits(self):
        settings = self.fit_w.getSettings()
        xmin = settings['xmin']
        xmax = settings['xmax']
        self.canvas.axes.set_xlim(left=xmin, right=xmax)
        ymin = settings['ymin']
        ymax = settings['ymax']
        self.canvas.axes.set_ylim(bottom=ymin, top=ymax)

    def y_intercept_title(self, inverse=False):
        xy_settings = self.xyplot_setting_w.getSettings()
        print('xy_settings', xy_settings)
        try:
            yname = xy_settings['Ychannel']
            print(yname)
        except KeyError: return
        try:
            m = re.search(
                r'([A-z]{1,2})([0-9]{1,3})_([A-z]{1,2})([0-9]{1,3})',
                yname)
        except TypeError: return
        try: elem1, n1, elem2, n2 = m.groups()
        except Exception as e:
            print('y_intercept_title', yname, e)
            return 'y$_0$'
        if not inverse:
            title = r'$^{' + n1 + '/' + n2 + '}$' + elem1 + '$_0$'
        else:
            title = r'$^{' + n2 + '/' + n1 + '}$' + elem1 + '$_0$'
        return title


    def calculate_mswd(self, x, y, k, b, xerr, yerr):
        """
        计算加权拟合的MSWD
        使用x和y的误差进行权重计算
        """
        if len(x) < 3:  # 至少需要3个点才能计算有意义的MSWD
            return float('nan')
        
        # 计算每个点的残差
        y_pred = k * x + b
        residuals = y - y_pred
        
        # 计算每个点的权重（使用误差的倒数平方）
        # 这里组合x和y的误差，通过误差传播得到y的预测误差
        weights = []
        for i in range(len(x)):
            # 预测值y_pred的误差：考虑x误差通过斜率传播 + y误差
            sigma_y_pred = np.sqrt((k * xerr[i])**2 + yerr[i]**2)
            if sigma_y_pred > 0:
                weight = 1 / (sigma_y_pred**2)
            else:
                weight = 0
            weights.append(weight)
        
        weights = np.array(weights)
        
        # 计算加权残差平方和
        weighted_residuals = residuals * np.sqrt(weights)
        sum_weighted_residuals_sq = np.sum(weighted_residuals**2)
        
        # 自由度 = n - 2 (拟合了斜率和截距两个参数)
        dof = len(x) - 2
        
        if dof <= 0:
            return float('nan')
        
        mswd = sum_weighted_residuals_sq / dof
        
        return mswd

    def calculate_mswd_with_different_weights(self, x, y, k, b, xerr, yerr, cov_xy=None):
        """
        计算三种不同权重假设下的MSWD
        
        参数:
        - xerr, yerr: 1σ误差
        - cov_xy: 协方差（如果存在）
        
        返回: (mswd_separate, mswd_combined, mswd_full) 三种MSWD
        """
        if len(x) < 3:
            return float('nan'), float('nan'), float('nan')
        
        y_pred = k * x + b
        residuals = y - y_pred
        n = len(x)
        dof = n - 2
        
        if dof <= 0:
            return float('nan'), float('nan'), float('nan')
        
        # 1. 只考虑Y误差（传统方案）
        weights_y = []
        for i in range(n):
            if yerr[i] > 0:
                weights_y.append(1 / yerr[i]**2)
            else:
                weights_y.append(0)
        weights_y = np.array(weights_y)
        mswd_y = np.sum((residuals * np.sqrt(weights_y))**2) / dof
        
        # 2. 组合误差（X误差通过斜率传播 + Y误差）
        weights_combined = []
        for i in range(n):
            sigma_combined = np.sqrt((k * xerr[i])**2 + yerr[i]**2)
            if sigma_combined > 0:
                weights_combined.append(1 / sigma_combined**2)
            else:
                weights_combined.append(0)
        weights_combined = np.array(weights_combined)
        mswd_combined = np.sum((residuals * np.sqrt(weights_combined))**2) / dof
        
        # 3. 完全协方差方案（使用误差椭圆信息）
        if cov_xy is not None:
            weights_full = []
            for i in range(n):
                # 构建协方差矩阵
                cov_matrix = np.array([[xerr[i]**2, cov_xy[i]], 
                                      [cov_xy[i], yerr[i]**2]])
                
                # 计算预测值的总方差，考虑X和Y的协方差
                # 通过误差传播：Var(y_pred) = k²·Var(x) + Var(y) + 2k·Cov(x,y)
                var_pred = (k**2 * xerr[i]**2 + yerr[i]**2 + 
                           2 * k * cov_xy[i])
                
                if var_pred > 0:
                    weights_full.append(1 / var_pred)
                else:
                    weights_full.append(0)
            
            weights_full = np.array(weights_full)
            mswd_full = np.sum((residuals * np.sqrt(weights_full))**2) / dof
        else:
            mswd_full = float('nan')
        
        return mswd_y, mswd_combined, mswd_full

    def draw_linear_fit(self, x, y):
        settings = self.fit_w.getSettings()
        if settings['fixed']:
            y_int = settings['y_intercept']
            k, b, r, k_err = self.linear_fit_boundary(x, y, y_int)
            b_err = 0
        else:
            k, b, r, k_err, b_err = self.linear_fit_no_boundary(x, y)

        # 计算MSWD
        xerr = []
        yerr = []
        cov_xy_list = []
        grp = data.activeSelectionGroup()

        if grp:
            selections = grp.selections()
            chnx = data.timeSeries(self.xyplot_setting_w.getSettings()['Xchannel'])
            chny = data.timeSeries(self.xyplot_setting_w.getSettings()['Ychannel'])
            for sel in selections:
                xerr.append(data.result(sel, chnx).uncertaintyAs2SE() / 2)  # 2SE转1SE
                yerr.append(data.result(sel, chny).uncertaintyAs2SE() / 2)

                x_data = chnx.dataForSelection(sel)
                y_data = chny.dataForSelection(sel)
                if len(x_data) > 1 and len(y_data) > 1:
                    cov = np.cov(x_data, y_data)[0, 1]
                    cov_xy_list.append(cov)
                else:
                    cov_xy_list.append(0)
##        mswd = self.calculate_mswd(x, y, k, b, np.array(xerr), np.array(yerr))
        mswd_y, mswd_combined, mswd_full = self.calculate_mswd_with_different_weights(
            x, y, k, b, np.array(xerr), np.array(yerr), np.array(cov_xy_list))

        fit_settings = self.fit_w.getSettings()
        if fit_settings['ellipse']:  # 带协方差的椭圆
            mswd_display = mswd_full
            mswd_type = "MSWD_full"
        elif fit_settings['ellipse2']:  # 不带协方差的椭圆
            mswd_display = mswd_combined
            mswd_type = "MSWD_comb"
        else:  # 只有errorbar
            mswd_display = mswd_y
            mswd_type = "MSWD_y"


        xmin = settings['xmin']
        xmax = settings['xmax']
        x1, x2 = xmin, xmax
        if x1 is None: x1 = x.min()
        if x2 is None: x2 = x.max()
        X = np.array([x1, x2])
        Y = X*k+b
        self.canvas.axes.plot(X, X*k+b, color=(0.7,0.7,0.7))
        x_intercept = -b/k

        method = self.xyplot_setting_w.getSettings()['age_method']
        t, t_err = float('nan'), float('nan')

        system = METHOD_TO_SYSTEM.get(method)
        if system and system in SYSTEM_LAMBDAS:
            Lambda = SYSTEM_LAMBDAS[system]
        else:
            Lambda = None

        try:
            if method == 'inversHf' and Lambda:
                t, t_err = self.calculate_invers(k, b, k_err, b_err, Lambda)
            elif method == 'isochHf' and Lambda:
                t, t_err = self.calculate_isoch(k, k_err, Lambda)
            elif method == 'Pb76':
                t = agePb76(b)
            elif method == 'isochSr' and Lambda:
                t, t_err = self.calculate_isoch(k, k_err, Lambda)
            elif method == 'inversSr' and Lambda:
                t, t_err = self.calculate_invers(k, b, k_err, b_err, Lambda)
            elif method == 'isochKCa' and Lambda:
                t, t_err = self.calculate_isoch(k, k_err, Lambda)
            elif method == 'inversKCa' and Lambda:
                t, t_err = self.calculate_invers(k, b, k_err, b_err, Lambda)
            elif method == 'isochMoDau' and Lambda:
                t, t_err = self.calculate_isoch(k, k_err, Lambda)
            elif method == 'inversMoDau' and Lambda:
                t, t_err = self.calculate_invers(k, b, k_err, b_err, Lambda)
            elif method == 'isochBa' and Lambda:
                t, t_err = self.calculate_isoch(k, k_err, Lambda)
            elif method == 'inversBa' and Lambda:
                t, t_err = self.calculate_invers(k, b, k_err, b_err, Lambda)
            else:
                t = float('nan')
        except Exception: t = float('nan')

        age = f'{t:0.2f}' if t else f'{t}'
        age_err = f'± {t_err:0.2f}' if t_err else ''
        age_msg = f'{age}' + age_err
##        mswd_str = f'MSWD = {mswd:.2f}' if not math.isnan(mswd) else ''
        mswd_str = f'{mswd_type} = {mswd_display:.2f}' if not math.isnan(mswd_display) else ''

        b_title = self.y_intercept_title('invers' in method)
        if b_title is None: b_title = 'b'
        if 'invers' in method: b_value = 1/b
        else: b_value = b

        self.canvas.axes.set_title(
            f'k={k:0.4f}, b={b:0.4f}, a={(-b/k):0.2f}, R^2={r**2:0.2f}\n'+
            b_title
            + f':{b_value:0.4f}, age_{method}/Ma:{age_msg}, {mswd_str}')

        return k, b

    def calculate_isoch(self, k, k_err, Lambda):
        t = t = math.log(k+1)/Lambda / 1e6
        t_err = (k_err / (k + 1)) / Lambda / 1e6
        return t, t_err

    def calculate_invers(self, k, b, k_err, b_err, Lambda):
        t = t = math.log(1 + 1 / (-b / k)) / Lambda / 1e6
        dt_dk = (1 / (k * (b + k))) / Lambda / 1e6
        dt_db = (-1 / (b * (b + k))) / Lambda / 1e6
        t_err = np.sqrt((dt_dk * k_err)**2 + (dt_db * b_err)**2)
        return t, t_err

    def linear_fit_boundary(self, x, y, y_int):
        def func(x, k): return k * x + y_int
        popt, pcov = optimize.curve_fit(func, x, y)
        k = popt[0]
        k_err = np.sqrt(pcov[0][0])
        corr_coef = np.corrcoef(x, y)[0][1]
        sqrt_N = x.shape[0] ** 0.5
        return k, y_int, corr_coef**2, k_err / sqrt_N

    def linear_fit_no_boundary(self, x, y):
        res = stats.linregress(x, y)
        k, b, r, k_err, b_err = res.slope, res.intercept, res.rvalue, res.stderr, res.intercept_stderr
        sqrt_N = x.shape[0] ** 0.5
        return k, b, r, k_err / sqrt_N, b_err / sqrt_N

    def draw_ref_lines(self, k, b, x):
        """画参考等时线 (Reference isochron lines)。

        反等时线与正等时线均固定 y 截距为当前拟合 b 值，按年龄计算斜率：
          反等时线 (method 含 'invers')：slope = -b * (exp(λ·t) - 1)  (负)
          正等时线 (method 含 'isoch') ：slope =       (exp(λ·t) - 1)  (正，与 b 无关)
        每条线横跨整个 x 范围 (xmin 到 xmax)，在线的可见末端标注年龄。
        年龄项可以是绝对值 (Ma)，也可以是相对百分数 (如 +5%、-10%)，
            后者按当前拟合年龄换算，标注为 "+5% (945 Ma)" 形式，
            并用点线 (:) 与绝对年龄的虚线 (--) 区分。
        标注位置随斜率符号自适应：负斜率线穿出底边时标注在上方，
            正斜率线穿出顶边时标注在下方，避免飘出绘图区。
        线形：淡蓝色；标注 <1 Ga 用 Ma，>=1 Ga（含）用 Ga。
        必须在 set_limits() 之后调用，此时 xlim/ylim 已确定，画线自动裁剪。
        """
        fit_settings = self.fit_w.getSettings()
        if not fit_settings.get('show_ref_lines'):
            return
        ages = fit_settings.get('ref_ages') or []
        if not ages:
            return

        method = self.xyplot_setting_w.getSettings().get('age_method')
        system = METHOD_TO_SYSTEM.get(method)
        Lambda = SYSTEM_LAMBDAS.get(system) if system else None
        if Lambda is None:
            return

        is_inverse = 'invers' in method
        if not is_inverse and 'isoch' not in method:
            return

        xmin = fit_settings['xmin']
        xmax = fit_settings['xmax']
        x1 = x.min() if xmin is None else xmin
        x2 = x.max() if xmax is None else xmax
        if x2 <= x1:
            return

        style = dict(linestyle='--', color=(0.4, 0.6, 0.8),
                     linewidth=1.0, zorder=1)
        ax = self.canvas.axes

        # 当前拟合年龄（相对百分数项的换算基准）
        t_fit = None
        if any(str(tok).endswith('%') for tok in ages):
            try:
                if is_inverse:
                    # 与 calculate_invers 一致: t = ln(1 + 1/(-b/k)) / λ
                    t_fit = math.log(1 + 1 / (-b / k)) / Lambda / 1e6
                else:
                    # 与 calculate_isoch 一致: t = ln(k+1) / λ
                    t_fit = math.log(k + 1) / Lambda / 1e6
            except (ValueError, ZeroDivisionError):
                t_fit = None

        for tok in ages:
            tok = str(tok)
            if tok.endswith('%'):
                if t_fit is None or not (t_fit > 0):
                    continue
                pct = float(tok[:-1])
                t_ma = t_fit * (1 + pct / 100.0)
                label = f'{tok} ({format_ref_age(t_ma)})'
                ls = ':'  # 相对参考线用点线与绝对线区分
            else:
                t_ma = float(tok)
                label = format_ref_age(t_ma)
                ls = '--'
            if not (t_ma > 0) or t_ma == float('inf'):
                continue
            t = t_ma * 1e6  # Ma -> years
            if is_inverse:
                slope = -b * (math.exp(Lambda * t) - 1)
            else:
                slope = math.exp(Lambda * t) - 1
            ax.plot([x1, x2], [slope*x1 + b, slope*x2 + b],
                    linestyle=ls, color=style['color'],
                    linewidth=style['linewidth'], zorder=style['zorder'])
            x_end, y_end = self._visible_line_end(slope, b, x1, x2)
            y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
            # 负斜率线穿出底边 -> 标注放线上方；正斜率线穿出顶边 -> 放线下方
            if slope < 0:
                va, y_text = 'bottom', y_end + 0.015 * y_range
            else:
                va, y_text = 'top', y_end - 0.015 * y_range
            ax.text(
                x_end, y_text, label,
                color=style['color'], fontsize=7,
                ha='right', va=va)

    def _visible_line_end(self, k, b, x1, x2):
        """返回直线 y=k*x+b 在 x∈[x1,x2] 内可见部分的右端坐标。

        负斜率（反等时线）若在到达右边界前已穿出 y 轴底边，则返回
        底边交点坐标，避免年龄标注飘到绘图区外；正斜率对称处理顶边。
        """
        ax = self.canvas.axes
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        xl = max(x1, xlim[0])
        xr = min(x2, xlim[1])
        x_end = xr
        if k < 0:
            if k*xr + b < ylim[0]:  # 右端在可见区下方 -> 线穿出底边
                x_cross = (ylim[0] - b) / k
                x_end = max(xl, min(xr, x_cross))
        elif k > 0:
            if k*xr + b > ylim[1]:  # 右端在可见区上方 -> 线穿出顶边
                x_cross = (ylim[1] - b) / k
                x_end = max(xl, min(xr, x_cross))
        y_end = k*x_end + b
        return x_end, y_end

    def draw_err_age_band(self, k, b, x):
        """画年龄误差带：在拟合年龄 ± 百分数处画两条参考线，中间填灰色阴影。

        百分数（如 10%, -10%）以当前拟合年龄为基准换算成实际年龄，
        再按等时线公式算斜率。两条线分别对应 pct_min 和 pct_max，
        用 fill_between 在它们之间画半透明灰色区域。
        必须在 set_limits() 之后调用。
        """
        fit_settings = self.fit_w.getSettings()
        if not fit_settings.get('show_err_age'):
            return
        err_tokens = fit_settings.get('err_age') or []
        if len(err_tokens) < 2:
            return

        method = self.xyplot_setting_w.getSettings().get('age_method')
        system = METHOD_TO_SYSTEM.get(method)
        Lambda = SYSTEM_LAMBDAS.get(system) if system else None
        if Lambda is None:
            return

        is_inverse = 'invers' in method
        if not is_inverse and 'isoch' not in method:
            return

        # 当前拟合年龄
        try:
            if is_inverse:
                t_fit = math.log(1 + 1 / (-b / k)) / Lambda / 1e6
            else:
                t_fit = math.log(k + 1) / Lambda / 1e6
        except (ValueError, ZeroDivisionError):
            return
        if not (t_fit > 0):
            return

        # 解析百分数
        pcts = []
        for tok in err_tokens:
            tok = str(tok)
            if not tok.endswith('%'):
                continue
            try:
                pcts.append(float(tok[:-1]))
            except ValueError:
                continue
        if len(pcts) < 2:
            return
        pct_lo, pct_hi = min(pcts), max(pcts)

        xmin = fit_settings['xmin']
        xmax = fit_settings['xmax']
        x1 = x.min() if xmin is None else xmin
        x2 = x.max() if xmax is None else xmax
        if x2 <= x1:
            return

        ax = self.canvas.axes
        xs = np.linspace(x1, x2, 100)

        def slope_for_pct(pct):
            t_ma = t_fit * (1 + pct / 100.0)
            t = t_ma * 1e6
            if is_inverse:
                return -b * (math.exp(Lambda * t) - 1)
            else:
                return math.exp(Lambda * t) - 1

        k_lo = slope_for_pct(pct_lo)
        k_hi = slope_for_pct(pct_hi)
        y_lo = k_lo * xs + b
        y_hi = k_hi * xs + b

        # 灰色阴影
        ax.fill_between(xs, y_lo, y_hi, facecolor='0.85',
                         edgecolor='none', alpha=0.5, zorder=0.5)
        # 两条边界线
        line_style = dict(linestyle='--', color=(0.4, 0.6, 0.8),
                          linewidth=1.0, zorder=1)
        ax.plot([x1, x2], [k_lo*x1+b, k_lo*x2+b], **line_style)
        ax.plot([x1, x2], [k_hi*x1+b, k_hi*x2+b], **line_style)

    def disp_sel_name(self, idx):
        grp = data.activeSelectionGroup()
        if not grp: return
        selection = grp.selections()[idx]
        name_msg = selection.name
        _,_,age_msg = self.time_result_for_selection(selection)
        self.label.setText(f'{name_msg}, {age_msg}')

class PlotCanvas(FigureCanvasQTAgg):
    got_argmin = Signal(int)
    def __init__(self, parent=None, xlim=None, ylim=None,
                 width=4, height=1.6, dpi=80):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        self.xlim, self.ylim = xlim, ylim
        FigureCanvasQTAgg.__init__(self, fig)
        self.setParent(parent)
        FigureCanvasQTAgg.updateGeometry(self)
        fig.canvas.mpl_connect('motion_notify_event', self.on_hover)

    def clear_figure(self):
        self.axes.cla()

    def clear_all_axes(self):
        for ax in self.figure.axes[1:]:
            ax.remove()
        self.figure.axes.clear()
        self.axes = None

    def set_ax_label(self, xlabel='', ylabel=''):
        self.axes.set_xlabel(xlabel)
        self.axes.set_ylabel(ylabel)

    def on_hover(self, event):
        x,y = event.xdata, event.ydata
        if x is None or y is None: return
        if not hasattr(self, 'XY'): return
        dist = np.sum((self.XY - np.array([x,y]))**2, axis=1)
        i = np.argmin(dist)
        self.got_argmin.emit(i)

    def setxy(self, xy):
        self.XY = np.array(xy).T

    def adjust_layout(self):
        self.figure.subplots_adjust(left=0.15, right=0.85, bottom=0.15, top=0.8)
        self.draw_idle()

class XYPlotSettingWidget(QWidget):
    channel_changed = Signal()
    got_y_interp = Signal(float)

    def __init__(self):
        QWidget.__init__(self)
        self.settings = {}
        formLayout = QtGui.QFormLayout()
        formLayout.setContentsMargins(0,0,0,0)
        self.setLayout(formLayout)
        self.all_cb = []

        channels = data.timeSeriesNames()
        self.x_cb = QtGui.QComboBox(self)
        self.y_cb = QtGui.QComboBox(self)
        self.all_cb.extend([self.x_cb, self.y_cb])

        age_methods = METHODS
        self.age_cb = QtGui.QComboBox(self)
        self.age_cb.addItems(age_methods)
        self.settings['age_method'] = age_methods[0]
        self.age_cb.setCurrentText(age_methods[0])
        self.age_cb.currentTextChanged.connect(partial(self.setSetting,'age_method'))
        self.age_cb.currentTextChanged.connect(self.update_channels_for_method)
        self.age_cb.currentTextChanged.connect(self.emit_channel_changed)
        formLayout.addRow('age method', self.age_cb)

        for cb, label in zip([self.x_cb, self.y_cb], ['Xchannel', 'Ychannel']):
            cb.addItems(channels)
            if channels:
                cb.setCurrentText(channels[0])
                self.setSetting(label, cb.currentText)
            cb.currentTextChanged.connect(partial(self.setSetting, label))
            cb.currentTextChanged.connect(self.emit_channel_changed)
            formLayout.addRow(f'{label} name', cb)

        # Lambda 输入行
        lambda_layout = QtGui.QHBoxLayout()
        self.lambda_label = QLabel("λ (1/y):")
        self.lambda_edit = QtGui.QLineEdit()
        self.order_label = QLabel('*1e-11')
##        self.lambda_edit.setEnabled(False)  # 默认禁用
        self.lambda_edit.returnPressed.connect(self.on_lambda_changed)
        lambda_layout.addWidget(self.lambda_label)
        lambda_layout.addWidget(self.lambda_edit)
        lambda_layout.addWidget(self.order_label)
        lambda_layout.addStretch()
        formLayout.addRow('Decay constant', lambda_layout)

        self.update_lambda_display(self.age_cb.currentText)

    def update_lambda_display(self, method):
        """根据当前方法更新lambda输入框的显示"""
        system = METHOD_TO_SYSTEM.get(method)
        if system and system in SYSTEM_LAMBDAS and SYSTEM_LAMBDAS[system] is not None:
            self.lambda_edit.setText(f"{SYSTEM_LAMBDAS[system]*1e11:.4f}")
            self.lambda_edit.setToolTip(f"Return to change λ")

    def on_lambda_changed(self):
        """用户编辑lambda值后更新全局字典"""
        method = self.age_cb.currentText
        system = METHOD_TO_SYSTEM.get(method)
        if not system:
            return
        
        try:
            new_lambda = float(self.lambda_edit.text)*1e-11
            old_lambda = SYSTEM_LAMBDAS.get(system)
            
            if old_lambda != new_lambda:
                SYSTEM_LAMBDAS[system] = new_lambda
                if old_lambda is None:
                    old_lambda = float('nan')
                print(f"Updated λ for {system}: {old_lambda:.4e} -> {new_lambda:.4e}")
                # 触发重绘以更新年龄计算
                self.channel_changed.emit()
        except ValueError:
            # 恢复原值
            if system in SYSTEM_LAMBDAS:
                self.lambda_edit.setText(f"{SYSTEM_LAMBDAS[system]*1e11:.4f}")
            QMessageBox.warning(self, "Invalid Input", 
                               "Please enter a valid number in scientific notation (e.g., 1.865e-11)")

        print(SYSTEM_LAMBDAS)

    def update_channels_for_method(self, method):
        self.update_lambda_display(method)

        method_channels = {
            'isochSr': ('final_Rb87_Sr86', 'final_Sr87_Sr86', 0.710),
            'inversSr': ('final_Rb87_Sr87', 'final_Sr86_Sr87', 1.408),
            'inversHf': ('final_Lu176_Hf176', 'final_Hf177_Hf176', 3.55),
            'isochHf': ('final_Lu176_Hf177', 'final_Hf176_Hf177', 0.282),
            'inversKCa': ('final_K40_Ca40', 'final_Ca44_Ca40', 0.0215),
            'isochKCa': ('final_K40_Ca44', 'final_Ca40_Ca44', 46.5),
            'Pb76': ('final_Pb204_Pb206', 'final_Pb207_Pb206', float('nan')),
            'isochMoDau': ('Mother_Nonradio', 'Daughter_Nonradio', float('nan')),
            'inversMoDau': ('Mother_Daughter', 'Nonradio_Daughter', float('nan')),
            'inversBa': ('final_La138_Ba138', 'final_Ba137_Ba138', 0.1565),
            'isochBa': ('final_La138_Ba137', 'final_Ba138_Ba137', 6.388),
        }
        
        if method not in method_channels:
            return
            
        base_x, base_y, y_interp = method_channels[method]
        channels = data.timeSeriesNames()
        
        def get_variants(base):
            variants = []
            pattern = re.compile(rf'^{base}(?:_[a-zA-Z]+)*$')
            for ch in channels:
                if pattern.match(ch):
                    variants.append(ch)
            return variants
            
        x_variants = get_variants(base_x)
        y_variants = get_variants(base_y)
        
        def update_combobox(cb, variants):
            current = cb.currentText
            cb.blockSignals(True)
            cb.clear()
            
            for variant in variants:
                cb.addItem(variant)
            
            if variants:
                cb.insertSeparator(len(variants))
            
            for ch in channels:
                if ch not in variants:
                    cb.addItem(ch)
            
            if variants:
                cb.setCurrentText(variants[0])
            cb.blockSignals(False)
        
        update_combobox(self.x_cb, x_variants)
        update_combobox(self.y_cb, y_variants)
        
        if x_variants:
            self.setSetting('Xchannel', self.x_cb.currentText)
        if y_variants:
            self.setSetting('Ychannel', self.y_cb.currentText)
        self.got_y_interp.emit(y_interp)

    def getSettings(self):
        return self.settings
    def setSetting(self, key, value):
        self.settings[key] = value
    def emit_channel_changed(self, text):
        self.channel_changed.emit()

    def update_cb(self):
        channels = data.timeSeriesNames()
        for cb in self.all_cb:
            cb.blockSignals(True)
            cb.clear()
            cb.addItems(channels)
            cb.blockSignals(False)
        self.update_channels_for_method(self.age_cb.currentText)

class FittingOptionWidget(QWidget):
    setting_changed = Signal()
    def __init__(self):
        QWidget.__init__(self)
        self.settings = {
            'fixed': True, 'y_intercept': 3.55,
            'xmin_enabled': False, 'xmin': 0,
            'xmax_enabled': False, 'xmax': float('inf'),
            'ymin_enabled': False, 'ymin': 0,
            'ymax_enabled': False, 'ymax': float('inf'),
            'errorbar': True, 'ellipse': False, 'ellipse2': False,
            'show_ref_lines': False,
            'ref_ages': ['100', '200', '500', '1000', '2000', '3000',
                         '4000', '4500'],
            'show_err_age': False,
            'err_age': ['10%', '-10%'],
            }

        hbox1 = QtGui.QHBoxLayout()
        self.fixed_w = QtGui.QCheckBox()
        self.fixed_w.setChecked(self.settings['fixed'])
        self.fixed_w.toggled.connect(lambda t: self.setSetting("fixed", bool(t)))
        self.y_intercept = QtGui.QLineEdit()
        self.y_intercept.setText(str(self.settings['y_intercept']))
        self.y_intercept.returnPressed.connect(
            partial(self.on_return, 'y_intercept', self.y_intercept))
        hbox1.addWidget(QLabel('y_intercept'))
        hbox1.addWidget(self.fixed_w)
        hbox1.addWidget(self.y_intercept)

        vbox = QtGui.QVBoxLayout()
        vbox.setContentsMargins(0,0,0,0)
        vbox.addLayout(hbox1)
        self.edits = {}
        hboxes = {}
        for name in ['xmin', 'xmax', 'ymin', 'ymax']:
            hbox = QtGui.QHBoxLayout()
            cb = QtGui.QCheckBox()
            cb.setChecked(self.settings[f'{name}_enabled'])
            cb.toggled.connect(partial(self.setBool, f"{name}_enabled"))
            edit = QtGui.QLineEdit()
            self.edits[name] = edit
            edit.setText(str(self.settings[name]))
            edit.textChanged.connect(partial(self.setFloat, name))
            hbox.addWidget(QLabel(name))
            hbox.addWidget(cb)
            hbox.addWidget(edit)
            hboxes[name] = hbox
        gbox = QtGui.QGridLayout()
        gbox.addLayout(hboxes['xmin'], 0,0)
        gbox.addLayout(hboxes['xmax'], 0,1)
        gbox.addLayout(hboxes['ymin'], 1,0)
        gbox.addLayout(hboxes['ymax'], 1,1)
        vbox.addLayout(gbox)

        plot_option_box = QtGui.QHBoxLayout()
        display_text = ['errorbar', 'ellipse with Rho', 'ellipse without Rho']
        for plot_style, cb_text in zip(
            ['errorbar', 'ellipse', 'ellipse2'], display_text):
            cb = QtGui.QCheckBox(cb_text)
            cb.setChecked(self.settings[plot_style])
            cb.toggled.connect(partial(self.setBool, plot_style))
            plot_option_box.addWidget(cb)
        vbox.addLayout(plot_option_box)

        ref_box = QtGui.QHBoxLayout()
        self.show_ref_w = QtGui.QCheckBox('Show ref lines')
        self.show_ref_w.setChecked(self.settings['show_ref_lines'])
        self.show_ref_w.toggled.connect(
            partial(self.setBool, 'show_ref_lines'))
        self.ref_ages_edit = QtGui.QLineEdit()
        self.ref_ages_edit.setText(', '.join(self.settings['ref_ages']))
        self.ref_ages_edit.setToolTip(
            'Comma-separated ages (Ma) and/or relative offsets like +5%, '
            '-10% (applied to the fitted age), press Enter to apply')
        self.ref_ages_edit.returnPressed.connect(self.on_ref_ages_changed)
        ref_box.addWidget(self.show_ref_w)
        ref_box.addWidget(self.ref_ages_edit)
        vbox.addLayout(ref_box)

        err_box = QtGui.QHBoxLayout()
        self.show_err_age_w = QtGui.QCheckBox('Show error of age')
        self.show_err_age_w.setChecked(self.settings['show_err_age'])
        self.show_err_age_w.toggled.connect(
            partial(self.setBool, 'show_err_age'))
        self.err_age_edit = QtGui.QLineEdit()
        self.err_age_edit.setText(', '.join(self.settings['err_age']))
        self.err_age_edit.setToolTip(
            'Two percentages relative to fitted age (e.g. 10%, -10%), '
            'press Enter to apply')
        self.err_age_edit.returnPressed.connect(self.on_err_age_changed)
        err_box.addWidget(self.show_err_age_w)
        err_box.addWidget(self.err_age_edit)
        vbox.addLayout(err_box)

        self.setLayout(vbox)
        for cb in [self.fixed_w,]:
            cb.toggled.connect(self.emit_change_signal)
        for edit in [self.y_intercept,]:
            edit.returnPressed.connect(self.emit_change_signal)

    def setSetting(self, key, value):
        self.settings[key] = value
    def setBool(self, key, text):
        try:
            value = bool(text)
            self.settings[key] = value
            self.emit_change_signal()
        except Exception: return
    def setFloat(self, key, text):
        try:
            value = float(text)
            self.settings[key] = value
            self.emit_change_signal()
        except Exception:
            return
    def on_return(self, key, edit):
        t = edit.text
        self.settings[key] = float(t)

    def on_ref_ages_changed(self):
        """解析逗号分隔的年龄列表，按回车生效；非法输入恢复原值。

        每项可以是绝对年龄（数字，单位 Ma，须 > 0），
        或相对百分数（如 +5%、-10%，须 > -100%，基准为当前拟合年龄）。
        统一规范化后存为字符串 token 列表，绘图时再换算。
        """
        text = self.ref_ages_edit.text
        tokens = []
        ok = True
        for part in re.split(r'[,，;；\s]+', text):
            part = part.strip()
            if not part:
                continue
            if part.endswith('%'):
                try:
                    pct = float(part[:-1].strip())
                except ValueError:
                    ok = False
                    break
                if pct <= -100:  # 年龄会 <= 0，无意义
                    ok = False
                    break
                tokens.append(f'{pct:+g}%')
            else:
                try:
                    v = float(part)
                except ValueError:
                    ok = False
                    break
                if v <= 0:
                    ok = False
                    break
                tokens.append(f'{v:g}')
        if ok and tokens:
            self.settings['ref_ages'] = tokens
            self.emit_change_signal()
        else:
            self.ref_ages_edit.setText(
                ', '.join(self.settings['ref_ages']))

    def on_err_age_changed(self):
        """解析年龄误差百分数（如 10%, -10%），须全部为百分数且 > -100%。
        取所有有效 token 的 min/max 作为两条边界线。按回车生效；非法恢复原值。
        """
        text = self.err_age_edit.text
        pcts = []
        ok = True
        for part in re.split(r'[,，;；\s]+', text):
            part = part.strip()
            if not part:
                continue
            if not part.endswith('%'):
                ok = False
                break
            try:
                pct = float(part[:-1].strip())
            except ValueError:
                ok = False
                break
            if pct <= -100:
                ok = False
                break
            pcts.append(pct)
        if ok and len(pcts) >= 2:
            pcts = [min(pcts), max(pcts)]
            self.settings['err_age'] = [f'{pcts[0]:+g}%', f'{pcts[1]:+g}%']
            self.emit_change_signal()
        else:
            self.err_age_edit.setText(', '.join(self.settings['err_age']))

    def getSettings(self):
        D = {}
        D['fixed'] = self.settings['fixed']
        D['y_intercept'] = self.settings['y_intercept']
        D['xmin'] = (self.settings['xmin'] if self.settings['xmin_enabled']
                     else None)
        D['xmax'] = (
            self.settings['xmax']
            if (self.settings['xmax_enabled'] and
                self.settings['xmax'] < float('inf'))
            else None)
        D['ymin'] = (self.settings['ymin'] if self.settings['ymin_enabled']
                     else None)
        D['ymax'] = (
            self.settings['ymax']
            if (self.settings['ymax_enabled'] and
                self.settings['ymax'] < float('inf'))
            else None)
        for key in ['errorbar', 'ellipse', 'ellipse2']: D[key] = self.settings[key]
        D['show_ref_lines'] = self.settings['show_ref_lines']
        D['ref_ages'] = list(self.settings['ref_ages'])
        D['show_err_age'] = self.settings['show_err_age']
        D['err_age'] = list(self.settings['err_age'])
        return D
    def emit_change_signal(self, *args, **kwargs):
        self.setting_changed.emit()

    def set_y_interp(self, y_interp):
        if y_interp == y_interp:
            self.fixed_w.setChecked(True)
            self.fixed_w.toggled.emit(True)
            self.y_intercept.setText(f'{y_interp}')
            self.y_intercept.returnPressed.emit()
        else:
            self.fixed_w.setChecked(False)
            self.fixed_w.toggled.emit(False)

l238 = 1.55125e-10
l235 = 9.8485e-10
l232 = 0.49475e-10
k = 137.818
lut = np.logspace(0, 23, 1000, base=2.7)
lu76 = (1/k)*(np.exp(l235*lut) - 1)/(np.exp(l238*lut) - 1)

def agePb76(y_intercept):
    try:
        return np.interp(y_intercept, lu76, lut) * 1e-6 + 1
    except Exception as e:
        print('agePb76', type(e), e)
        return float('nan')

class ColorChannelWidget(QWidget):
    color_channel_changed = Signal()
    
    def __init__(self):
        QWidget.__init__(self)
        self.settings = {'color_channel': 'fixed_orange'}
        
        thebox = QtGui.QHBoxLayout()
        thebox.setContentsMargins(0, 0, 0, 0)
        self.setLayout(thebox)
        
        self.color_cb = QtGui.QComboBox(self)
        self.color_cb.currentTextChanged.connect(self.on_color_channel_changed)
        thebox.addWidget(QLabel("Color Source:"))
        thebox.addWidget(self.color_cb)
        
        self.range_label = QLabel("Color Range: -")
        thebox.addWidget(self.range_label)
        
        self.refresh_btn = QtGui.QPushButton("Refresh Color Range")
        self.refresh_btn.clicked.connect(self.emit_refresh)
        thebox.addWidget(self.refresh_btn)
        
        self.update_channels()

    def update_channels(self):
        self.color_cb.blockSignals(True)
        current = self.color_cb.currentText
        self.color_cb.clear()
        
        channels = data.timeSeriesNames()
        
        self.color_cb.addItem("fixed_orange", "固定橙色")
        
        self.color_cb.insertSeparator(1)
        
        for ch in channels:
            self.color_cb.addItem(ch)
        
        if current and self.color_cb.findText(current) >= 0:
            self.color_cb.setCurrentText(current)
        else:
            self.color_cb.setCurrentIndex(0)
        
        self.color_cb.blockSignals(False)

    def update_for_method(self, method):
        self.color_cb.blockSignals(True)
        current = self.color_cb.currentText
        
        priority_elements = PPM_FOR_METHODS.get(method, [])
        
        priority_channels = []
        for elem in priority_elements:
            ppm_n_channel = f"{elem}_ppm_n"
            if ppm_n_channel in data.timeSeriesNames():
                priority_channels.append(ppm_n_channel)
            
            ppm_channel = f"{elem}_ppm"
            if ppm_channel in data.timeSeriesNames() and ppm_channel not in priority_channels:
                priority_channels.append(ppm_channel)
        
        self.color_cb.clear()
        
        self.color_cb.addItem("fixed_orange", "固定橙色")
        
        if priority_channels:
            for ch in priority_channels:
                self.color_cb.addItem(ch)
            self.color_cb.insertSeparator(len(priority_channels) + 1)
            
            all_channels = data.timeSeriesNames()
            other_channels = [ch for ch in all_channels 
                             if ch not in priority_channels and 
                             not ch.startswith(tuple(priority_elements))]
            
            for ch in other_channels:
                self.color_cb.addItem(ch)
            
            if current not in priority_channels:
                self.color_cb.setCurrentText(priority_channels[0])
                self.settings['color_channel'] = priority_channels[0]
        else:
            self.color_cb.addItem("fixed_orange", "固定橙色")
            self.color_cb.insertSeparator(1)
            for ch in data.timeSeriesNames():
                self.color_cb.addItem(ch)
            
            self.color_cb.setCurrentIndex(0)
            self.settings['color_channel'] = 'fixed_orange'
        
        if (current and current != 'fixed_orange' and 
            self.color_cb.findText(current) >= 0):
            self.color_cb.setCurrentText(current)
        
        self.color_cb.blockSignals(False)
        self.color_channel_changed.emit()

    def on_color_channel_changed(self, text):
        if text == 'fixed_orange':
            self.settings['color_channel'] = None
            self.refresh_btn.setVisible(False)
        else:
            self.settings['color_channel'] = text
            self.refresh_btn.setVisible(True)
        
        self.color_channel_changed.emit()

    def getSettings(self):
        return self.settings.copy()

    def set_range_label(self, text):
        self.range_label.setText(f"Color Range: {text}")

    def emit_refresh(self):
        self.color_channel_changed.emit()
