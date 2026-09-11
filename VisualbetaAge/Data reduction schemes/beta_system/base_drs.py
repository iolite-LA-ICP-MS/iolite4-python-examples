from iolite import QtGui
from iolite.QtGui import QSizePolicy, QWidget
from iolite.types import Result
import numpy as np

from functools import partial
import itertools

class BaseDRS:
    """所有定年体系的基类"""
    
    def __init__(self, data, drs):
        self.data = data
        self.drs = drs
        self.logger = None
        self.ALIAS = {}  # 各体系可以覆盖这个属性
        self.name_map = {}
        self.CALCULATED_ISOTOPES = {}
        self.NO_MATCH_STR = 'No match' + '-'*8
        self.system_name = 'base'
        self.RATIOS = []
        self.RATIOS_MORE = []
        self.Names_expected = []
        self.font_size = 12

    def run(self):
        """执行DRS计算"""
        raise NotImplementedError

    def settingsWidget(self):
        """返回该体系的设置界面"""
        raise NotImplementedError

    def set_logger(self, logger):
        self.logger = logger

    def add_form_row(self, layout, label_text, field_widget, font):
        if isinstance(font, int):
            font_obj = QtGui.QFont()
            font_obj.setPointSize(font)
            font = font_obj

        label = QtGui.QLabel(label_text)
        label.setFont(font)
        layout.addRow(label, field_widget)
        if hasattr(field_widget, 'setFont'):
            field_widget.setFont(font)


    def make_widget(self):
        widget = QtGui.QWidget()
        self.dynamic_widget = QtGui.QWidget(widget)
        vbox = QtGui.QVBoxLayout()
        hbox = QtGui.QHBoxLayout()
        vbox.addLayout(hbox)
        vbox.addWidget(self.dynamic_widget)
        widget.setLayout(vbox)
        widget.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)

        #### font size
        font = QtGui.QFont()
        font.setPointSize(self.font_size)
        widget.setFont(font)


        # make some formLayout and add them to places
        inner_vbox = QtGui.QVBoxLayout()
        self.dynamic_widget.setLayout(inner_vbox)
        ref_form = QtGui.QFormLayout()
        ratio_form = QtGui.QFormLayout()
        interfere_form = QtGui.QFormLayout()
        params_form = QtGui.QFormLayout()
        self.params_form = params_form

        ratio_group_box = QtGui.QGroupBox("Isochron ratios")
        interference_group_box = QtGui.QGroupBox("Interferences and corrections")
        normalization_group_box = QtGui.QGroupBox("Normalization")
        ratio_group_box.setLayout(ratio_form)
        interference_group_box.setLayout(interfere_form)
        normalization_group_box.setLayout(params_form)
        for gbox in [ratio_group_box, interference_group_box,
                     normalization_group_box]:
            gbox.setFont(font)

        inner_vbox.addLayout(ref_form)
        inner_vbox.addWidget(ratio_group_box)
        inner_vbox.addWidget(interference_group_box)
        inner_vbox.addWidget(normalization_group_box)

        input_channels = list(self.data.timeSeriesNames(self.data.Input))
        input_channels.insert(0, self.NO_MATCH_STR)
        try: input_channels.remove('TotalBeam')
        except ValueError: pass

        system_enabled = f'{self.system_name}_enabled'
        if system_enabled not in self.drs.settings():
            self.set_setting(system_enabled, bool, False)
        checkbox = QtGui.QCheckBox(widget)
        checkbox.setChecked(self.drs.settings()[system_enabled])
        checkbox.toggled.connect(partial(
            self.set_setting, system_enabled, bool))
        checkbox.toggled.connect(self.dynamic_widget.setVisible)
        self.dynamic_widget.setVisible(self.drs.settings()[system_enabled])
        if hasattr(self, 'system_name_to_display'):
            label_text = self.system_name_to_display
        else:
            label_text = self.system_name
        lb = QtGui.QLabel(label_text)
        lb.setFont(font)
        hbox.addWidget(lb)
        hbox.addWidget(checkbox)
        hbox.addStretch()


        rmNames = self.data.selectionGroupNames(self.data.ReferenceMaterial)
        rmcbName = f'{self.system_name}_RefName'
        rmcbText = 'RefMaterial'
        rmcbValue = self.drs.settings().get(rmcbName, '')
        rmComboBox = QtGui.QComboBox(widget)
        rmComboBox.addItems(rmNames)
        if rmcbName in self.drs.settings():
            rmcbValue = self.drs.settings().get(rmcbName)
            if rmcbValue not in rmNames:
                if rmNames:
                    rmcbValue = rmNames[0]
                else:
                    raise Exception('No rm available')
        else:
            for rm in rmNames:
                if self.default_rm.lower() in rm.lower():
                    rmcbValue = rm
                    break
            else:
                if rmNames:
                    rmcbValue = rmNames[0]
                else:
                    raise Exception('No rm available')
        rmComboBox.setCurrentText(rmcbValue)
        self.set_setting(rmcbName, str, rmcbValue)

        rmComboBox.currentTextChanged.connect(partial(
            self.set_setting, rmcbName, str))
        self.add_form_row(ref_form, rmcbText, rmComboBox, font)

        self.all_channel_cb = []
        for name in self.Names_expected:
            if name in self.ratio_names:
                the_form = ratio_form
            else:
                the_form = interfere_form

            cb = QtGui.QComboBox(widget)
            self.all_channel_cb.append(cb)
            cb.addItems(input_channels)

            if name in self.CALCULATED_ISOTOPES:
                source, factor = self.get_calc_iso_info(name)

            if name in self.drs.settings():
                cb.setCurrentText(self.drs.settings()[name])
            else:
                if name in self.ALIAS:
                    candi = self.ALIAS[name] + [name]
                elif name in self.CALCULATED_ISOTOPES:
                    candi = [source]
                else:
                    candi = [name]
                for alias in candi:
                    if alias in input_channels:
                        cb.setCurrentText(alias)
                        break
                else:
                    cb.setCurrentText(input_channels[0])
            self.set_setting(name, str, cb.currentText)
            cb.currentTextChanged.connect(partial(
                self.set_setting, name, str))
            self.add_form_row(the_form, f'{name} name', cb, font)

            if name in self.CALCULATED_ISOTOPES:
                edit = QtGui.QLineEdit()
                if factor > 1e-2:
                    edit.setText(f'{factor: .4f}')
                else:
                    edit.setText(f'{factor}')
                edit.textChanged.connect(partial(
                    self.set_setting, f'{self.system_name}_{name}_factor', float))
                self.add_form_row(the_form, f'{name} factor', edit, font)

        matrixFactorName = f'{self.system_name}MatrixFactor'
        matrixFactorEdit = QtGui.QLineEdit(widget)
        matrixFactorEdit.textChanged.connect(
            partial(self.set_setting, matrixFactorName, float))
        matrix_factor_value = self.drs.settings().get(matrixFactorName, 1)
        matrixFactorEdit.setText(matrix_factor_value)
        self.add_form_row(params_form, matrixFactorName, matrixFactorEdit, font)

        return widget
        
    def common_initialization(self):
        """公共初始化流程"""
        # 这里放原来runDRS()开头的公共代码
        self.make_name_map()
        
##    def error_corrcoef(self, channel1, channel2, sel):
##        result = Result()
##        try:
##            chn1 = self.data.timeSeries(channel1)
##            chn2 = self.data.timeSeries(channel2)
##        except RuntimeError as e:
##            self.print_error('error_corrcoef', e)
##            return result
##        array_1 = chn1.dataForSelection(sel)
##        array_2 = chn2.dataForSelection(sel)
##        result.setValue(np.corrcoef(array_1, array_2)[0,1])
##        return result

    def error_corrcoef(self, channel1, channel2, sel):
        """计算相关系数"""
        result = Result()
        try:
            chn1 = self.data.timeSeries(channel1)
            chn2 = self.data.timeSeries(channel2)
            array_1 = chn1.dataForSelection(sel)
            array_2 = chn2.dataForSelection(sel)
            result.setValue(np.corrcoef(array_1, array_2)[0,1])
        except RuntimeError as e:
            self.print_error('error_corrcoef', e)
        return result

    def get_calc_iso_info(self, elem):
        calc_info = self.CALCULATED_ISOTOPES[elem]
##        source = self.drs.settings().get(f'{self.system_name}_{elem}_source', calc_info['source'])
        source = self.drs.settings().get(elem, calc_info['source'])
        factor_setting = self.drs.settings().get(f'{self.system_name}_{elem}_factor', calc_info['factor'])
        return source, factor_setting

    def get_cps(self, elem):
        CPSname = elem+'_CPS'

        if elem in self.CALCULATED_ISOTOPES:
            source, factor_setting = self.get_calc_iso_info(elem)
            source_cps = self.get_cps(source)
            return source_cps * float(factor_setting)

        elif elem in self.name_map:
            chn = self.name_map[elem]+'_CPS'
            if chn in self.data.timeSeriesNames():
                return self.data.timeSeries(chn).data().copy()
        elif CPSname in self.data.timeSeriesNames():
            return self.data.timeSeries(CPSname).data()
        elif elem in self.data.timeSeriesNames():
            return self.data.timeSeries(elem).data()
        else:
            raise Exception(f'{elem} do not exist')

    def get_cps_spline(self, elem, group_name):
        CPSname = elem+'_CPS'

        if elem in self.CALCULATED_ISOTOPES:
            source, factor_setting = self.get_calc_iso_info(elem)
            source_cps_spl = self.get_cps_spline(source, group_name)
            return source_cps_spl * float(factor_setting)
        
        elif elem in self.data.timeSeriesNames():
            return self.data.spline(group_name, elem).data().copy()
        elif elem in self.name_map:
            chn = self.name_map[elem]+'_CPS'
            if chn in self.data.timeSeriesNames():
                return self.data.spline(group_name, chn).data().copy()
        elif CPSname in self.data.timeSeriesNames():
            return self.data.spline(group_name, CPSname).data().copy()
        raise Exception(f'{elem} do not exist')

    def get_ref_value(self, rm_name, numerator, denominator):
        try:
            ref_data = self.data.referenceMaterialData(rm_name)
            all_numerators = [numerator, numerator.strip('or')]
            if numerator in self.name_map:
                all_numerators.append(self.name_map[numerator])
            all_denominators = [denominator, denominator.strip('or')]
            if denominator in self.name_map:
                all_denominators.append(self.name_map[denominator])
            for n, d in zip(all_numerators, all_denominators):
                if f'{n}/{d}' in ref_data:
                    return ref_data[f'{n}/{d}'].value()
                elif f'{d}/{n}' in ref_data:
                    return 1 / ref_data[f'{d}/{n}'].value()
##            if f'{numerator}/{denominator}' in ref_data:
##                return ref_data[f'{numerator}/{denominator}'].value()
##            elif f'{denominator}/{numerator}' in ref_data:
##                return 1 / ref_data[f'{denominator}/{numerator}'].value()
            else:
                raise Exception(f'{numerator}/{denominator} not found')
        except Exception as e:
            self.print_error(f'getting ref {numerator}/{denominator}', e)
            raise

    def set_setting(self, key, value_type, value):
        self.drs.setSetting(key, value_type(value))

    def make_name_map(self):
        Names = [self.drs.settings()[k] for k in self.Names_expected]
        self.name_map.update(zip(self.Names_expected, Names))
        skips = []
        for k,v in list(self.name_map.items()):
            if v == self.NO_MATCH_STR:
                self.name_map.pop(k)
                if not k in ['Lu176', 'Rb87']: skips.append(k)

    def should_run(self):
        return self.drs.settings().get(f'{self.system_name}_enabled', False)

    def do_ratios(self):
        self.drs.message(f'Doing {self.system_name} ratios...')
        indexChannel = self.data.timeSeriesList()[0]

        for n1,n2 in self.RATIOS:
            rmName = self.drs.settings()[f'{self.system_name}_RefName']
            name = f'{n1}_{n2}'
            try:
                ratio = self.get_cps(n1) / self.get_cps(n2)
            except Exception as e:
                self.print_error(f'doing ratio n1 {n1} n2 {n2}', e)
                ratio = indexChannel.data() * np.nan

##            try:
##                x1 = self.get_cps(n1)
##            except Exception as e:
##                self.print_error(f'doing ratio n1 {n1} n2 {n2}', e)
##                x1 = indexChannel.data() * 0
##            try:
##                x2 = self.get_cps(n2)
##            except Exception as e:
##                self.print_error(f'doing ratio n1 {n1} n2 {n2}', e)
##                x2 = indexChannel.data() * 0
##            ratio = x1 / x2

            self.data.createTimeSeries(
                name, self.data.Intermediate,
                indexChannel.time(), ratio.copy())
##            spline = self.data.spline(rmName, name).data()
            spline = self.get_cps_spline(name, rmName)
            try:
                ref = self.get_ref_value(rmName, n1, n2)
##                ref_data = self.data.referenceMaterialData(rmName)
##                ref = ref_data[f'{n1}/{n2}'].value()
            except Exception as e:
                self.print_error(f'getting ref {n1}, {n2}', e)
                continue
            try:
                corrected = ratio / spline * ref
            except Exception as e:
                print(n1, n2, type(e), e)
                continue
            self.data.createTimeSeries(
                f'final_{name}', self.data.Output,
                indexChannel.time(), corrected.copy())

        for n1,n2 in self.RATIOS_MORE:
            name = f'{n1}_{n2}'
            try:
                ratio = self.get_cps(n1) / self.get_cps(n2)
            except Exception as e:
                self.print_error(f'ratios more {n1}, {n2}', e)
                ratio = indexChannel.data() * np.nan
            self.data.createTimeSeries(
                name, self.data.Intermediate,
                indexChannel.time(), ratio.copy())

    def baseline_sub(self):
        allInputChannels = self.data.timeSeriesList(self.data.Input)
        blGrp = None

        all_bs = self.data.selectionGroupList(self.data.Baseline)
        if len(all_bs) > 1:
            raise Exception('There are multiple baseline groups.')
        elif len(all_bs) < 1:
            raise Exception("No baselines.")
        else:
            blGrp = all_bs[0]
        if len(blGrp.selections()) < 1:
            raise Exception('No baseline selections.')
        indexChannel = self.data.timeSeriesList()[0]
        mask = np.ones_like(indexChannel.time())
        self.drs.message("Baseline subtracting")
        self.drs.baselineSubtract(
            blGrp, allInputChannels, mask, 25, 50)

    def update_channels(self):
        for cb in self.all_channel_cb:
            old_text = cb.currentText
            cb.clear()
            timeSeriesNames = list(self.data.timeSeriesNames(self.data.Input))
            cb.addItems([self.NO_MATCH_STR] + timeSeriesNames)
            cb.setCurrentText(old_text)

    def print_error(self, msg, e):
        text = f'{self.system_name} {msg}, {e.__class__.__name__}, {e}'
        print(text)
        if self.logger:
            self.logger.warning(text)

    def do_matrix_corr(self):
        try:
            # matrix factor correction
            F = self.drs.settings()[f'{self.system_name}MatrixFactor']
            for chn in self.data.timeSeriesList(self.data.Output):
                for name in self.MATRIX_CORRS:
##                    if chn.name.startswith(f'final_{name}'):
                    if name in chn.name:
                        correction = chn.data() * F
                        self.data.createTimeSeries(
##                            chn.name + '_corr', self.data.Output,
                            chn.name + '_n', self.data.Output,
                            chn.time(), correction)
        except Exception as e:
            self.print_error('doing matrix correction', e)

    def do_rho(self):
        for chn1, chn2 in self.RHOS:
            try:
                asso_name = f'{chn1} - {chn2} Rho'
                self.data.registerAssociatedResult(
                    asso_name,
                    partial(self.error_corrcoef, chn1, chn2))
                self.freeze_associated(asso_name)
                self.data.removeAssociatedResult(asso_name)
            except Exception as e:
                self.print_error('doing rho {chn1},{chn2}', e)

    def freeze_associated(self, associated_name):
        print('entering freeze_associated')
        indexChannel = self.data.timeSeriesList()[0]
        X = indexChannel.data() * np.nan
        groups = self.data.selectionGroupNames()
##        print('groups', len(groups), groups)
        all_sel = itertools.chain(
            *[self.data.selectionGroup(grp).selections() for grp in groups])
        for sel in all_sel:
            asso_value = self.data.associatedResult(sel, associated_name).value()
            idx = indexChannel.selectionIndices(sel)
            X[idx] = asso_value
        self.data.createTimeSeries(
            f'{associated_name} frozen', self.data.Output,
            indexChannel.time(), X)

    def do_ppm(self):
        rmName = self.drs.settings()[f'{self.system_name}_RefName']
        try:
            area_factor = self.drs.settings()[f'{self.system_name}_area_factor']
        except KeyError:
            area_factor = 1
        print('area_factor', area_factor)
        the_time = self.data.timeSeriesList()[0].time()
        for chn, elem in self.PPMS:
            try:
##                spline = self.data.spline(rmName, f'{chn}_CPS')
                spline = self.get_cps_spline(chn, rmName)
                signal = self.get_cps(chn) / spline
                ref = self.data.referenceMaterialData(rmName)[elem].value()
                ppm = signal * ref
                ppm_n = signal * ref / area_factor
                self.data.createTimeSeries(
                    elem+'_ppm', self.data.Output, the_time, ppm)
                self.data.createTimeSeries(
                    elem+'_ppm_n', self.data.Output, the_time, ppm_n)
            except Exception as e:
                    self.print_error('doing ppm', e)
