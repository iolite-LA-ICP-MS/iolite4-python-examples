from .base_drs import BaseDRS

from iolite import QtGui
import numpy as np
from functools import partial

class RbSrDRS(BaseDRS):
    """Rb-Sr定年体系实现"""
    
    def __init__(self, data, drs):
        super().__init__(data, drs)
        self.ALIAS = {
            'Rb88': ['Rb104'],
            'Sr86': ['Sr105'],
            'Sr87': ['Sr106'],
            'Sr86F': ['Sr.F86', 'Sr105'],
            'Sr86r': ['Sr.F86', 'Sr105'],
            'Rb85r': ['Rb85.F', 'Rb104'],
            'Rb85o': ['Rb85'],
            'Sr88o': ['Sr88'],
            'Sr88r': ['Sr.F88', 'Sr107'],
        }
        isotope_abundance = {
            'Rb85': 1/1.3860, 'Rb87': 0.3860/1.3860}
        f = isotope_abundance['Rb87'] / isotope_abundance['Rb85']
        self.CALCULATED_ISOTOPES = {
            'Rb87': {'source': 'Rb85', 'factor': f},
            }
        self.RATIOS = [
            ['Rb87', 'Sr86'], ['Sr87', 'Sr86'],
            ['Rb87', 'Sr87'], ['Sr86', 'Sr87'],
        ]

        self.ratio_names = ['Rb85', 'Rb87', 'Sr86', 'Sr87', 'Sr88',]
        self.interfere_names = ['Rb85o', 'Rb85r', 'Sr86r', 'Sr88o', 'Sr88r',]
        self.Names_expected = self.ratio_names + self.interfere_names
        self.system_name = 'RbSr'
        self.default_rm = 'G_NIST610'

        self.PPMS = [['Rb85', 'Rb'], ['Sr88', 'Sr']]
        self.RHOS = [
            ['final_Rb87_Sr86', 'final_Sr87_Sr86'],
            ['final_Rb87_Sr87', 'final_Sr86_Sr87'],
            ]

        self.MATRIX_CORRS = ['Rb87_Sr87', 'Rb87_Sr86']

        self.all_columns = '''
Rrb%
Rsr%
ComSr%
Rb87_Sr86
Sr87_Sr86
Rb87_Sr87
Sr86_Sr87
final_Rb87_Sr86
final_Sr87_Sr86
final_Rb87_Sr87
final_Sr86_Sr87
final_Rb87_Sr86_n
final_Rb87_Sr87_n
final_Rb87_Sr87_ComCorr
final_Rb87_Sr87_ComCorr_n
ModelAge
ModelAge_n
Rb_ppm
Rb_ppm_n
Sr_ppm
Sr_ppm_n
final_Rb87_Sr86 - final_Sr87_Sr86 Rho frozen
final_Rb87_Sr87 - final_Sr86_Sr87 Rho frozen
'''.strip().splitlines()

    def run(self):
        """执行Rb-Sr计算"""
        self.common_initialization()

        indexChannel = self.data.timeSeriesList()[0]

        # 86SrF/88SrF
        try:
            Sr86F88F = self.get_cps('Sr86') / self.get_cps('Sr88')
            self.data.createTimeSeries(
                '86Sr_88Sr_ratio', self.data.Intermediate,
                indexChannel.time(), Sr86F88F)
        except Exception as e:
            self.print_error('doing 86Sr/88Sr', e)

        # Rrb% = 85RbF/(85Rb + 85RbF) * 100
        try:
            RrbPerc = (self.get_cps('Rb85r') /
                       (self.get_cps('Rb85o') +
                        self.get_cps('Rb85r'))) * 100
            self.data.createTimeSeries(
                'Rrb%', self.data.Intermediate, indexChannel.time(), RrbPerc)
        except Exception as e:
            self.print_error('doing Rrb%', e)

        # Rsr% = 88SrF/(88Sr + 88SrF) * 100
        try:
            RsrPerc = (self.get_cps('Sr88r')/ (self.get_cps('Sr88o') + self.get_cps('Sr88r'))) * 100
            self.data.createTimeSeries(
                'Rsr%', self.data.Intermediate, indexChannel.time(), RsrPerc)
        except Exception as e:
            print('doing Rsr%', type(e), e)

        self.do_ratios()
        self.do_matrix_corr()
        self.do_ppm()
        self.do_rho()

        try:
            self.calc_model_age()
        except Exception as e:
            self.print_error('calc_model_age', e)

        self.reset_columns()
        print('RbSr calculation finished')

    def calc_model_age(self):
##        print(self.drs.settings().get('RbSr_modelage'))
##        print(self.drs.settings().get('RbSr_initSr87Sr86'))
        if not self.drs.settings().get('RbSr_modelage', False):
            return
        initSr87Sr86 = self.drs.settings().get('RbSr_initSr87Sr86')
##        print(type(initSr87Sr86), initSr87Sr86)
        RbSrMatrixFactor = self.drs.settings().get('RbSrMatrixFactor')
##        print(type(RbSrMatrixFactor), RbSrMatrixFactor)
        if initSr87Sr86 is None or RbSrMatrixFactor is None:
            return

        try:
            final_Sr86_Sr87 = self.get_cps('final_Sr86_Sr87')
        except Exception as e:
            self.print_error('calc_model_age', e)
            return
        try:
            final_Rb87_Sr87 = self.get_cps('final_Rb87_Sr87')
        except Exception as e:
            self.print_error('calc_model_age', e)
            return

##        ComSrPercent = (1 - final_Sr86_Sr87 / (1/initSr87Sr86)) * 100
##        final_Rb87_Sr87_ComCorr = final_Rb87_Sr87 / (1 - ComSrPercent)
##        ModelAge = np.log(final_Rb87_Sr87_ComCorr + 1) * 1.3972 * 1e5
##        final_Rb87_Sr87_ComCorr_n = final_Rb87_Sr87_ComCorr * RbSrMatrixFactor
##        ModelAge_n = np.log(final_Rb87_Sr87_ComCorr_n + 1) / 1.3972 * 1e5

        ComSrPercent = (final_Sr86_Sr87 / (1/initSr87Sr86)) * 100
        final_Rb87_Sr87_ComCorr = final_Rb87_Sr87 / (1 - ComSrPercent/100)
        ModelAge = np.log(1 + (1/final_Rb87_Sr87_ComCorr)) / 1.3972 * 1e5
        final_Rb87_Sr87_ComCorr_n = final_Rb87_Sr87_ComCorr * RbSrMatrixFactor
        ModelAge_n = np.log(1+ (1/final_Rb87_Sr87_ComCorr_n)) / 1.3972 * 1e5

        indexChannel = self.data.timeSeriesList()[0]
        for name, signal in [
            ['ComSr%', ComSrPercent],
            ['final_Rb87_Sr87_ComCorr', final_Rb87_Sr87_ComCorr],
            ['ModelAge', ModelAge],
            ['final_Rb87_Sr87_ComCorr_n', final_Rb87_Sr87_ComCorr_n],
            ['ModelAge_n', ModelAge_n],
            ]:
                self.data.createTimeSeries(
                    name, self.data.Output, indexChannel.time(), signal)


    def settingsWidget(self):
        """Rb-Sr设置界面"""
        widget = self.make_widget()
        formLayout = self.params_form

        model_age_checkbox = QtGui.QCheckBox(widget)
        model_age_checkbox.toggled.connect(
            partial(self.set_setting, 'RbSr_modelage', bool))
        model_age_checkbox.setChecked(True)
        self.add_form_row(
            formLayout, 'ModelAge', model_age_checkbox, self.font_size)

        self.init_sr87_sr86_edit = QtGui.QLineEdit(widget)
        self.init_sr87_sr86_edit.textChanged.connect(partial(
            self.set_setting, 'RbSr_initSr87Sr86', float))
        self.init_sr87_sr86_edit.setText('0.700')
        self.add_form_row(
            formLayout, 'initial Sr87/Sr86', self.init_sr87_sr86_edit,
            self.font_size)


        return widget

##    def get_cps(self, elem):
##        """获取计数率"""
##        isotope_abundance = {
##            'Rb85': 1/1.3860, 'Rb87': 0.3860/1.3860}
##        if elem == 'Rb87':
##            f = isotope_abundance['Rb87'] / isotope_abundance['Rb85']
##            return self.get_cps('Rb85') * f
##        else:
##            return super().get_cps(elem)

    def reset_columns(self):
        all_channels = self.data.timeSeriesNames()
        for col in self.all_columns:
            col = col.strip()
            if col in all_channels:
                try:
                    chn = self.data.timeSeries(col)
                    x = chn.data().copy()
                    t = chn.time().copy()
                    self.data.removeTimeSeries(col)
                    self.data.createTimeSeries(col, self.data.Output, t, x)
                except Exception as e:
                    self.print_error('reset_columns', e)















