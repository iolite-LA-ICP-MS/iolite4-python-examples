from .base_drs import BaseDRS

from iolite import QtGui
import numpy as np
from functools import partial

class LuHfDRS(BaseDRS):
    """Lu-Hf定年体系实现"""
    
    def __init__(self, data, drs):
        super().__init__(data, drs)
        self.ALIAS = {
            'Lu177': ['Lu207'],
            'Yb177': ['Yb204', 'Yb205'],
            'Lu175o': ['Lu175'],
            'Lu175r': ['Lu177', 'Lu207'],
            'Yb172o': ['Yb172'],
            'Yb172r': ['Yb177', 'Yb204', 'Yb205'],
        }
        self.CALCULATED_ISOTOPES = {
            'Lu176': {'source': 'Lu175', 'factor': 0.0259/0.9741},
            }
        self.RATIOS = [
            ['Hf176', 'Lu175'],
            ['Lu175', 'Hf176'],
            ['Lu176', 'Hf176'],
            ['Hf177', 'Hf176'],
            ['Lu176', 'Hf177'],
            ['Hf176', 'Hf177'],
            ['Lu175', 'Hf177'],
            ['Yb172o', 'Hf177'],
        ]
        self.RATIOS_MORE = [
            ['Lu175r', 'Lu175o'],
            ['Yb172r', 'Yb172o']
            ]

        self.ratio_names = ['Lu175', 'Lu176', 'Hf176', 'Hf177', 'Hf178',]
        self.interfere_names = ['Lu175o', 'Lu175r', 'Yb172o', 'Yb172r',]
        self.Names_expected = self.ratio_names + self.interfere_names
        self.system_name = 'LuHf'
        self.default_rm = 'G_NIST610'

        self.PPMS = [['Lu175', 'Lu'], ['Hf178', 'Hf']]
        self.RHOS = [
            ['final_Lu176_Hf176', 'final_Hf177_Hf176'],
            ['final_Lu176_Hf177', 'final_Hf176_Hf177'],
            ]
        self.MATRIX_CORRS = []

    def run(self):
        """执行Lu-Hf计算"""
        self.common_initialization()

        self.do_ratios()

        # Lu-Hf特有计算
        self.drs.message(f'{self.system_name} Interference substraction...')
        try:
            try:
                Lu175r = self.get_cps('Lu175r')
                Yb172r = self.get_cps('Yb172r')
                has_Lu175r_Yb172r = True
            except Exception as e:
                Lu175r = self.get_cps('Hf176') * np.nan
                Yb172r = self.get_cps('Hf176') * np.nan
                has_Lu175r_Yb172r = False
                print('failed getting Lu175r, Yb172r', type(e), e)
            # ILu = (Lu177*2.60/97.40)/Hf176
            ILu = Lu175r * 2.6 / 97.4 / self.get_cps('Hf176')
            # IYb = (Yb177*12.73/21.82)/Hf176
            IYb = Yb172r * 12.73 / 21.82 / self.get_cps('Hf176')

            # RLu%
            try:
                RLu_percent = Lu175r / (Lu175r - self.get_cps('Lu175o')) * 100
            except Exception as e:
                print('RLu% failed', e)
                RLu_percent = Lu175r * np.nan
            # RYb%
            try:
                RYb_percent = Yb172r / (Yb172r - self.get_cps('Yb172o')) * 100
            except Exception as e:
                print('RYb% failed', e)
                RYb_percent = Yb172r * np.nan

            # final_Lu176_Hf176_interCor = final_Lu176_Hf176 / (1- Ilu - IYb)
            r66 = self.data.timeSeries('final_Lu176_Hf176').data()
            r66_c = r66 / (1 - ILu - IYb)
            # final_Lu176_Hf177_interCor = final_Lu176_Hf177 / (1-0)
            r67 = self.data.timeSeries('final_Lu176_Hf177').data()
            r67_c = r67 / (1 - 0)
            # final_Hf177_Hf176_interCor = final_Hf177_Hf176 / (1- Ilu - IYb)
            r76 = self.data.timeSeries('final_Hf177_Hf176').data()
            r76_c = r76 / (1-ILu-IYb)
            # final_Hf176_Hf177_interCor = final_Hf176_Hf177 * (1- Ilu - IYb)
            rHf67 = self.data.timeSeries('final_Hf176_Hf177').data()
            rHf67_c = rHf67 * (1-ILu-IYb)
            # cHf = (Hf177 * 5.23 /18.55) / Hf176 *100
            cHf = self.get_cps('Hf177') *5.23/18.55 / self.get_cps('Hf176')*100
            # final_Lu176_Hf176_interComCor = final_Lu176_Hf176_interCor  / (1- cHf/100)
            r66_cc = r66_c / (1-cHf/100)
            # t = Ln (1 + (1/ final_Lu176_Hf176_interComCor) / (1.865 * 10^-11) * (1* 10^-6)
            t = np.log(1+(1/r66_cc))/1.865*1e5

            indexChannel = self.data.timeSeriesList()[0]
            for name, signal, kind in [
                ['ILu', ILu, self.data.Intermediate],
                ['IYb', IYb, self.data.Intermediate],
                ['RLu%', RLu_percent, self.data.Intermediate],
                ['RYb%', RYb_percent, self.data.Intermediate],
                ['final_Lu176_Hf176_interCor', r66_c, self.data.Output],
                ['final_Hf177_Hf176_interCor', r76_c, self.data.Output],
                ['final_Lu176_Hf177_interCor', r67_c, self.data.Output],
                ['final_Hf176_Hf177_interCor', rHf67_c, self.data.Output],
                ['cHf', cHf, self.data.Intermediate],
                ['final_Lu176_Hf176_interComCor', r66_cc, self.data.Output],
                ['time', t, self.data.Output],
                ]:
                self.data.createTimeSeries(
                    name, kind, indexChannel.time(), signal)

            # matrix factor correction
            F = self.drs.settings()[f'{self.system_name}MatrixFactor']
            # final_Lu176_Hf176_interCor_n = F * final_Lu176_Hf176_interCor
            r66_cn = F * r66_c
            # final_Lu176_Hf177_interCor_n = F * final_Lu176_Hf177_interCor
            r67_cn = F * r67_c
            # final_Lu176_Hf176_n = F * final_Lu176_Hf176
            r66_n = F * r66
            # final_Lu176_Hf177_n = F * final_Lu176_Hf177
            r67_n = F * r67
            # final_Lu176_Hf176_interComCor_n = f * final_Lu176_Hf176_interComCor
            r66_ccn = F * r66_cc
            # time_n = Ln (1 + (1/ final_Lu176_Hf176_interComCor_n) / (1.865 * 10^-11) * (1* 10^-6) 
            t_n = np.log(1+(1/r66_ccn))/1.865*1e5
            for name, signal, kind in [
                ['final_Lu176_Hf176_interCor_n', r66_cn, self.data.Output],
                ['final_Lu176_Hf177_interCor_n', r67_cn, self.data.Output],
                ['final_Lu176_Hf176_n', r66_n, self.data.Output],
                ['final_Lu176_Hf177_n', r67_n, self.data.Output],
                ['final_Lu176_Hf176_interComCor_n', r66_ccn, self.data.Output],
                ['time_n', t_n, self.data.Output],
                ]:
                self.data.createTimeSeries(
                    name, kind, indexChannel.time(), signal)
        except Exception as e:
            self.print_error('doing matrix corr', e)

        self.do_ppm()
        self.do_rho()
        print('luhf calculation finished')

    def settingsWidget(self):
        """Lu-Hf特有设置界面"""
        widget = self.make_widget()
        formLayout = self.params_form
        intercor_checkbox = QtGui.QCheckBox(widget)
        should_intercor = self.drs.settings().get('LuHf_intercor', True)
        if 'LuHf_intercor' not in self.drs.settings():
            self.drs.settings()['LuHf_intercor'] = should_intercor
        intercor_checkbox.setChecked(should_intercor)
        intercor_checkbox.toggled.connect(
            partial(self.set_setting, 'LuHf_intercor', bool))
        self.add_form_row(
            formLayout, 'LuHf_intercor', intercor_checkbox, self.font_size)

        self.area_factor_edit = QtGui.QLineEdit(widget)
        area_factor = self.drs.settings().get('LuHf_area_factor', 1)
        if 'LuHf_area_factor' not in self.drs.settings():
            self.drs.settings()['LuHf_area_factor'] = area_factor
        self.area_factor_edit.setText(str(area_factor))
        self.area_factor_edit.returnPressed.connect(self.set_area_factor)
        self.area_factor_edit.textChanged.connect(partial(
            self.set_setting, 'LuHf_area_factor', float))
        self.add_form_row(
            formLayout, 'area_factor(sample/rm)', self.area_factor_edit,
            self.font_size)

        return widget


    def set_area_factor(self, text=None):
        if text is None:
            text = self.area_factor_edit.text
        self.set_setting('LuHf_area_factor', float, text)
