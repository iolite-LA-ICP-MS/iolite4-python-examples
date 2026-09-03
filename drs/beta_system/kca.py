from .base_drs import BaseDRS

from iolite import QtGui
import numpy as np
from functools import partial

class KCaDRS(BaseDRS):
    """K-Ca定年体系实现"""
    
    def __init__(self, data, drs):
        super().__init__(data, drs)
        self.ALIAS = {f'Ca{mass}': [f'Ca{mass+19}'] for mass in [40,42,44]}
        self.ALIAS.update({
            'K39o': ['K39'],
            })
        self.CALCULATED_ISOTOPES = {
            'K40': {'source': 'K39', 'factor': 0.0001289},
            }
        self.RATIOS = [
            ['K40', 'Ca42'], ['Ca40', 'Ca42'], ['K40', 'Ca40'],
            ['Ca42', 'Ca40'], ['K40', 'Ca44'], ['Ca40', 'Ca44'],
            ['Ca44', 'Ca40'],
        ]

        self.ratio_names = ['K39', 'K40', 'Ca40', 'Ca42', 'Ca44',]
        self.interfere_names = ['K39o', 'K39r', 'Ca44o', 'Ca44r',]
        self.Names_expected = self.ratio_names + self.interfere_names
        self.system_name = 'KCa'
        self.default_rm = 'ARM-1'

        self.PPMS = [['K39', 'K'], ['Ca44', 'Ca']]
        self.RHOS = [
            ['final_K40_Ca40', 'final_Ca44_Ca40'],
            ['final_K40_Ca44', 'final_Ca40_Ca44'],
            ]

        self.MATRIX_CORRS = ['K40_Ca40', 'K40_Ca44']

    def run(self):
        self.common_initialization()
        indexChannel = self.data.timeSeriesList()[0]
        self.do_ratios()
        self.do_matrix_corr()
        self.do_ppm()
        self.do_rho()
        print('KCa calculation finished')

    def settingsWidget(self):
        widget = self.make_widget()
        return widget

    def get_cps(self, elem):
        """获取计数率"""
        if elem == 'K40':
            f = 0.0001289
            return self.get_cps('K39') * f
        else:
            return super().get_cps(elem)

