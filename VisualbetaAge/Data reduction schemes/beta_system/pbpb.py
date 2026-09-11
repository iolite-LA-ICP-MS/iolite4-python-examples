from .base_drs import BaseDRS

from iolite import QtGui
import numpy as np
from functools import partial

class PbPbDRS(BaseDRS):
    """Pb-Pb定年体系实现"""
    
    def __init__(self, data, drs):
        super().__init__(data, drs)
        self.ALIAS = {
        }
        self.RATIOS = [
            ['Pb207', 'Pb206'], ['Pb204', 'Pb206'],
        ]

        self.Names_expected = sorted(set(sum(
            self.RATIOS + self.RATIOS_MORE, start=[])))
        self.ratio_names = self.Names_expected[:]
        self.interfere_names = []

        self.system_name = 'PbPb'
        self.default_rm = 'G_NIST610'
        self.PPMS = []
        self.RHOS = []
        self.MATRIX_CORRS = []

    def run(self):
        self.common_initialization()
        indexChannel = self.data.timeSeriesList()[0]
        self.do_ratios()
        print('PbPb calculation finished')

    def settingsWidget(self):
        """Rb-Sr设置界面"""
        widget = self.make_widget()
        return widget
