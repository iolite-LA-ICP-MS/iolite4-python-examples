from .base_drs import BaseDRS

from iolite import QtGui
import numpy as np
from functools import partial

class MoDauDRS(BaseDRS):
    def __init__(self, data, drs):
        super().__init__(data, drs)
        self.RATIOS = [
            ['Mother', 'Nonradio'], ['Daughter', 'Nonradio'],
            ['Mother', 'Daughter'], ['Nonradio', 'Daughter'],
        ]
        self.CALCULATED_ISOTOPES = {
            'Mother': {'source': 'Mother', 'factor': 1},
            }

        self.ratio_names = ['Mother', 'Daughter', 'Nonradio']
        self.interfere_names = []
        self.Names_expected = self.ratio_names + self.interfere_names
        self.default_rm = 'G_NIST610'
        self.system_name = 'MoDau'
        self.system_name_to_display = 'Other'

        self.PPMS = []
        self.RHOS = []
        self.MATRIX_CORRS = ['Mother_Daughter', 'Mother_Nonradio']

    def run(self):
        self.common_initialization()
        self.do_ratios()
        self.do_matrix_corr()
        self.do_ppm()
        self.do_rho()

        print('MoDau calculation finished')

    def settingsWidget(self):
        widget = self.make_widget()
        return widget
