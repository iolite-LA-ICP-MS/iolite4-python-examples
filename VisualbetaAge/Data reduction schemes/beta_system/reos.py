from .base_drs import BaseDRS

from iolite import QtGui
import numpy as np
from functools import partial

class ReOsDRS(BaseDRS):
    """Re-Os定年体系实现"""
    
    def __init__(self, data, drs):
        super().__init__(data, drs)
        self.ALIAS = {
            'Re185o': ['Re185'],
            'Re185r': ['Re205'],
            'Os192o': ['Os192'],
            'Os192r': ['Os212'],
        }
        self.CALCULATED_ISOTOPES = {
            'Re187': {'source': 'Re185', 'factor': 1.6738},
            }
        self.RATIOS = [
            ['Re187', 'Os187'], ['Os192', 'Os187'], ['Re187', 'Os192'],
            ['Os187', 'Os192'],
        ]

        self.ratio_names = ['Re185', 'Re187', 'Os187', 'Os192',]
        self.interfere_names = ['Re185o', 'Re185r', 'Os192o', 'Os192r',]
        self.Names_expected = self.ratio_names + self.interfere_names
        self.system_name = 'ReOs'
        self.default_rm = 'G_NIST610'

        self.PPMS = [['Re185', 'Re'], ['Os192', 'Os']]
        self.RHOS = []

        self.MATRIX_CORRS = ['Os187StarStar_Re187', 'AgeOs187Re187StarStar']

    def run(self):
        self.common_initialization()

        indexChannel = self.data.timeSeriesList()[0]

        # Rc ReOs
        try:
            RcRe = (self.get_cps('Re185r') /
                       (self.get_cps('Re185r') +
                        self.get_cps('Re185o') ))
            self.data.createTimeSeries(
                'RcRe', self.data.Intermediate, indexChannel.time(), RcRe)
        except Exception as e:
            self.print_error('doing RcRe', e)
        try:
            RcOs = (self.get_cps('Os192r') /
                       (self.get_cps('Os192r') +
                        self.get_cps('Os192o') ))
            self.data.createTimeSeries(
                'RcOs', self.data.Intermediate, indexChannel.time(), RcOs)
        except Exception as e:
            self.print_error('doing RcOs', e)

        self.do_ratios()

        try:
            os187 = self.get_cps('Os187')
            re205 = self.get_cps('Re185r')
            os192 = self.get_cps('Os192')
            re187 = self.get_cps('Re187')
            os187_star = os187 - re205 * 62.6 / 37.4
            os187_star_star = os187_star - os192 * 1.6/41.0
            Rstar = os187_star / re187
            Rstar_star = os187_star_star / re187
            age_star = np.log(Rstar + 1) / 0.00001666
            age_star_star = np.log(Rstar_star + 1) / 0.00001666
            self.data.createTimeSeries(
                'Os187Star', self.data.Output, indexChannel.time(),
                os187_star)
            self.data.createTimeSeries(
                'Os187StarStar', self.data.Output, indexChannel.time(),
                os187_star_star)
            self.data.createTimeSeries(
                'Os187Star_Re187', self.data.Output, indexChannel.time(),
                Rstar)
            self.data.createTimeSeries(
                'Os187StarStar_Re187', self.data.Output, indexChannel.time(),
                Rstar_star)
            self.data.createTimeSeries(
                'AgeOs187Re187Star', self.data.Output, indexChannel.time(),
                age_star)
            self.data.createTimeSeries(
                'AgeOs187Re187StarStar', self.data.Output, indexChannel.time(),
                age_star_star)
        except Exception as e:
            self.print_error('ReOs correction', e)

        self.do_matrix_corr()
        self.do_ppm()

        print('ReOs calculation finished')

    def settingsWidget(self):
        widget = self.make_widget()
        return widget


