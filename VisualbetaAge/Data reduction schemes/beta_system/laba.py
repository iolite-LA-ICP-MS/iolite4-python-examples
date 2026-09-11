from .base_drs import BaseDRS

from iolite import QtGui
import numpy as np
from functools import partial

class LaBaDRS(BaseDRS):
    """La-Ba定年体系实现"""
    
    def __init__(self, data, drs):
        super().__init__(data, drs)
        self.ALIAS = {
            'La139o': ['La139'],
            'Ce140o': ['Ce140'],
        }
        self.CALCULATED_ISOTOPES = {
            'La138': {'source': 'La139', 'factor': 0.0008889},
            }
        self.RATIOS = [
            ['La138', 'Ba137'], ['Ba138', 'Ba137'],
            ['Ba137', 'Ba138'], ['La138', 'Ba138'],
        ]

        self.ratio_names = ['La139', 'La138', 'Ba138', 'Ba137',]
        self.interfere_names = ['La139o', 'La139r', 'Ce140o', 'Ce140r',]
        self.Names_expected = self.ratio_names + self.interfere_names
        self.system_name = 'LaBa'
        self.default_rm = 'G_NIST610'

        self.PPMS = [['La139', 'La'], ['Ba137', 'Ba']]
        self.RHOS = [
            ['final_La138_Ba138', 'final_Ba137_Ba138'],
            ['final_La138_Ba137', 'final_Ba138_Ba137'],
            ]
        self.MATRIX_CORRS = []

    def run(self):
        self.common_initialization()

        indexChannel = self.data.timeSeriesList()[0]

        self.do_ratios()

        try:
            Ba138 = self.get_cps('Ba138')
            try:
                La139r = self.get_cps('La139r')
            except Exception as e:
                self.print_error('La139r undefined, using 0', e)
                La139r = Ba138 * 0
            try:
                Ce140r = self.get_cps('Ce140r')
            except Exception as e:
                self.print_error('Ce140r undefined, using 0', e)
                Ce140r = Ba138 * 0
            Ba138_star = Ba138 - La139r * 0.0008881/0.9991119 - Ce140r*0.00251/0.88449
            self.data.createTimeSeries(
                'Ba138Star', self.data.Output, indexChannel.time(),
                Ba138_star)
        except Exception as e:
            self.print_error('Ba138 correction', e)

        rmName = self.drs.settings()[f'{self.system_name}_RefName']
        for n1, n2 in self.RATIOS:
            name = f'{n1}_{n2}_interCor'
            try:
                N1 = 'Ba138Star' if n1=='Ba138' else n1
                N2 = 'Ba138Star' if n2=='Ba138' else n2
                ratio = self.get_cps(N1) / self.get_cps(N2)
            except Exception as e:
                self.print_error(f'doing ratio N1 {N1} N2 {N2}', e)
                ratio = indexChannel.data() * np.nan
            self.data.createTimeSeries(
                name, self.data.Intermediate, indexChannel.time(), ratio.copy())
            spline = self.get_cps_spline(name, rmName)
            try:
                ref = self.get_ref_value(rmName, n1, n2)
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


        self.do_matrix_corr()
        self.do_ppm()
        self.do_rho()

        print('LaBa calculation finished')

    def settingsWidget(self):
        widget = self.make_widget()

        formLayout = self.params_form
        intercor_checkbox = QtGui.QCheckBox(widget)
        should_intercor = self.drs.settings().get('LaCe_intercor', True)
        if 'LaCe_intercor' not in self.drs.settings():
            self.set_setting('LaCe_intercor', bool, should_intercor)
        intercor_checkbox.setChecked(should_intercor)
        intercor_checkbox.toggled.connect(
            partial(self.set_setting, 'LaCe_intercor', bool))
        self.add_form_row(
            formLayout, 'LaCe_intercor', intercor_checkbox, self.font_size)

        self.area_factor_edit = QtGui.QLineEdit(widget)
        area_factor = self.drs.settings().get('LaBa_area_factor', 1)
        if 'LaBa_area_factor' not in self.drs.settings():
            self.set_setting('LaBa_area_factor', float, area_factor)
        self.area_factor_edit.setText(str(area_factor))
        self.area_factor_edit.returnPressed.connect(self.set_area_factor)
        self.area_factor_edit.textChanged.connect(partial(
            self.set_setting, 'LaBa_area_factor', float))
        self.add_form_row(
            formLayout, 'area_factor(sample/rm)', self.area_factor_edit,
            self.font_size)

        return widget


    def set_area_factor(self, text=None):
        if text is None:
            text = self.area_factor_edit.text
        self.set_setting('LaBa_area_factor', float, text)
