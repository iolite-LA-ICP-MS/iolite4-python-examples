#/ Name: Compare with iolite 3
#/ Authors: Joe Petrus and Bence Paul
#/ Description: Compares iolite 3 and 4 results
#/ References: None
#/ Version: 1.0
#/ Contact: support@iolite-software.com

from iolite.QtGui import QLabel
from igor import packed
import pandas as pd
import numpy as np

MEAN_INDEX = 2
UNC_INDEX = 3

def update():	
	qaqc.clearReport()
	f = data.importedFiles()[0].filePath()
	if not f.endswith('pxp'):
		qaqc.pushHtml('This module must be used after importing an iolite 3 experiment. Aborting.')
		qaqc.finished(qaqc.Error)
		return

	d = packed.load(f)[1]

	int_folder = d['root'][b'Packages'][b'iolite'][b'integration']

	matrix_names = [s for s in int_folder.keys() if s.decode('utf-8').startswith('m_') and s != b'm_Baseline_1']

	dfi3 = pd.DataFrame()

	for m in matrix_names:		
		print('Processing: %s'%(m))
		aim = int_folder[m].wave['wave']['wData']
		labels = int_folder[m].wave['wave']['labels'][1][2:]
		tim = int_folder[m.decode('utf-8').replace('m_', 't_').encode('utf-8')].wave['wave']['wData']
		sns = np.arange(1, aim.shape[0])

		for sn in sns:
			s = pd.Series(aim[sn, 1:, MEAN_INDEX], index=labels)
			s.name = m.decode('utf-8').replace('m_', '')+str(sn)
			print('   %i %s %s'%(sn, s.name, str(s.shape)))
			dfi3 = dfi3.append(s)
		
	#print('Final df shape %s'%(str(dfi3.shape)))
	#qaqc.pushHtml(dfi3.to_html(col_space=100))

	dfdiff = pd.DataFrame()

	for sg_name in data.selectionGroupNames():
		if 'Baseline' in sg_name or 'NIST610' in sg_name:
			continue

		sg = data.selectionGroup(sg_name)
		seln = 1
		for sel in sg.selections():
			sel_name = sg_name+str(seln)
			s = pd.Series(index=dfi3.columns)
			s.name = sel_name
			for ch_name in dfi3.columns:
				i4_ch_name = ch_name.decode('utf-8')
				if 'SQ' in i4_ch_name:
					i4_ch_name = i4_ch_name.split("_")[0] + i4_ch_name.split("_")[3][1:] + "_ppm"
				elif 'ppm' in i4_ch_name:
					i4_ch_name = i4_ch_name.split("_")[0] + i4_ch_name.split("_")[2][1:] + "_ppm"
				elif i4_ch_name == 'Final206_238':
					i4_ch_name = 'Final Pb206/U238'
				elif i4_ch_name == 'Final207_235':
					i4_ch_name = 'Final Pb207/U235'
				elif i4_ch_name == 'Final207_206':
					i4_ch_name = 'Final Pb207/Pb206'
				elif i4_ch_name == 'Raw_206_238':
					i4_ch_name = 'Pb206/U238'
				elif i4_ch_name == 'Raw_207_235':
					i4_ch_name = 'Pb207/U235'

				try:
					ch = data.timeSeries(i4_ch_name)
				except RuntimeError:
					continue

				s[ch_name] = data.result(sel, ch).value()

			dfdiff = dfdiff.append(100*(dfi3.loc[sel_name, :] - s)/dfi3.loc[sel_name, :])		
			seln += 1

	qaqc.pushHtml('<h3>Percent difference between iolite 3 and 4:</h3>')
	qaqc.pushHtml(dfdiff.to_html(col_space=100))

	qaqc.finished(qaqc.Success)


def settingsWidget():
	qaqc.setSettingsWidget(QLabel('No settings for this module yet.'))
