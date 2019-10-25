from iolite.QtGui import QInputDialog
from iolite.QtCore import QDateTime
#from sklearn.cluster import MeanShift, DBSCAN, OPTICS, SpectralClustering
#from sklearn.cluster import MeanShift, DBSCAN, SpectralClustering
from sklearn import preprocessing
import pandas as pd
import numpy as np

from iolite import BoolResult

okObj = BoolResult()

input_group_name = QInputDialog.getItem(None, "Auto Selection Adjuster", "Group:", data.selectionGroupNames(), 0, False, okObj)

if not okObj:
	raise RuntimeError('User cancelled -- abort!')

if not input_group_name:
	raise RuntimeError('No group supplied -- abort!')

group = data.selectionGroup(input_group_name)
output_group = data.createSelectionGroup(input_group_name + '_auto_adj', group.type)

# Specify which channels to use for the clustering:
channel_names = ['Si29_ppm', 'Mg24_ppm', 'Al27_ppm', 'Fe57_ppm', 'Sr88_ppm', 'Zr90_ppm']
#channel_names = data.timeSeriesNames(data.Output)

for selection in group.selections():
	# Create a data frame of the selection's data for channels specified above:
	d = {}
	for channel_name in channel_names:
		d[channel_name] = data.timeSeries(channel_name).dataForSelection(selection)
	df = pd.DataFrame(d)

	df['Zr_norm'] = df.Zr90_ppm / df.Zr90_ppm.iloc[0]
	df['Zr_roll_std'] = df.Zr_norm.rolling(5).std()

	if '1859' in selection.name:
		df.to_csv('/Users/bence/ownCloud/testing/Rio Clustering/1859_data.csv')

	endIndex = -1

	for p in df['Zr_roll_std'] :
		if np.isnan(p) or p < 2:
			endIndex += 1
		else:
			break

	t = data.timeSeries('Zr90_ppm').timeForSelection(selection)

	start_ms = int(1000*(t[0]))
	end_ms =  int(1000*(t[endIndex]))
	new_s = data.createSelection(output_group, QDateTime.fromMSecsSinceEpoch(start_ms), QDateTime.fromMSecsSinceEpoch(end_ms), selection.name)

print('Done!')
