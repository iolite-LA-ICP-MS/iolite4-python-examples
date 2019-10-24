# Note: depending on the number of channels and selections, this can take quite a while!

from iolite.QtGui import QInputDialog
from iolite.QtCore import QDateTime
from sklearn.cluster import MeanShift, DBSCAN, OPTICS, SpectralClustering
from sklearn import preprocessing
import pandas as pd
import numpy as np

# Run length encoding from stack overflow :) 
def rle(inarray):
	ia = np.asarray(inarray)
	n = len(ia)
	if n == 0: 
		return (None, None, None)
	else:
		y = np.array(ia[1:] != ia[:-1])     # pairwise unequal (string safe)
		i = np.append(np.where(y), n - 1)   # must include last element posi
		z = np.diff(np.append(-1, i))       # run lengths
		p = np.cumsum(np.append(0, z))[:-1] # positions
		return(z, p, ia[i])

input_group_name = None
input_group_name = QInputDialog.getItem(None, "Auto Selection Adjuster", "Group:", data.selectionGroupNames(), 0, False)

if not input_group_name:
	raise RuntimeError('No group supplied -- abort!')

group = data.selectionGroup(input_group_name)
output_group = data.createSelectionGroup(input_group_name + '_auto_adj', group.type)

# Specify which channels to use for the clustering:
#channel_names = ['Si29_ppm', 'Mg24_ppm', 'Al27_ppm', 'Fe57_ppm']
#channel_names = data.timeSeriesNames(data.Output)
channel_names = data.timeSeriesNames(data.Input)

for selection in group.selections():	
	# Create a data frame of the selection's data for channels specified above:
	d = {}    
	for channel_name in channel_names:
		d[channel_name] = data.timeSeries(channel_name).dataForSelection(selection)
	df = pd.DataFrame(d)

	# Scale the data? Not sure if this is required.
	x = df.values
	min_max_scalar = preprocessing.MinMaxScaler()
	x_scaled = min_max_scalar.fit_transform(x)
	df = pd.DataFrame(x_scaled)

	cc = MeanShift()
	clustering = cc.fit(df)
	n_clusters = len(np.unique(clustering.labels_))

	# Calculate longest stretches of same label
	z, p, _ = rle(clustering.labels_)
	max_rl = np.amax(z)
	start_index = p[np.argmax(z)]

	t = data.timeSeries(channel_names[0]).timeForSelection(selection)

	start_ms = int(1000*(t[start_index]))
	delta = np.median(np.diff(t))
	end_ms =  start_ms + max_rl * delta * 1000
	new_s = data.createSelection(output_group, QDateTime.fromMSecsSinceEpoch(start_ms), QDateTime.fromMSecsSinceEpoch(end_ms), selection.name)

print('Done!')
