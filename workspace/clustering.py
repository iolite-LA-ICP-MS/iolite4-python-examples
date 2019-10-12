from iolite.QtGui import QInputDialog
from sklearn.cluster import MeanShift, DBSCAN, OPTICS
from sklearn import preprocessing
import pandas as pd
import numpy as np

input_group_name = QInputDialog.getItem(None, "Clustering Example", "Group:", data.selectionGroupNames(), 0, False)

if not input_group_name:
	raise RuntimeError('No group supplied -- abort!')

group = data.selectionGroup(input_group_name)
channel_names = ['Co59_ppm', 'Mn55_ppm', 'Zn66_ppm', 'Rb85_ppm']
#channel_names = data.timeSeriesNames(data.Output)

df = pd.DataFrame()

# Collect all the data for the specified group:
for i, selection in enumerate(group.selections()):
    s = pd.Series(index=channel_names)
    s.name = selection.name
    if not s.name:
        s.name = group.name + str(i)

    for channel_name in channel_names:
       s[channel_name] = data.result(selection, data.timeSeries(channel_name)).value()

    df = df.append(s)

x = df.values
min_max_scalar = preprocessing.MinMaxScaler()
x_scaled = min_max_scalar.fit_transform(x)
df = pd.DataFrame(x_scaled)

#cc = MeanShift(cluster_all=False)
#cc = DBSCAN(eps=0.05, min_samples=5)
cc = OPTICS(metric='cityblock')
clustering = cc.fit(df)
n_clusters = len(np.unique(clustering.labels_))

# Create new groups with selections according to their cluster labels
for i, label in enumerate(clustering.labels_):
    cluster_group_name = input_group_name + '_cluster' + str(label)
    try:
        new_group = data.selectionGroup(cluster_group_name)
    except RuntimeError:
        new_group = data.createSelectionGroup(cluster_group_name, group.type)

    s = group.selection(i)
    new_s = data.createSelection(new_group, s.startTime, s.endTime, s.name)

print('Done!')
