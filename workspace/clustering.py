from iolite.QtGui import QInputDialog
from sklearn.cluster import MeanShift
import pandas as pd
import numpy as np

input_group_name = QInputDialog.getItem(None, "Clustering Example", "Group:", data.selectionGroupNames(), 0, False)

if not input_group_name:
	raise RuntimeError('No group supplied -- abort!')

group = data.selectionGroup(input_group_name)
channel_names = data.timeSeriesNames()

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

# Do mean-shift clustering:
ms = MeanShift()
clustering = ms.fit(df)
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
