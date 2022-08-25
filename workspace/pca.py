import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler

#channels = [data.timeSeries('Sr88'), data.timeSeries('Ba138')]
channels = data.timeSeriesList(data.Input)
#groups = data.selectionGroupList(data.ReferenceMaterial | data.Sample)
groups = data.selectionGroupList(data.Sample)


cdl = []

for channel in channels:
	cd = np.array([])
	for group in groups:
		for sel in group.selections():
			d = channel.dataForSelection(sel)
			cd = np.append(cd, d)
			#cd = np.append(cd, d)
	cdl.append(cd)
			
scaler = MinMaxScaler()

					
nc = 5

A = np.column_stack(cdl)
A = A[~np.isnan(A).any(axis=1)]
scaler.fit(A)
A = scaler.transform(A)

pca = PCA(n_components=nc)
pca.fit(A)

print(pca.explained_variance_ratio_)

all = np.column_stack([c.data() for c in channels])
all[np.isnan(all).any(axis=1)] = 0
all = scaler.transform(all)
t = pca.transform(all)

for i in range(nc):
	pcai = np.copy(t[:, i])
	data.createTimeSeries(f'PCA{i}', data.Intermediate, data.timeSeries('TotalBeam').time(), pcai)

