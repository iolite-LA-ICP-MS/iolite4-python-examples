import numpy as np

for c in data.timeSeriesList(data.Input):
        if c.name == 'TotalBeam':
                continue

        di = np.diff(c.data())
        di = np.abs(di[~np.isnan(di)])
        di = di[np.nonzero(di)]
        guess = 1000/np.min(di)
        if guess != guess:
                guess = 10

        print('%s = %i'%(c.name, guess))
        c.setProperty('Dwell Time (ms)', guess)