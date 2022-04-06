#/ Type: DRS
#/ Name: 3D Trace Elements
#/ Authors: Joe Petrus and Bence Paul
#/ Description: Trace elements with multiple reference materials
#/ References: None
#/ Version: 1.1
#/ Contact: support@iolite-software.com

"""
Copyright (c) 2021 iolite software

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

from iolite.QtGui import QLabel, QAbstractTableModel, QTableView, QTabWidget, QWidget, QSortFilterProxyModel, QLineEdit
from iolite.QtGui import QToolButton, QMenu, QWidgetAction, QCheckBox,QAbstractItemView, QColor, QComboBox, QDoubleSpinBox, QColorDialog
from iolite.QtGui import QWidget, QVBoxLayout, QSizePolicy, QDialog, QHBoxLayout, QShortcut, QKeySequence, QGroupBox, QPen
from iolite.QtGui import QHeaderView, QSplitter, QPlainTextEdit, QDialogButtonBox, QMessageBox, QTableWidget, QTableWidgetItem
from iolite.QtGui import QInputDialog, QFileDialog, QFont, QBrush, QAction, QStyledItemDelegate, QFrame, QFormLayout, QSpinBox, QListWidget
from iolite.QtGui import QItemSelectionModel, QPushButton, QApplication, QMargins
from iolite.QtCore import Signal, QAbstractListModel, QModelIndex, QPoint, QFile, QIODevice, QTimer, QEvent, QSettings, QDir, QRegularExpression
from iolite.Qt import Qt
from iolite.QtUiTools import QUiLoader

from iolite.ui import IolitePlotPyInterface as Plot
from iolite.ui import Iolite3DPlotPyInterface as Plot3d
from iolite.ui import OverlayButton, QCPColorGradient, QCPErrorBars, PythonSyntaxHighlighter, OverlayMessage, QCPRange
from iolite.ui import CommonUIPyInterface as CUI
from iolite.types import Result

import numpy as np
import pandas as pd
import re
import sys
import itertools
import ast

from datetime import datetime
from math import sqrt, log, ceil
from functools import partial
from types import SimpleNamespace
from enum import Flag, auto

from scipy.interpolate import UnivariateSpline
from scipy.odr import Model, RealData, ODR
from scipy import stats
from scipy.signal import savgol_filter

import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std
from iolite_helpers import fitLine, formatResult

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# Replace this with one of the colorful QCPColorGradients?
import random
colors = {
    i: QColor(random.randrange(50, 200), random.randrange(50, 200), random.randrange(50, 200), 125)
    for i in range(100)
}

'''
Custom Exceptions to help track down issues
'''
class NoRMSelectionsError(RuntimeError):
    pass


class MissingRMGroupError(RuntimeError):
    pass


'''
Block functions and classes
'''

def linear(B, x):
    return B[0]*x + B[1]


def lineartz(B, x):
    return B[0]*x


def makeBeamSeconds():
    method = drs.setting('BeamSecondsMethod')
    channelName = drs.setting('BeamSecondsChannel')
    value = drs.setting('BeamSecondsValue')

    try:
        bs = data.timeSeries('BeamSeconds')
        if bs.property('BeamSecondsMethod') == method and bs.property('BeamSecondsChannel') == channelName and bs.property('BeamSecondsValue') == value:
            return
    except Exception as e:
        print(e)
        pass

    if 'log' in method.lower():
        drs.createBeamSecondsFromLaserLog()        
    elif 'gap' in method.lower():
        drs.createBeamSecondsFromSamples()
    elif 'cutoff' in method.lower():
        drs.createBeamSecondsFromCutoff(data.timeSeries(channelName), value)
    elif 'jump' in method.lower():
        drs.createBeamSecondsFromJump(data.timeSeries(channelName), value)
    else:
        print('Tried to make beam seconds with an unknown method')

    data.timeSeries('BeamSeconds').setProperty('BeamSecondsMethod', method)
    data.timeSeries('BeamSeconds').setProperty('BeamSecondsChannel', channelName)
    data.timeSeries('BeamSeconds').setProperty('BeamSecondsValue', value)


class Block(object):

    def __init__(self, selections, label='Unlabeled'):
        self.selections = selections
        self.label = label
        self.fits = {}
        self.lastDFHash = 0
        self.df = None

    def hash(self):
        shash = np.sum([s.midTimeInSec for s in self.selections])
        thash = np.sum([c.hash() for c in data.timeSeriesList(data.Intermediate, {'DRSType': 'BaselineSubtracted'}) if 'TotalBeam' not in c.name])
        drshash = drs.setting('NormalizeExternals') + np.sum([ord(a) for a in drs.setting('MasterExternal')])
        return shash+thash+drshash

    def fitHash(self, channel):
        c = data.timeSeries(channel)
        ftz = 0 if not c.property('FitThroughZero') else c.property('FitThroughZero')
        model = np.sum([ord(a) for a in c.property('Model')])
        return self.hash() + ftz + model

    def midTime(self):
        return np.mean([s.midTimeInSec for s in self.selections])

    def dataFrame(self):
        if self.df is None or self.lastDFHash != self.hash():
            self.updateDataFrame()

        return self.df

    def dataFrameForChannel(self, name):
        df = self.dataFrame()
        cols = [col for col in self.df.columns if col.startswith(name)]
        cols += ['sel_mid_time', 'sel_duration', 'group']
        df = df[cols]
        ext = data.timeSeries(name).property('External standard')
        if type(ext) == str :
            df = df[df['group'].isin(ext.split(','))]
        else:
            raise RuntimeError(f'No externals set for {name}')

        return df

    def updateDataFrame(self):
        self.df = pd.DataFrame()
        ryields, abdf = calculateRelativeYields()

        channelNames = [n for n in data.timeSeriesNames(data.Input) if 'TotalBeam' not in n]
        cpsChannels = [data.timeSeries(f'{c}_CPS') for c in channelNames]

        for sel in self.selections:
            uuid = sel.property('UUID')
            rmdata = data.referenceMaterialData(sel.group().name)

            # Original approach was to use the nanmedian of ablation factors:
            norm = ryields[sel.group().name] if drs.setting('NormalizeExternals') else 1
            # There can be noticable Tc dependent variation of af, so going to try fitting that...
            #abdf = abdf.dropna()
            #smx = sm.add_constant(abdf['Tc'], prepend=False)
            #smy = abdf[sel.group().name]
            #res = sm.RLM(smy, smx).fit()

            self.df.loc[uuid, 'sel_mid_time'] = float(sel.midTimeInSec)
            self.df.loc[uuid, 'sel_duration'] = float(sel.duration)
            self.df.loc[uuid, 'group'] = sel.group().name

            for name, channel in zip(channelNames, cpsChannels):
                try:
                    #norm = linear(res.params, data.elements[channel.property('Element')]['Tcond_Lodders'])
                    self.df.loc[uuid, name] = data.result(sel, channel).value()/norm
                    self.df.loc[uuid, f'{name}_Uncert'] = data.result(sel, channel).uncertaintyAs2SE()/(norm)
                    rmValue = rmdata[channel.property('Element')].valueInPPM()
                    rmUncert = rmdata[channel.property('Element')].uncertainty()
                    self.df.loc[uuid, f'{name}_RMppm'] = rmValue
                    self.df.loc[uuid, f'{name}_RMppm_Uncert'] = rmUncert if rmUncert else rmValue*0.02
                except:
                    self.df.loc[uuid, f'{name}_RMppm'] = np.nan
                    self.df.loc[uuid, f'{name}_RMppm_Uncert'] = np.nan
                    continue

        self.df = self.df.sort_values(by=['sel_mid_time'])
        self.lastDFHash = self.hash()

    def fit(self, name):
        if data.timeSeries(name).property('External standard') == 'Model':
            # Note: when the property is set with a dict from python
            # it is converted to QVariantMap, therefore the keys become strings
            # so need to use the string version of label to have it match!
            try:
                sens = data.timeSeries(name).property('ModelSensitivities')[str(self.label)]
                self.fits[name] = {
                    'slope': sens,
                    'slope_uncert': sens*0.01,
                    'intercept': 0.,
                    'intercept_uncert': 0.,
                    'hash': self.fitHash(name),
                    'sm_res': 0
                }
            except:
                self.fits[name] = {
                    'slope': np.nan,
                    'slope_uncert': np.nan,
                    'intercept': np.nan,
                    'intercept_uncert': np.nan,
                    'hash': self.fitHash(name),
                    'sm_res': np.nan
                }
        elif name not in self.fits:
            self.fits[name] = self.updateFit(name)

        elif not self.fits[name] is None and self.fits[name]['hash'] != self.fitHash(name):
            self.fits[name] = self.updateFit(name)

        return self.fits[name]

    def updateFit(self, name):
        channel = data.timeSeries(name)
        fitThroughZero = channel.property('FitThroughZero')
        model = channel.property('Model') if channel.property('Model') else 'ODR'
        df = self.dataFrameForChannel(name).copy()

        if not df[f'{channel.name}_RMppm'].notna().values.any():
            print(f'No data for {channel.name} for chosen RMs. No fit.')
            return

        df = df.dropna()
        func = lineartz if fitThroughZero or len(df['group'].unique()) == 1 else linear
        if len(df['group'].unique()) == 1:
            slope = df[name].mean()/df[f'{name}_RMppm'].mean()
            slope_uncert = (df[name].mean()/df[f'{name}_RMppm'].mean())*(df[f'{name}_Uncert'].mean()/df[name].mean())
            intercept = 0.
            intercept_uncert = 0.
            res = None
        else:
            smy = df[f'{name}']
            smx = df[f'{name}_RMppm'] if fitThroughZero else sm.add_constant(df[f'{name}_RMppm'], prepend=False)
            smw = 1/(df[f'{name}_Uncert']/df[name])**2

            if model == 'WLS':
                res = sm.WLS(smy, smx, weights=smw).fit()
            elif model == 'OLS':
                res = sm.OLS(smy, smx).fit()
            elif model == 'RLM':
                res = sm.RLM(smy, smx).fit()
            elif model == 'ODR':
                m = Model(func)
                md = RealData(df[f'{name}_RMppm'], df[f'{name}'], sx=df[f'{name}_RMppm_Uncert']/df[f'{name}_RMppm'], sy=df[f'{name}_Uncert']/df[name])
                b0 = [1000.] if fitThroughZero else [1000., 100.]
                od = ODR(md, m, beta0=b0)
                oo = od.run()
                #oo.pprint()
                res = SimpleNamespace()
                res.params = [oo.beta[0]] if fitThroughZero else [oo.beta[0], oo.beta[1]]
                res.bse = [oo.sd_beta[0]] if fitThroughZero else [oo.sd_beta[0], oo.sd_beta[1]]
            elif model == 'York':
                #fit = fitLine(df[f'{name}_RMppm'], df[f'{name}_RMppm_Uncert']/df[f'{name}_RMppm'], df[f'{name}'], df[f'{name}_Uncert']/df[name], np.zeros(len(df[name])))
                fit = fitLine(df[f'{name}_RMppm'], df[f'{name}_RMppm_Uncert'], df[f'{name}'], df[f'{name}_Uncert'], np.zeros(len(df[name])))
                res = SimpleNamespace()
                res.params = [fit['m'], fit['b']]
                res.bse = [fit['sigma_m'], fit['sigma_b']]

            slope = res.params[0]
            slope_uncert = res.bse[0]

            if fitThroughZero:
                intercept = 0.
                intercept_uncert = 0.
            else:
                intercept = res.params[1]
                intercept_uncert = res.bse[1]

        return {
            'slope': slope,
            'slope_uncert': slope_uncert,
            'intercept': 0. if fitThroughZero else intercept,
            'intercept_uncert': 0. if fitThroughZero else intercept_uncert,
            'hash': self.fitHash(name),
            'sm_res': res
        }

    def slope(self, name):
        if not self.fit(name) is None:
            return self.fit(name)['slope']
        else:
            return None

    def slopeUncert(self, name):
        return self.fit(name)['slope_uncert']

    def intercept(self, name):
        return self.fit(name)['intercept']

    def interceptUncert(self, name):
        return self.fit(name)['intercept_uncert']


def calculateRelativeYields():
    groupNames = data.selectionGroupNames(data.ReferenceMaterial)
    masterGroupName = drs.setting('MasterExternal')
    if not masterGroupName:
        masterGroupName = groupNames[0]
    masterGroup = data.selectionGroup(masterGroupName)

    cpsChannelNames = [n for n in data.timeSeriesNames(data.Intermediate, {'DRSType': 'BaselineSubtracted'}) if 'TotalBeam' not in n]
    cpsChannels = [data.timeSeries(c) for c in cpsChannelNames]

    M = np.empty( (len(cpsChannels), len(groupNames)) )
    M.fill(np.nan)

    for col, groupName in enumerate(groupNames):
        group = data.selectionGroup(groupName)
        for row, channel in enumerate(cpsChannels):
            channelElement = channel.property('Element')
            try:
                groupResult = data.groupResult(group, channel).value()
                groupRMValue = data.referenceMaterialData(groupName)[channelElement].valueInPPM()
                masterGroupResult = data.groupResult(masterGroup, channel).value()
                masterGroupRMValue = data.referenceMaterialData(masterGroupName)[channelElement].valueInPPM()
                M[row][col] = (groupResult/groupRMValue) / (masterGroupResult/masterGroupRMValue)
            except Exception as e:
                pass
                #print(f'Problem: {group.name} {channelElement} {groupResult} {groupRMValue} {masterGroupResult} {masterGroupRMValue}')

    # Todo: investigate whether there are any trends with the relative yield and mass or Tc?
    df = pd.DataFrame(M, columns=groupNames)
    df['Tc'] = [data.elements[c.property('Element')]['Tcond_Lodders'] for c in cpsChannels]
    df['Element'] = [c.property('Element') for c in cpsChannels]
    ablationFactors = dict(zip(groupNames, np.nanmedian(M, axis=0)))
    return ablationFactors, df

def assignExternalAffinities(sels=None):
    from sklearn.neighbors import NearestCentroid
    from sklearn.preprocessing import StandardScaler

    externalsInUse = list(set(list(itertools.chain(*[c.property('External standard').split(',') for c in data.timeSeriesList(data.Input) if 'TotalBeam' not in c.name and 'External standard' in c.properties().keys()]))))
    if sels is None:
        sels = list(itertools.chain(*[sg.selections() for sg in data.selectionGroupList(data.ReferenceMaterial | data.Sample)]))

    if len(externalsInUse) == 0:
        return externalsInUse
    if len(externalsInUse) == 1:
        [s.setProperty('External affinity', externalsInUse[0]) for s in sels]
        return externalsInUse

    isElements = list(set([s.property('Internal element') for s in sels]))
    affinityElements = list(set([s.property('Affinity elements') for s in sels]))

    print(f'assignExternalAffinities with {isElements} {affinityElements}')

    if len(isElements) == 0 or len(affinityElements) == 0:
        print('Could not update external affinities for selections...')
        return externalsInUse

    combinations = list(itertools.product(isElements, affinityElements))
    print(f'Combinations are {combinations}')

    def sumForSelection(sel, channels):
        n = len(data.timeSeries(channels[0]).dataForSelection(sel))
        s = np.zeros(n)
        for c in channels:
            s += data.timeSeries(f'{c}_CPS').dataForSelection(sel)

        return s

    aff = {}
    for comb in combinations:
        print(comb)
        ise = comb[0]
        afe = comb[1]
        if not ise or not afe or afe in externalsInUse:
            print(f"Continuing due to unset IntStd or Affinity not being specified")
            continue

        X = np.empty( (0, len(afe.split(','))) )
        y = np.array([])
        for gi, gn in enumerate(externalsInUse):
            g = data.selectionGroup(gn)
            for sel in g.selections():
                norm = sumForSelection(sel, ise.split(','))
                sd = np.column_stack( (data.timeSeries(f'{c}_CPS').dataForSelection(sel)/norm for c in afe.split(',')) )
                X = np.vstack( (X, sd) )
                gy = np.array([gi]*len(norm))
                y = np.concatenate( (y, gy) )

        scaler = StandardScaler().fit(X, y)
        X = scaler.transform(X)

        nc = NearestCentroid()
        nc.fit(X,y)
        aff[comb] = {'classifier': nc, 'scaler': scaler}

    for sel in sels:
        ise = sel.property('Internal element')
        afe = sel.property('Affinity elements')
        if not ise or not afe:
            continue

        if afe in externalsInUse:
            # If a selection has a rm for affinity elements rather than a list of elements
            # set it to have that external affinity without using the classifier.
            sel.setProperty('External affinity', afe)
        else:
            channels = afe.split(',')
            norm = sumForSelection(sel, ise.split(','))
            seld = np.array([np.median(data.timeSeries(f'{c}_CPS').dataForSelection(sel)/norm) for c in channels]).reshape(1, -1)
            seld = np.nan_to_num(aff[(ise, afe)]['scaler'].transform(seld))
            i = int(aff[(ise,afe)]['classifier'].predict(seld)[0])
            sel.setProperty('External affinity', externalsInUse[i])

    return externalsInUse


def findBlocks(method=None):
    '''
    Have 4 modes:
        1. Assigned - Only use selections that have been assigned explicitly (fastest)
        2. Simple - Use the fast method looking at selection time jumps to find new blocks (faster)
        3. Clustering - Use clustering with a specific number of blocks (fast)
        4. Auto Clustering - Use clustering that searches for the *best* number of blocks (slow)

    We also want to be able to use one of the more automatic methods with some (not all) selections being specified.
    '''
    externalsInUse = set(list(itertools.chain(*[c.property('External standard').split(',') for c in data.timeSeriesList(data.Input) if 'TotalBeam' not in c.name and type(c.property('External standard')) == str])))
    try:
        externalsInUse.remove('Model')
    except:
        pass

    # the user may have selected an external that doesn't exist. If so, it will raise a RuntimeError.
    # That exception should be caught in calling functions
    try:
        groups = [data.selectionGroup(ext) for ext in externalsInUse if ext]
    except RuntimeError:
        raise MissingRMGroupError

    selections = list(itertools.chain(*[sg.selections() for sg in groups]))
    selections.sort(key=lambda s: s.midTimeInSec)
    selMidTimes = [s.midTimeInSec for s in selections]

    def assignedBlock(s):
        if s.property('Block') is not None:
            return s.property('Block')

        return -1

    assignedBlocks = np.array([assignedBlock(s) for s in selections])
    allAssigned = np.all(assignedBlocks>=0)

    diffs = np.column_stack( (range(len(selMidTimes)), np.insert(np.diff(selMidTimes), 0, 0)))
    cutoff = np.mean(diffs[:, 1])

    if not method:
        method = 'Simple'

    if method == 'Simple' and not allAssigned:
        current_label = 1
        labels = []
        for i, v in enumerate(diffs[:, 1]):
            if v > cutoff:
                current_label += 1
            labels.append(current_label)
    elif method == 'Auto Clustering' and not allAssigned:
        from sklearn.cluster import KMeans
        from sklearn.metrics import silhouette_score
        scores = []
        a = np.array(selMidTimes).reshape(-1,1)
        for nc in range(2, len(selMidTimes)):
            km = KMeans(n_clusters=nc).fit(a)
            ss = silhouette_score(a, km.labels_)
            scores.append(ss)

        nc = np.argmax(scores)+2 # 1 (bc above starts at nc=2) + 1 (because starting at 0)
        print(f'findBlocks decided to use {nc} clusters')
        km = KMeans(n_clusters=int(nc)).fit(a)
        labels = km.labels_ + 1
    elif method == 'Clustering' and not allAssigned:
        from sklearn.cluster import KMeans
        nc = drs.setting('NClusters')
        if not nc:
            nc = 5
        print(f'findBlocks decided to use {nc} clusters')
        a = np.array(selMidTimes).reshape(-1,1)
        km = KMeans(n_clusters=int(nc)).fit(a)
        labels = km.labels_ + 1
    else:
        labels = assignedBlocks

    block_selections = {}

    print(f'Used method = {method}, got labels = {labels}')
    specified = [si for si, s in enumerate(selections) if s.property('Block') is not None]
    for k in np.unique(labels):
        matches = np.where(np.array(labels).astype(int) == k)[0]        
        ind = [i for i in matches if i not in specified]
        block_selections[k] = list(np.array(selections)[ind])

    blocks = []

    for k in block_selections:
        specified = [s for si, s in enumerate(selections) if s.property('Block') == k]
        if len(block_selections[k] + specified) == 0:
            continue
        blocks.append(Block(block_selections[k] + specified, k))

    blocks.sort(key=lambda b: b.midTime())
    # Relabel according to time order:
    for i, block in enumerate(blocks):
        block.label = i+1

    return blocks


def fitSurface(blocks, channelName):

    slopes = np.array([block.slope(channelName) for block in blocks])
    slopes_unc = np.array([block.slopeUncert(channelName) for block in blocks])
    intercepts = np.array([block.intercept(channelName) for block in blocks])
    intercepts_unc = np.array([block.interceptUncert(channelName) for block in blocks])
    times = np.array([block.midTime() for block in blocks])

    cpsChannel = data.timeSeries(f'{channelName}_CPS')
    splineType = drs.setting('SplineType')

    if len(blocks) < 2 and splineType != "MeanMean":
        IoLog.information("Changing spline type to MeanMean due to fewer than 2 blocks found")
        splineType = 'MeanMean'
    elif len(blocks) < 3 and splineType != "StepLinear":
        IoLog.information("Changing spline type to StepLinear due to fewer than 3 blocks found")
        splineType = 'StepLinear'

    if np.all(np.isnan(slopes)):
        return None, None

    slopes_unc[slopes_unc == 0] = np.nanmean(slopes)*0.05
    intercepts_unc[intercepts_unc == 0] = 1e-6

    slope_spl = data.spline(times, slopes, slopes_unc, splineType, cpsChannel.time())
    intercept_spl = data.spline(times, intercepts, intercepts_unc, splineType, cpsChannel.time())

    data.createTimeSeries(f'{channelName}_slope', data.Intermediate, cpsChannel.time(), slope_spl, {'DRS': '3D Trace Elements'})
    data.createTimeSeries(f'{channelName}_intercept', data.Intermediate, cpsChannel.time(), intercept_spl, {'DRS': '3D Trace Elements'})

    def surface(t, c):
        i = np.searchsorted(cpsChannel.time(), t)
        m = slope_spl[i]
        b = intercept_spl[i]
        return m*c + b

    def surfaceInv(t, I):
        i = np.searchsorted(cpsChannel.time(), t)
        m = slope_spl[i]
        b = intercept_spl[i]
        # I = m*c + b, so c = (I - b)/m
        return (I - b)/m

    return surface, surfaceInv


class Calibration(object):

    normal = 0
    inverse = 1

    def __init__(self):
        self.surfaces = {}
        self.blocks = []
        self.frac = {}

    def block(self, bn):
        if bn >= len(self.blocks):
            raise RuntimeError(f'Block number {bn} out of range {len(self.blocks)}')

        return self.blocks[bn]

    def updateBlocks(self):
        self.blocks = findBlocks(drs.setting('BlockFindingMethod'))

    def surface(self, name, update=False, inv=False):
        if name not in self.surfaces or update:
            self.updateSurface(name)

        s = self.surfaces[name][self.inverse] if inv else self.surfaces[name][self.normal]
        return s

    def updateSurface(self, name):
        self.surfaces[name] = fitSurface(self.blocks, name)

    def semiquant(self, name):
        cps = data.timeSeries(f'{name}_CPS')
        return self.surface(name, inv=True)(cps.time(), cps.data())

    def clearFractionationCache(self):
        self.frac = {}

    def clearSurfaceCache(self):
        self.surfaces = {}

    def fractionation(self, name, update=False):
        if name not in self.frac or update:
            print(f"Updating fractionation for {name}")
            self.updateFractionation(name)

        return self.frac[name]

    def updateFractionation(self, name):
        channel = data.timeSeries(name)

        try:
            externals = [es for es in channel.property('External standard').split(',') if es]
        except:
            self.frac[name] = pd.DataFrame()
            return

        if not externals:
            self.frac[name] = pd.DataFrame()
            return

        groups = [data.selectionGroup(ext) for ext in externals]
        cpsChannel = data.timeSeries(f'{name}_CPS')

        try:
            makeBeamSeconds()
            bs = data.timeSeries('BeamSeconds')
        except Exception as e:
            print(e)
            print('There was a problem making beam seconds and therefore fractionation cannot be determined.')
            return

        allSels = [[s for s in sg.selections()] for sg in data.selectionGroupList(data.ReferenceMaterial | data.Sample)]
        allSels = list(itertools.chain.from_iterable(allSels))
        allIS = [sel.property('Internal element') for sel in allSels]
        isElementsList = [ie for ie in set(allIS) if ie != 'None' and ie != '' and ie]
        isElementsList.sort()

        def sumsForSel(sel, isElements):
            if not isElements:
                raise RuntimeError('No internals set')

            selInd = cpsChannel.selectionIndices(sel)
            selPPMSum = np.zeros(len(bs.timeForSelection(sel)))
            rmPPMSum = 0
            for el in [el for el in isElements.split(',') if el]:
                selPPMSum += self.semiquant(el)[selInd]
                rmPPMSum += data.referenceMaterialData(sel.group().name)[data.timeSeries(f'{el}_CPS').property('Element')].valueInPPM()
            return selPPMSum, rmPPMSum

        fdf = pd.DataFrame()

        for sg in groups: # Loop through each external for this channel
            for sel in sg.selections(): # Loop through each selection of each group
                for isElements in isElementsList: # Collect ratio data for each of the IS combinations in use
                    try:
                        selPPMSum, rmPPMSum = sumsForSel(sel, isElements)
                        thisPPM = data.referenceMaterialData(sg.name)[cpsChannel.property('Element')].valueInPPM()

                        t = bs.dataForSelection(sel)
                        selInd = cpsChannel.selectionIndices(sel)
                        sq = self.semiquant(name)[selInd]
                        r = (sq/selPPMSum)*(rmPPMSum/thisPPM)

                        if isElements == name:
                            r = np.ones(len(r))

                        df = pd.DataFrame({'t': t, 'r': r, 'IS': isElements, 'group': sg.name}, index=pd.Series([pd.Timedelta(milliseconds=int(tt*1000)) for tt in t]))
                        fdf = fdf.append(df) if len(fdf) > 0 else df
                    except Exception as e:
                        continue

        self.frac[name] = fdf

    def fitFractionation(self, name, isElements=None, td=None, k=None, group=None):
        ft = data.timeSeries(name).property('FractionationFitType')
        fc = data.timeSeries(name).property('FractionationCorrection')

        if not k:            
            k = 1 if not ft or ft == 'Linear' else 3

        if not isElements:
            isdf = self.fractionation(name)
        else:
            fdf = self.fractionation(name)
            isdf = fdf[fdf['IS'] == isElements]

        if group is not None:
            isdf = isdf[isdf['group'] == group]

        if not td:
            # Aim for 30 points
            td = (isdf.index.max()-isdf.index.min())/30
        elif isinstance(td, str):
            td = pd.Timedelta(td)

        if not td or len(isdf) == 0:
            print(f'Could not fit fractionation for {name} {isElements} {td} {k} {group}')
            return None, None, None, None

        rsd = isdf.resample(td).sem()['r']

        isdf = isdf.resample(td).median()
        t = isdf['t']
        r = isdf['r']

        rsd[rsd==0] = 0.02*np.nanmean(r)
        rsd[rsd!=rsd] = 0.02*np.nanmean(r)
        rsd[rsd<0.001*np.nanmean(r)] = 0.02*np.nanmean(r)
        rsd[rsd>1] = 0.02*np.nanmean(r)

        def ones(x):
            try:
                return np.ones(len(x))
            except:
                return np.ones(1)

        if fc and name != isElements:
            sx = np.linspace(0, t.max(), 100)
            if k == 3:
                spline = UnivariateSpline(t[1:-1], r[1:-1], w=1/rsd[1:-1], k=k, s=len(t)*2)
            else:
                spline = UnivariateSpline(t[1:-1], r[1:-1], w=1/rsd[1:-1], k=k, s=1e9)
        elif fc and name == isElements:
            spline = ones
        else:
            spline = None

        return t, r, rsd, spline


#    def materialFractionation(self, selection, masterExt):
#        '''
#        Want to calculate:
#            (master_IS_SQ/selection_channel_SQ) * (selection_IS_ppm / master_IS_ppm)
#        '''

#        try:
#            bs = data.timeSeries('BeamSeconds')
#        except:
#            drs.createBeamSecondsFromLaserLog()
#            bs = data.timeSeries('BeamSeconds')

#        isChannels = selection.property('Internal element').split(',')
#        isValue = selection.property('Internal value') # TODO: get units and convert this to ppm if it isn't
#        t = data.timeSeries('TotalBeam').time()

#        mdf = pd.DataFrame()
#        masterPPM = np.sum([data.referenceMaterialData(masterExt)[data.timeSeries(channelName).property('Element')].valueInPPM() for channelName in isChannels])

#        sq = {}
#        for channelName in isChannels:
#            sq[channelName] = self.semiquant(channelName)

#        for channelName in isChannels:
#            dht, dhd = data.compileDownholeFromArray(data.selectionGroup(masterExt), sq[channelName])
#            df = pd.DataFrame({'t': dht, f'Master_{channelName}':dhd}, index=pd.Series([pd.Timedelta(milliseconds=int(tt*1000)) for tt in dht]))
#            mdf = mdf.append(df) if len(mdf) > 0 else df

#        selectionSQ = np.zeros(len(t))
#        ind = data.timeSeries('TotalBeam').selectionIndices(selection)
#        selt = data.timeSeries('BeamSeconds').data()[ind]
#        for channelName in isChannels:
#            df = pd.DataFrame({'t': selt, f'Selection_{channelName}': sq[channelName][ind]}, index=pd.Series([pd.Timedelta(milliseconds=int(tt*1000)) for tt in selt]))
#            mdf = mdf.append(df)


#        td = (mdf.index.max()-mdf.index.min())/30
#        sdf = mdf.resample(td).median()
#        sdf = sdf[1:-1]
#        sdf = sdf.dropna()

#        masterSQ = sdf[ [f'Master_{channelName}' for channelName in isChannels]].sum(axis=1)
#        selectionSQ = sdf[ [f'Selection_{channelName}' for channelName in isChannels]].sum(axis=1)
#        t = sdf['t']
#        #fd = (masterSQ/selectionSQ) * (isValue/masterPPM)
#        fd = (selectionSQ/masterSQ) * (masterPPM/isValue)

#        fit = np.polyfit(t, fd, 1)
#        return (t, fd), partial(linear, (fit[0], fit[1]))

#        #return (masterSQ/selectionSQ) * (isValue/masterPPM)

    def uncertainty(self, name):
        return np.mean([block.slopeUncert(name)/block.slope(name) for block in self.blocks])


def runDRS():
    drs.message("Starting 3D Trace Elements DRS...")
    drs.progress(0)
    drs.setProperty('isRunning', True)

    commonProps = {'DRS': drs.name()}

    # Get settings
    settings = drs.settings()

    print(settings)
    drs.clearSelectionProperties(QRegularExpression("Sensitivity.+"))

    indexChannel = data.timeSeries(settings["IndexChannel"])
    drs.setIndexChannel(indexChannel)

    # Setup index time
    drs.message("Setting up index time...")
    drs.progress(5)
    drs.setIndexChannel(indexChannel)

    try:
        blGrp = data.selectionGroupList(data.Baseline)[0]
    except:
        IoLog.error("There are no baseline groups. 3D Trace Elements cannot proceed...")
        drs.message("Error. See Messages")
        drs.progress(100)
        drs.finished()
        return

    if len(blGrp.selections()) < 1:
        IoLog.error("No baseline selections. Please select some baselines. 3D Trace Elements cannot proceed...")
        drs.message("Error. See Messages")
        drs.progress(100)
        drs.finished()
        return

    # Setup the mask
    maskOption = settings["Mask"]
    maskChannel = data.timeSeries(settings["MaskChannel"])
    maskMethod = settings['MaskMethod']
    cutoff = settings["MaskCutoff"]
    trim = settings["MaskTrim"]

    if maskOption:
        drs.message("Making mask...")
        drs.progress(10)
        if 'Laser' in maskMethod:
            mask = drs.createMaskFromLaserLog(trim)
        else:
            mask = drs.createMaskFromCutoff(maskChannel, cutoff, trim)
        data.createTimeSeries('mask', data.Intermediate, indexChannel.time(), mask)
    else:
        mask = np.ones_like(indexChannel.data())
        data.createTimeSeries('mask', data.Intermediate, indexChannel.time(), mask)

    # Setup total Pb
    if len(data.timeSeriesNames(data.Input, {'Element': 'Pb'})) >= 3:
        totalPb = np.zeros(len(indexChannel.time()))

        for pb in [c for c in data.timeSeriesList(data.Input, {'Element': 'Pb'}) if 'Total' not in c.name]:
            totalPb += pb.data()

        totalPbChannel = data.createTimeSeries('PbTotal', data.Input, indexChannel.time(), totalPb, {**commonProps, 'Element': 'Pb', 'Mass': 208})
        totalPbChannel.setProperty('External standard', pb.property('External standard'))
        totalPbChannel.setProperty('FitThroughZero', pb.property('FitThroughZero'))
        totalPbChannel.setProperty('Model', pb.property('Model'))

    # Baseline Subtraction
    drs.baselineSubtract(blGrp, data.timeSeriesList(data.Input), mask, 10, 20)

    # Check that some Externals have been selected:
    atLeastOneExt = False
    for channel in data.timeSeriesList(data.Input):
        if type(channel.property("External standard")) == str:
            atLeastOneExt = True
            break

    if not atLeastOneExt:
        IoLog.error(f"No external standards are set. Please set at least one RM before continuing.")
        drs.message.emit("DRS did not finish.")
        drs.progress.emit(100)
        drs.finished.emit()
        return

    # Find blocks
    drs.message.emit('Finding blocks')
    drs.progress.emit(23)
    cal = Calibration()
    cal.updateBlocks()

    # Calculate SQ channels
    for ii, input in enumerate(data.timeSeriesList(data.Input)):
        if 'TotalBeam' in input.name:
            continue

        drs.progress.emit(25 + 25*float(ii)/len(data.timeSeriesList(data.Input)))
        drs.message.emit(f'Applying surface for {input.name}')

        try:
            surface = cal.surface(input.name, inv=True)
        except:
            IoLog.warning(f"There was an issue calculating the surface for {input.name}")
            continue

        if not surface:
            IoLog.warning(f'No surface for {input.name}')
            continue

        cps = data.timeSeries(f'{input.name}_CPS')
        ppm = surface(cps.time(), cps.data())

        props = { **commonProps,
            'Element': input.property('Element'),
            'Mass': input.property('Mass'),
            'Reference Material': input.property('External standard'),
            'Reference Material uncertainty': cal.uncertainty(input.name), # This needs some thinking....
            'Units': 'µg.g-1'
        }

        data.createTimeSeries(f'{input.name}_ppm', data.Output, indexChannel.time(), ppm, props)

    sels = list(itertools.chain(*[sg.selections() for sg in data.selectionGroupList(data.ReferenceMaterial | data.Sample)]))

    mfc = np.ones(len(indexChannel.time()))

    externalsInUse = assignExternalAffinities()
    affIndex = data.createTimeSeriesFromMetadata('ExtAffinityIndex', 'External affinity')

    # Optionally apply a fractionation correction
    #
    # TODO: Should only do this if we have internal standards set, AND not sum normalisation... ?
    #
    if np.any([c.property('FractionationCorrection') for c in data.timeSeriesList(data.Input)]):
        print('Attempting fractionation correction...')
        allSels = [[s for s in sg.selections()] for sg in data.selectionGroupList(data.ReferenceMaterial | data.Sample)]
        allSels = list(itertools.chain.from_iterable(allSels))
        allIS = [sel.property('Internal element') for sel in allSels]
        isElementsList = list(set(allIS))
        isElementsList.sort()

        params = {}
        for c in [c for c in data.timeSeriesList(data.Input) if 'TotalBeam' not in c.name]:
            params[c.name] = cal.fractionation(c.name)

        bs = data.timeSeries('BeamSeconds').data()
        isIndex = data.createTimeSeriesFromMetadata('ISElementIndex', 'Internal element')

        for c in [c for c in data.timeSeriesList(data.Input) if 'TotalBeam' not in c.name]:
            fdf = params[c.name]
            if len(fdf) == 0:
                print(f'No frac for {c.name}')
                continue

            for isi, ise in enumerate(isElementsList):
                for gi, g in enumerate(externalsInUse):
                    t, r, rsd, sp = cal.fitFractionation(c.name, ise, group=g)
                    if not sp:
                        continue
                    ind = np.where( (isIndex.data() == isi) & (affIndex.data() == gi))[0]
                    ppm = data.timeSeries(f'{c.name}_ppm')
                    ppmd = ppm.data()
                    ppmd[ind] = ppmd[ind]/(sp(bs[ind]))
                    ppm.setData(ppmd)

    if False:
        # This is a bit of code for experimenting with forcing certain elements
        # in unknown samples to be homogeneous prior to normalization.
        # This can flatten out any remaining DH fractionation in
        # unknowns, but there may still be fractionation for other reasons
        for ppmc in data.timeSeriesList(data.Output, {'Units': 'µg.g-1'}):
            if ppmc.name not in ['Al27_ppm', 'Ca43_ppm', 'Fe57_ppm', 'Mg25_ppm']:
                continue
            for si, sel in enumerate(sels):
                selIS = data.timeSeries(sel.property('Internal element')+'_ppm').dataForSelection(sel)
                selIS = savgol_filter(selIS, 21, 3)
                selC = ppmc.dataForSelection(sel)
                selC = savgol_filter(selC, 21, 3)
                r = (selC/np.nanmedian(selC))/(selIS/np.nanmedian(selIS))
                ppmcd = ppmc.data()
                ppmcd[ppmc.selectionIndices(sel)] /= r
                ppmc.setData(ppmcd)

    # Calculate FQ
    if bool(settings['UseIntStds']):

        if not np.any([bool(sel.property('Internal value')) for sel in sels]) and not np.any([sel.property('Internal element') == 'Criteria' for sel in sels]):
            IoLog.error("No internal standard values found. Have you set any?")
            drs.message.emit("DRS did not finish.")
            drs.progress.emit(100)
            drs.finished.emit()
            return

        data.createTimeSeriesFromMetadata('ISValue', 'Internal value', sels)

        norm = np.ones(len(indexChannel.time()))
        crit_index = np.empty(len(indexChannel.time()))
        crit_index[:] = np.nan

        criteria_globals = {'where': np.where}
        for c in data.timeSeriesList(data.Input):
            try:
                criteria_globals[c.name] = data.timeSeries(f'{c.name}_CPS').data()
            except:
                pass

        try:
            criteria_list = settings['ISCriteria']
            for i in range(len(criteria_list)):
                criteria_list[i]['criteria'] = eval(f"lambda: where({criteria_list[i]['criteria']})", criteria_globals)
                criteria_list[i]['indicies'] = criteria_list[i]['criteria']()
        except Exception as e:
            pass
            #print(f'Not using criteria or there was a problem parsing the criteria...')

        ppmd = lambda s: data.timeSeries(f'{s}_ppm').data()

        for ii, sel in enumerate(sels):
            isvalue = sel.property('Internal value')
            iselement = sel.property('Internal element')
            isunit = sel.property('Internal units')
            si = indexChannel.selectionIndices(sel)

            if ii%20 == 0:
                drs.progress.emit(50 + 40*float(ii)/len(sels))
                drs.message.emit('Applying internal standards')

            if iselement == 'Criteria':
                for ci, criteria in enumerate(criteria_list):
                    csi = np.intersect1d(si, criteria['indicies'])
                    if len(csi) == 0:
                        continue
                    sum = np.zeros(len(csi))

                    for analyte in criteria['analytes'].split(','):
                        if not analyte:
                            continue
                        f = 1
                        if criteria['oxides']:
                            el = data.timeSeries(analyte).property('Element')
                            if el in ','.join(criteria['oxide_forms']):
                                matches = [f for f in criteria['oxide_forms'] if el in f]
                                if matches:
                                    f = data.oxideToElementFactor(matches[0])
                            else:
                                f = data.oxideToElementFactor(el)

                        sum += ppmd(analyte)[csi]*f

                    norm[csi] = (criteria['value']*1e4)/sum
                    crit_index[csi] = ci

            else:
                if not isvalue:
                    continue
                sum = np.zeros(len(si))
                for el in iselement.split(','):
                    if 'oxide' in isunit:
                        actual_el = data.timeSeries(el).property('Element')
                        try:
                            oxide_forms = sel.property('OxideForms')
                            matches = [f for f in oxide_forms.split(',') if actual_el in f]
                            f = data.oxideToElementFactor(matches[0])
                        except:
                            f = data.oxideToElementFactor(actual_el)
                    else:
                        f = 1

                    f *= 0.0001 if 'wtpc' in isunit else 1
                    sum += data.timeSeries(f'{el}_ppm').data()[si]*f

                norm[si] = float(isvalue)/(sum)

        # Replace parts that end up inf with 1
        norm[np.abs(norm) == np.inf] = 1

        data.createTimeSeries('Norm', data.Intermediate, indexChannel.time(), norm, commonProps)
        data.createTimeSeries('CriteriaIndex', data.Intermediate, indexChannel.time(), crit_index, commonProps)

        for c in data.timeSeriesList(data.Output):
            d = c.data()*norm
            c.setData(d)


    data.updateResults()

    # This bit of code uses the residual of its nearest RM to apply an additional correction
    # Todo: make it configurable
    if drs.setting('AffinityCorrection'):
        drs.message.emit('Doing affinity correction')
        drs.progress.emit(92)
        extFactors = {}
        for ext in externalsInUse:
            extFactors[ext] = {}
            for ppmc in data.timeSeriesList(data.Output, {'Units': 'µg.g-1'}):
                try:
                    rm = data.referenceMaterialData(ext)[ppmc.property('Element')].valueInPPM()
                    meas = data.groupResult(data.selectionGroup(ext), ppmc).value()
                    extFactors[ext][ppmc.name] = meas/rm
                    if abs(extFactors[ext][ppmc.name]-1) > 0.15:
                        print(f'Not going to correct {ppmc.name} for {ext} due to its large relative difference from the accepted value')
                        extFactors[ext][ppmc.name] = 1.
                        continue
                except Exception as e:
                    print(e)

        print(extFactors)

        for sel in sels:
            if sel.group().name in externalsInUse:
                continue
            for ppmc in data.timeSeriesList(data.Output, {'Units': 'µg.g-1'}):
                d = ppmc.data()
                try:
                    f = extFactors[sel.property('External affinity')][ppmc.name]
                    ind = ppmc.selectionIndices(sel)
                    d[ind] =  d[ind]/f
                except Exception as e:
                    print(e)
                ppmc.setData(d)

    data.updateResults()

    # Store sensitivities so LODs can be determined by results manager
    drs.message.emit('Storing sensitivities')
    drs.progress.emit(95)
    badSels = []
    for sel in sels:
        badChannels = []
        for cps in data.timeSeriesList(data.Intermediate):
            if 'CPS' not in cps.name or 'TotalBeam' in cps.name:
                continue
            try:
                ppm = data.timeSeries(cps.name.replace('CPS', 'ppm'))
                s = data.result(sel, cps).value()/data.result(sel, ppm).value()
                sel.setProperty('Sensitivity %s'%(cps.name.replace('_CPS', '')), s)
            except:
                badChannels.append(cps.name.replace('CPS', 'ppm'))
                continue

        badSels.append(sel)

    badChannelsString = ', '.join(badChannels)
    badChannelsString = badChannelsString if len(badChannelsString) < 25 else badChannelsString[0:23]+'...'
    badSelsString = ', '.join([s.name for s in badSels])
    badSelsString = badSelsString if len(badSelsString) < 25 else badSelsString[0:23]+'...'
    IoLog.warning(f'There was a problem calculating the sensitivity of {badChannelsString} for selection(s) {badSelsString}')

    # Need to update results again to get LODs calculated
    data.updateResults()
    drs.message.emit("Finished!")
    drs.progress.emit(100)
    drs.finished.emit()


'''
GUI-related classes
'''
class ExternalsModel(QAbstractTableModel):

    throughZeroChanged = Signal()
    fractionationChanged = Signal()

    def __init__(self, parent):
        super().__init__(parent)
        self.refreshChannels()

    def refreshChannels(self):
        self.beginResetModel()
        self.channels = [c for c in data.timeSeriesList(data.Input) if 'TotalBeam' not in c.name]

        try:
            drs.baselineSubtract(data.selectionGroupList(data.Baseline)[0], data.timeSeriesList(data.Input), None, 0, 0)
            for c in data.timeSeriesList(data.Input):
                if not c.property('Model'):
                    c.setProperty('Model', 'ODR')
        except:
            print('You must import data and create selections before using the 3DTE DRS.')


        self.endResetModel()

    def updateData(self):
        self.dataChanged.emit(self.index(0, 0), self.index(self.rowCount()-1, 1))

    def rowCount(self, index=QModelIndex()):
        return len(self.channels)

    def columnCount(self, index):
        return 5

    def headerData(self, section, orientation, role):
        if role != Qt.DisplayRole:
            return None

        if orientation == Qt.Horizontal:
            if section == 0:
                return 'Channel'
            elif section == 1:
                return 'External standards'
            elif section == 2:
                return 'Model'
            elif section == 3:
                return 'Zero'
            elif section == 4:
                return 'Frac. correction'

        return None

    def data(self, index, role):
        if not index.isValid():
            return None

        if role == Qt.UserRole:
            return self.channels[index.row()]

        if role == Qt.CheckStateRole and index.column() >= 3:
            if index.column() == 3 and self.channels[index.row()].property('FitThroughZero'):
                return Qt.Checked
            elif index.column() == 4 and self.channels[index.row()].property('FractionationCorrection'):
                return Qt.Checked
            return Qt.Unchecked

        if role != Qt.DisplayRole:
            return None

        if index.column() == 0:
            return self.channels[index.row()].name
        elif index.column() == 1:        
            return self.channels[index.row()].property('External standard')
        elif index.column() == 2:
            return self.channels[index.row()].property('Model')
        elif index.column() == 4:
            ft = self.channels[index.row()].property('FractionationFitType')
            if not ft:
                ft = 'None'
            return ft

        return None

    def setData(self, index, value, role = Qt.EditRole):
        if role == Qt.CheckStateRole and index.column() == 3:
            self.channels[index.row()].setProperty('FitThroughZero', value == Qt.Checked)
            self.throughZeroChanged.emit()
            self.dataChanged.emit(index, index)
        if role == Qt.CheckStateRole and index.column() == 4:
            self.channels[index.row()].setProperty('FractionationCorrection', value == Qt.Checked)
            ft = self.channels[index.row()].property('FractionationFitType')
            if value == Qt.Checked and (not ft or ft == 'None'):
                self.channels[index.row()].setProperty('FractionationFitType', 'Linear')
            elif value == Qt.Unchecked:
                self.channels[index.row()].setProperty('FractionationFitType', 'None')

            self.dataChanged.emit(index, index)
            self.fractionationChanged.emit()

    def flags(self, index):
        if index.column() in [1, 2, 3, 4]:
            return QAbstractTableModel.flags(self, index) | Qt.ItemIsEditable | Qt.ItemIsUserCheckable
        return QAbstractTableModel.flags(self, index)

class ExternalsDelegate(QStyledItemDelegate):

    def __init__(self, parent=None):
        QStyledItemDelegate.__init__(self, parent)

    def createEditor(self, parent, option, index):
        if index.column() == 1:
            return ReferenceMaterialComboBox(parent, index.data(Qt.UserRole))
        elif index.column() == 2:
            cb = QComboBox(parent)
            cb.addItems(['ODR', 'OLS', 'WLS', 'RLM', 'York'])
            return cb
        elif index.column() == 4:
            cb = QComboBox(parent)
            cb.addItems(['None', 'Linear', 'Spline'])
            return cb

        return QStyledItemDelegate.createEditor(self, parent, option, index)

    def setEditorData(self, editor, index):
        if index.column() == 1:
            editor.clear()
            editor.addItem(index.data(Qt.UserRole).property('External standard'))
            editor.currentText = index.data(Qt.UserRole).property('External standard')
        elif index.column() == 2:
            editor.currentText = index.data(Qt.UserRole).property('Model')
        elif index.column() == 4:
            editor.currentText = index.data(Qt.UserRole).property('FractionationFitType')

    def setModelData(self, editor, model, index):
        if index.column() == 2:
            index.data(Qt.UserRole).setProperty('Model', editor.currentText)
        elif index.column() == 4:
            index.data(Qt.UserRole).setProperty('FractionationCorrection', 'None' != editor.currentText)
            index.data(Qt.UserRole).setProperty('FractionationFitType', editor.currentText)
        try:
            model.sourceModel.dataChanged.emit(index, index)
        except:
            model.dataChanged.emit(index, index)


class ReferenceMaterialComboBox(QComboBox):

    def __init__(self, parent, channel):
        QComboBox.__init__(self, parent)
        self.channel = channel

    def updateText(self):
        self.clear()
        self.addItem(self.channel.property('External standard'))
        self.currentText = self.channel.property('External standard')

    def showPopup(self):
        menu = ReferenceMaterialsMenu(self)
        menu.activeChannels = [self.channel.name]
        menu.rmsForActiveChannels = []
        if type(self.channel.property('External standard')) == str:
            menu.rmsForActiveChannels = self.channel.property('External standard').split(',')
        menu.rmsChanged.connect(self.updateText)
        menu.setFixedWidth(self.width)
        p = self.pos
        p += QPoint(0, self.height)
        menu.exec_(self.parent().mapToGlobal(p))

class ReferenceMaterialsMenu(QMenu):

    rmsChanged = Signal()

    def __init__(self, parent):
        super().__init__(parent)
        self.activeChannels = []
        self.rmsForActiveChannels = []

        for rm in data.referenceMaterialNames():
            a = QWidgetAction(self)
            cb = QCheckBox(rm, self)
            cb.setStyleSheet('QCheckBox { padding-left: 5px; margin: 3px; }')
            a.setDefaultWidget(cb)
            self.addAction(a)
            cb.clicked.connect(partial(self.updateChannels, rm))

        self.addSeparator()
        modelAction = self.addAction('Model')
        modelAction.triggered.connect(self.setChannelsToModel)

        self.aboutToShow.connect(self.updateMenu)

    def updateMenu(self):
        for a in self.actions():
            try:
                cb = a.defaultWidget()
                cb.setChecked(False)
                if cb.text in list(itertools.chain.from_iterable([rms.split(',') for rms in self.rmsForActiveChannels])):
                    cb.setChecked(True)
            except Exception as e:
                print(e)

    def updateChannels(self, rmName, b):
        for c in [data.timeSeries(cn) for cn in self.activeChannels]:
            try:
                rms_for_ch = c.property('External standard').split(',') if c.property('External standard') else []
            except:
                rms_for_ch = []

            if 'Model' in rms_for_ch:
                rms_for_ch = []

            if b:
                rms_for_ch = list(set(rms_for_ch + [rmName]))
            else:
                rms_for_ch = list(filter(lambda rm: rm != rmName, rms_for_ch))

            c.setProperty('External standard', ','.join(rms_for_ch))

        self.rmsChanged.emit()

    def setChannelsToModel(self):
        for c in [data.timeSeries(cn) for cn in self.activeChannels]:
            c.setProperty('External standard', 'Model')

        self.rmsChanged.emit()


class InternalsDelegate(QStyledItemDelegate):

    def __init__(self, parent=None):
        QStyledItemDelegate.__init__(self, parent)

    def createEditor(self, parent, option, index):
        if index.column() == 2:
            return ChannelsComboBox(parent, index.data(Qt.UserRole))
        elif index.column() == 4:
            cb = QComboBox(parent)
            cb.addItems(['ppm', 'ppb', 'wtpc', 'wtpc_oxide'])
            return cb

        return QStyledItemDelegate.createEditor(self, parent, option, index)

    def setEditorData(self, editor, index):
        if index.column() == 2:
            editor.clear()
            editor.addItem(index.data(Qt.UserRole).property('Internal element'))
            editor.currentText = index.data(Qt.UserRole).property('Internal element')
        elif index.column() == 3:
            editor.text = str(index.data(Qt.UserRole).property('Internal value'))
        elif index.column() == 4:
            editor.currentText = index.data(Qt.UserRole).property('Internal units')

    def setModelData(self, editor, model, index):
        if index.column() in [2,4]:
            index.model().setData(index, editor.currentText) # These are QComboBox
        elif index.column() in [3]:
            index.model().setData(index, float(editor.text)) # These are QLineEdit

class InternalsModel(QAbstractTableModel):

    def __init__(self, parent):
        super().__init__(parent)
        self.refreshSelections()

    def refreshSelections(self):
        self.beginResetModel()
        self.selections = [[s for s in sg.selections()] for sg in data.selectionGroupList(data.ReferenceMaterial | data.Sample)]
        self.selections = list(itertools.chain.from_iterable(self.selections))
        self.endResetModel()

    def updateData(self, selections=None):
        # This is triggered when the channels are changed
        if selections is None:
            self.dataChanged.emit(self.index(0, 0), self.index(self.rowCount()-1, self.columnCount()-1))
            return

        for s in selections:
            try:
                is_names = s.property('Internal element').split(',')
            except:
                pass

            row = self.selections.index(s)

            if 'Criteria' in is_names:
                self.dataChanged.emit(self.index(row, 0), self.index(row, self.columnCount()-1))
                continue

            elements = [data.timeSeries(name).property('Element') for name in is_names if name]
            units = s.property('Internal units')

            if not units:
                units = 'ppm'
                s.setProperty('Internal units', units)

            try:
                if s.group().type == data.Sample:
                    sum = np.sum([Result(s.property(e), 0, 'wtpc', e).valueInUnits(units) for e in elements])
                elif s.group().type == data.ReferenceMaterial:
                    sum = np.sum([data.referenceMaterialData(s.group().name)[e].valueInUnits(units) for e in elements])
            except:
                sum = 0

            s.setProperty('Internal value', sum)
            self.dataChanged.emit(self.index(row, 0), self.index(row, self.columnCount()-1))

    def rowCount(self, parent=QModelIndex()):
        return len(self.selections)

    def columnCount(self, parent=QModelIndex()):
        return 6

    def headerData(self, section, orientation, role):
        if role != Qt.DisplayRole:
            return None

        if orientation == Qt.Horizontal:
            if section == 0:
                return 'Group'
            elif section == 1:
                return 'Selection'
            elif section == 2:
                return 'Element'
            elif section == 3:
                return 'Value'
            elif section == 4:
                return 'Units'
            elif section == 5:
                return 'Affinity'

        return None

    def data(self, index, role):
        if not index.isValid():
            return None

        if role == Qt.UserRole:
            return self.selections[index.row()]

        if role == Qt.ToolTipRole:
            return self.data(index, Qt.DisplayRole)

        if role != Qt.DisplayRole:
            return None

        if index.column() == 0:
            return self.selections[index.row()].group().name
        elif index.column() == 1:
            return self.selections[index.row()].name
        elif index.column() == 2:
            return self.selections[index.row()].property('Internal element')
        elif index.column() == 3:
            if self.selections[index.row()].property('Internal element') == 'Criteria':
                return '-'
            return self.selections[index.row()].property('Internal value')
        elif index.column() == 4:
            if self.selections[index.row()].property('Internal element') == 'Criteria':
                return '-'
            return self.selections[index.row()].property('Internal units')
        elif index.column() == 5:
            s = self.selections[index.row()]
            afe = s.property('Affinity elements')
            ext = s.property('External affinity')
            if ext:
                return f'{ext} (from {afe})'

            return f'{afe}'

        return None

    def setData(self, index, value, role = Qt.EditRole):
        if role == Qt.EditRole and index.column() == 2:
            index.data(Qt.UserRole).setProperty('Internal element', value)
            self.dataChanged.emit(index, index)
        elif role == Qt.EditRole and index.column() == 3:
            index.data(Qt.UserRole).setProperty('Internal value', value)
            self.dataChanged.emit(index, index)
        elif role == Qt.EditRole and index.column() == 4:
            index.data(Qt.UserRole).setProperty('Internal units', value)
            self.dataChanged.emit(index, index)

    def flags(self, index):
        if index.column() > 1:
            return QAbstractTableModel.flags(self, index) | Qt.ItemIsEditable
        return QAbstractTableModel.flags(self, index)


class ChannelsComboBox(QComboBox):

    def __init__(self, parent, selection):
        QComboBox.__init__(self, parent)
        self.selection = selection

    def updateText(self):
        self.clear()
        self.addItem(self.selection.property('Internal element'))
        self.currentText = self.selection.property('Internal element')

    def showPopup(self):
        menu = ChannelsMenu(self)
        menu.setSelections([self.selection])
        menu.setFixedWidth(self.width)
        menu.channelsChanged.connect(self.updateText)
        p = self.pos
        p += QPoint(0, self.height)
        menu.exec_(self.parent().mapToGlobal(p))

class ChannelsMenu(QMenu):

    channelsChanged = Signal(list)

    def __init__(self, parent, propName='Internal element'):
        super().__init__(parent)
        self.propName = propName

        for channel in [c for c in data.timeSeriesList(data.Input) if 'TotalBeam' not in c.name]:
            a = QWidgetAction(self)
            cb = QCheckBox(channel.name, self)
            cb.setStyleSheet('QCheckBox { padding-left: 5px; margin: 3px; }')
            a.setDefaultWidget(cb)
            self.addAction(a)
            cb.clicked.connect(partial(self.updateSelections, channel.name))

        self.addSeparator()
        selectAll = self.addAction('Select all')
        selectAll.triggered.connect(self.selectAll)
        selectNone = self.addAction('Select none')
        selectNone.triggered.connect(self.selectNone)
        self.addSeparator()
        crit = self.addAction('Use criteria')
        crit.triggered.connect(self.setSelectionsForCriteria)

        self.aboutToShow.connect(self.updateMenu)
        self.selections = []
        self.channels = []

    def setSelections(self, sels):
        self.selections = sels

    def updateMenu(self):
        try:
            s = self.selections[0]
            chs = s.property(self.propName).split(',')
            self.setChannels(chs)
        except Exception as e:
            print(e)

    def selectAll(self):
        self.channels = [c for c in data.timeSeriesNames(data.Input) if 'TotalBeam' not in c]
        for s in self.selections:
            s.setProperty(self.propName, ','.join(self.channels))
        self.channelsChanged.emit(self.selections)

    def selectNone(self):
        for s in self.selections:
            s.setProperty(self.propName, '')
        self.channels = []

        # If this is resetting the external affinities, also reset the Ext Affinity
        if self.propName == 'Affinity elements':
            for s in self.selections:
                s.setProperty('External affinity', '')

        self.channelsChanged.emit(self.selections)

    def setSelectionsForCriteria(self):
        for s in self.selections:
            s.setProperty(self.propName, 'Criteria')

        self.channelsChanged.emit(self.selections)

    def setChannels(self, channels):
        self.channels = channels

        for a in self.actions():
            try:
                a.defaultWidget().setChecked(a.defaultWidget().text in self.channels)
            except Exception as e:
                continue

    def updateSelections(self, channelName, b):

        self.channels = []
        for a in self.actions():
            try:
                if a.defaultWidget().isChecked():
                    self.channels.append(a.defaultWidget().text)
            except Exception as e:
                continue

        for s in self.selections:
            if s.property(self.propName):
                sie = s.property(self.propName).split(',')
                sie[:] = (v for v in sie if v != 'Criteria' and v not in data.referenceMaterialNames())
            else:
                sie = []
            if channelName not in sie and b: sie.append(channelName)
            if channelName in sie and not b: sie.remove(channelName)
            s.setProperty(self.propName, ','.join(sie))

        self.channelsChanged.emit(self.selections)


class CriteriaModel(QAbstractTableModel):

    def __init__(self, criteria, parent):
        QAbstractTableModel.__init__(self, parent)
        self.criteria = criteria

    def setCriteria(self, criteria, b=None): # note: the b is just to make connections from action.triggered happy
        if not criteria:
            criteria = []

        self.beginResetModel()
        self.criteria = list(criteria)
        self.endResetModel()

    def rowCount(self, index=QModelIndex()):
        return len(self.criteria)

    def addRow(self):
        self.beginInsertRows(QModelIndex(), self.rowCount(), self.rowCount())
        self.criteria.append({
            'name': 'Unnamed',
            'criteria': '',
            'analytes': '',
            'value': 100.,
            'oxides': True,
            'oxide_forms': []
        })
        self.endInsertRows()

    def removeRow(self, row):
        self.beginRemoveRows(QModelIndex(), row, row)
        del self.criteria[row]
        self.endRemoveRows()

    def moveRow(self, row, dir):
        if row == 0 and dir < 0:
            return
        elif row == self.rowCount() - 1 and dir > 0:
            return

        dirmod = 0 if dir < 0 else 1

        if self.beginMoveRows(QModelIndex(), row, row, QModelIndex(), row+dir+dirmod):
            self.criteria[row], self.criteria[row+dir] = self.criteria[row+dir], self.criteria[row]
            self.endMoveRows()

    def columnCount(self, index):
        return 6

    def headerData(self, section, orientation, role):
        if role != Qt.DisplayRole:
            return None

        if orientation == Qt.Horizontal:
            if section == 0:
                return 'Name'
            elif section == 1:
                return 'Criteria'
            elif section == 2:
                return 'Analytes'
            elif section == 3:
                return 'Value'
            elif section == 4:
                return 'Oxides'
            elif section == 5:
                return 'Oxide forms'

        return None

    def data(self, index, role):
        if not index.isValid():
            return None

        if role == Qt.DisplayRole or role == Qt.EditRole:
            if index.column() == 0:
                return self.criteria[index.row()]['name']
            elif index.column() == 1:
                return self.criteria[index.row()]['criteria']
            elif index.column() == 2:
                return self.criteria[index.row()]['analytes']
            elif index.column() == 3:
                return float(self.criteria[index.row()]['value'])
            elif index.column() == 4:
                return self.criteria[index.row()]['oxides']
            elif index.column() == 5:
                return ','.join(self.criteria[index.row()]['oxide_forms'])
        elif role == Qt.CheckStateRole and index.column() == 4:
            return Qt.Checked if self.criteria[index.row()]['oxides'] else Qt.Unchecked


        return None

    def setData(self, index, value, role = Qt.EditRole):
        if role == Qt.CheckStateRole and index.column() == 4:
            self.criteria[index.row()]['oxides'] = value == Qt.Checked

        if index.column() == 0:
            self.criteria[index.row()]['name'] = value
        elif index.column() == 1:
            self.criteria[index.row()]['criteria'] = value
        elif index.column() == 2:
            self.criteria[index.row()]['analytes'] = value
        elif index.column() == 3:
            self.criteria[index.row()]['value'] = float(value)
        elif index.column() == 5:
            self.criteria[index.row()]['oxide_forms'] = str(value).split(',')

        self.dataChanged.emit(index, index)

    def flags(self, index):
        if index.column() == 4:
            return QAbstractTableModel.flags(self, index) | Qt.ItemIsUserCheckable

        if index.column() != 4:
            return QAbstractTableModel.flags(self, index) | Qt.ItemIsEditable

        return QAbstractTableModel.flags(self, index)

class CriteriaDialog(QDialog):

    def __init__(self, parent=None):
        QDialog.__init__(self, parent)
        self.setWindowTitle('Internal standard criteria')
        self.table = QTableView(self)
        self.table.setModel(CriteriaModel([],self.table))
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.contextMenuPolicy = Qt.CustomContextMenu

        def makeButton(text, icon):
            b = QToolButton(self)
            b.text = text
            b.setIcon(CUI().icon(icon))
            b.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
            return b

        self.addButton = makeButton('Add', 'plus')
        self.removeButton = makeButton('Remove', 'minus')
        self.upButton = makeButton('Up', 'arrowup')
        self.downButton = makeButton('Down', 'arrowdown')
        self.presetsButton = makeButton('Presets', 'book')
        self.bb = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel, self)

        topLayout = QHBoxLayout()
        topLayout.addWidget(self.addButton)
        topLayout.addWidget(self.removeButton)
        topLayout.addWidget(self.upButton)
        topLayout.addWidget(self.downButton)
        topLayout.addStretch(1)
        topLayout.addWidget(self.presetsButton)

        self.setLayout(QVBoxLayout())
        self.layout().addLayout(topLayout)
        self.layout().addWidget(self.table)
        self.layout().addWidget(self.bb)

        self.bb.accepted.connect(self.accept)
        self.bb.rejected.connect(self.reject)
        self.addButton.clicked.connect(self.table.model().addRow)
        self.removeButton.clicked.connect(lambda: self.table.model().removeRow(self.table.selectionModel().selectedRows()[0].row()))
        self.upButton.clicked.connect(lambda: self.table.model().moveRow(self.table.selectionModel().selectedRows()[0].row(), -1))
        self.downButton.clicked.connect(lambda: self.table.model().moveRow(self.table.selectionModel().selectedRows()[0].row(), 1))
        self.table.customContextMenuRequested.connect(self.showAnalytesMenu)

        self.setupPresetsMenu()

        self.resize(700, 500)

    def setupPresetsMenu(self):
        self.presetsMenu = QMenu(self)
        self.presetsButton.setMenu(self.presetsMenu)
        self.presetsButton.setPopupMode(QToolButton.InstantPopup)

        settings = QSettings()
        presets = dict(settings.value('3DTE/CriteriaPresets', {}))

        for k in presets:
            a = self.presetsMenu.addAction(k)
            a.setProperty('criteria', presets[k])
            a.triggered.connect(partial(self.table.model().setCriteria, presets[k]))

        self.presetsMenu.addSeparator()
        saveAction = self.presetsMenu.addAction('Save')
        saveAction.triggered.connect(self.savePreset)
        removeAction = self.presetsMenu.addAction('Remove')
        removeAction.triggered.connect(self.removePreset)

    def savePreset(self):
        name = QInputDialog.getText(self, 'Save criteria preset', 'Preset name')

        if not name:
            return

        presets = QSettings().value('3DTE/CriteriaPresets', {})
        presets[name] = self.table.model().criteria
        QSettings().setValue('3DTE/CriteriaPresets', presets)
        self.setupPresetsMenu()

    def removePreset(self):
        names = list(QSettings().value('3DTE/CriteriaPresets', {}).keys())
        name = QInputDialog.getItem(self, 'Remove criteria preset', 'Preset name', names)

        if not name:
            return

        presets = QSettings().value('3DTE/CriteriaPresets', {})
        del presets[name]
        QSettings().setValue('3DTE/CriteriaPresets', presets)
        self.setupPresetsMenu()

    def showAnalytesMenu(self, point):
        index = self.table.indexAt(point)
        if index.column() != 2:
            return

        menu = ChannelsMenu(self)
        menu.removeAction(menu.actions()[-1])
        menu.setChannels(index.data().split(','))
        menu.channelsChanged.connect(lambda: index.model().setData(index, ','.join(menu.channels)))
        menu.exec_(self.mapToGlobal(point))

    def criteria(self):
        return self.table.model().criteria


class BlockAssignmentsModel(QAbstractTableModel):

    def __init__(self, selections, parent):
        super().__init__(parent)
        self.selections = selections
        self.blocks = []
        self.selToBlockMap = {}

    def setBlocks(self, blocks):
        self.beginResetModel()
        self.blocks = blocks
        self.selToBlockMap = {}

        for block in self.blocks:
            for sel in block.selections:
                self.selToBlockMap[sel.property('UUID')] = block.label

        self.endResetModel()

    def reset(self):
        self.beginResetModel()
        self.setBlocks(self.blocks)
        self.endResetModel()

    def rowCount(self, index=QModelIndex()):
        return len(self.selections)

    def columnCount(self, index=QModelIndex()):
        return 4

    def headerData(self, section, orientation, role):
        if role != Qt.DisplayRole:
            return None

        if orientation == Qt.Horizontal:
            if section == 0:
                return 'Group'
            elif section == 1:
                return 'Selection'
            elif section == 2:
                return 'Time'
            elif section == 3:
                return 'Block'

        return None

    def data(self, index, role):
        if not index.isValid():
            return None

        s = self.selections[index.row()]

        if role == Qt.UserRole:
            return s

        try:
            blockIndex = self.selToBlockMap[s.property('UUID')]
        except:
            blockIndex = None

        if role == Qt.BackgroundColorRole:
            if not blockIndex:
                return None

            return colors[blockIndex]

        if role == Qt.DecorationRole and index.column() == 3:
            if s.property('Block') is not None:
                return CUI().icon('bullseye')

        if role != Qt.DisplayRole:
            return

        if index.column() == 0:
            return s.group().name
        elif index.column() == 1:
            return s.name
        elif index.column() == 2:
            return s.startTime.toString('yyyy-MM-dd hh:mm:ss.zzz')
        elif index.column() == 3:
            return blockIndex

        return None

class BlockPlot(QWidget):

    def __init__(self, parent):
        if parent is None:
            raise Exception("The 3D Traces Block Plot was initialised without a parent")

        super().__init__(parent)
        self.settingsWidget = parent
        self.setLayout(QVBoxLayout())

        self.buttonLayout = QHBoxLayout()
        self.blockDownButton = QToolButton(self)
        self.blockDownButton.setIcon(CUI().icon('arrowleft'))
        self.blockLabel = QLabel("")
        self.blockUpButton = QToolButton(self)
        self.blockUpButton.setIcon(CUI().icon('arrowright'))
        self.configButton = QToolButton(self)
        self.configButton.setIcon(CUI().icon('wrench'))
        self.configButton.setFixedSize(25, 25)
        self.configButton.clicked.connect(self.configBlocks)

        self.buttonLayout.addWidget(self.blockDownButton)
        self.buttonLayout.addWidget(self.blockLabel)
        self.buttonLayout.addWidget(self.configButton)
        self.buttonLayout.addWidget(self.blockUpButton)
        self.buttonLayout.insertStretch(1)
        self.buttonLayout.insertStretch(4)
        self.layout().addLayout(self.buttonLayout)
        self.blockDownButton.clicked.connect(lambda: self.setBlock(self.bn - 1))
        self.blockUpButton.clicked.connect(lambda: self.setBlock(self.bn + 1))
        self.blockDownButton.setDisabled(True)
        #self.setAutoFillBackground(True)
        self.plot = Plot(self)
        self.plot.setBackground(CUI().tabBackgroundColor())
        self.plot.setToolsVisible(False)
        self.plot.left().label = 'Intensity'
        self.plot.bottom().label = 'Concentration'
        self.layout().addWidget(self.plot)
        
        self.logButton = OverlayButton(self.plot, 'BottomLeft', 0, 0)
        self.logButton.setCheckable(True)
        self.logButton.setText('10ⁿ')
        self.logButton.setFixedSize(25, 25)
        self.logButton.clicked.connect(self.toggleLog)
        self.plot.installEventFilter(self.logButton)
        self.isLog = False

        self.showLegend = False
        menu = self.plot.contextMenu()
        menu.addSeparator()
        la = menu.addAction('Legend')
        la.setCheckable(True)
        la.setChecked(self.showLegend)
        la.triggered.connect(self.toggleLegend)

        self.bn = 0

    def toggleLegend(self, b):
        self.showLegend = b
        self.updatePlot()

    def configBlocks(self):
        d = QDialog(self)
        d.setWindowTitle('Block assignments')
        d.setLayout(QVBoxLayout())
        d.resize(600, 600)

        topLayout = QHBoxLayout()
        setButton = QPushButton(d)
        setButton.setText('Set selected')
        topLayout.addWidget(setButton)
        clearButton = QPushButton(d)
        clearButton.setText('Clear selected')
        topLayout.addWidget(clearButton)
        saveAssignmentsButton = QPushButton(d)
        saveAssignmentsButton.setText('Save as assignments')
        topLayout.addWidget(saveAssignmentsButton)
        topLayout.addStretch()
        topLayout.addWidget(QLabel('Method', d))
        methodComboBox = QComboBox(d)
        methodComboBox.addItems(['Assigned', 'Simple', 'Clustering', 'Auto Clustering'])
        if drs.setting('BlockFindingMethod'):
            methodComboBox.setCurrentText(drs.setting('BlockFindingMethod'))
        else:
            methodComboBox.setCurrentText('Simple')
        nClustersSpinBox = QSpinBox(d)
        nClustersSpinBox.setMinimum(1)
        nClustersSpinBox.setMaximum(1000)
        try:
            nClusters = int(drs.setting('NClusters')) if int(drs.setting('NClusters')) > 0 else 5
            nClustersSpinBox.setValue(nClusters)
        except Exception as e:
            print(e)

        topLayout.addWidget(methodComboBox)
        topLayout.addWidget(nClustersSpinBox)
        nClustersSpinBox.setVisible(drs.setting('BlockFindingMethod') and drs.setting('BlockFindingMethod') == 'Clustering')
        methodComboBox.activated.connect(lambda: nClustersSpinBox.setVisible(methodComboBox.currentText == 'Clustering'))
        d.layout().addLayout(topLayout)

        table = QTableView(d)
        d.layout().addWidget(table)

        plot = Plot(d)
        g = plot.addGraph()
        tb = data.timeSeries('TotalBeam')
        g.setData(tb.time(), tb.data())
        plot.bottom().label = 'Time (s)'
        plot.left().label = 'TotalBeam'
        plot.setFixedHeight(200)
        plot.setToolsVisible(False)
        plot.bottom().setDateTime(True)
        plot.left().setLogarithmic(True)
        d.layout().addWidget(plot)

        bb = QDialogButtonBox(QDialogButtonBox.Close, d)
        bb.accepted.connect(d.accept)
        bb.rejected.connect(d.reject)
        d.layout().addWidget(bb)

        table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        table.setSortingEnabled(True)
        table.setSelectionMode(table.ExtendedSelection)
        table.setSelectionBehavior(table.SelectRows)

        externalsInUse = set(list(itertools.chain(*[c.property('External standard').split(',') for c in data.timeSeriesList(data.Input) if 'TotalBeam' not in c.name])))
        try:
            externalsInUse.remove('Model')
        except KeyError:
            pass

        groups = [data.selectionGroup(ext) for ext in externalsInUse if ext]
        selections = list(itertools.chain(*[sg.selections() for sg in groups]))
        selections.sort(key=lambda s: s.midTimeInSec)
        selMidTimes = [s.midTimeInSec for s in selections]

        blockModel = BlockAssignmentsModel(selections, d)
        table.setModel(blockModel)

        def update():
            print('update...')
            blocks = findBlocks(drs.setting('BlockFindingMethod'))
            blockModel.setBlocks(blocks)

            # Update plot
            plot.clearItems()
            blockSelections = []

            for block in blocks:
                minTime = np.min([s.startTime.toMSecsSinceEpoch()/1000.0 for s in block.selections])
                maxTime = np.max([s.endTime.toMSecsSinceEpoch()/1000.0 for s in block.selections])
                rect = plot.addRect( (minTime, 1), (maxTime, 0), 'ptPlotCoords', 'ptPlotCoords', plot.bottom(), plot.right())
                c = colors[block.label]
                c.setAlpha(90)
                rect.brush = QBrush(c)

            plot.right().range = QCPRange(0, 1)
            plot.rescaleAxes()
            plot.replot()

        update()

        def setBlocks():
            block = QInputDialog.getInt(d, 'Block number', 'Block number:', 1, 0, 100, 1)

            if block is None:
                return

            sels = [index.data(Qt.UserRole) for index in table.selectionModel().selectedRows()]
            for sel in sels:
                print(f'Setting block number for {sel.name} to {block}')
                sel.setProperty('Block', block)

            update()

        setButton.clicked.connect(setBlocks)

        def clearAssignments():
            sels = [index.data(Qt.UserRole) for index in table.selectionModel().selectedRows()]
            for sel in sels:
                sel.setProperty('Block', None)

            update()

        clearButton.clicked.connect(clearAssignments)

        def saveAssignments():
            for sel in blockModel.selections:
                sel.setProperty('Block', blockModel.selToBlockMap[sel.property('UUID')])

        saveAssignmentsButton.clicked.connect(saveAssignments)

        nClustersSpinBox.valueChanged.connect(lambda v: drs.setSetting('NClusters', v))
        nClustersSpinBox.valueChanged.connect(update)
        methodComboBox.activated.connect(lambda: drs.setSetting('BlockFindingMethod', methodComboBox.currentText))
        methodComboBox.activated.connect(update)

        d.exec_()
        self.settingsWidget.calibration.updateBlocks()
        self.settingsWidget.processExtSelection()

    def setBlock(self, bn):
        if bn < 0 or bn >= len(self.settingsWidget.calibration.blocks):
            print(f'Tried to set block to an invalid number... {bn}')
            return

        self.blockDownButton.setEnabled(bn > 0 and len(self.settingsWidget.calibration.blocks) > 0)
        self.blockUpButton.setEnabled(bn < (len(self.settingsWidget.calibration.blocks) - 1) and len(self.settingsWidget.calibration.blocks) > 0)

        self.bn = bn
        self.updatePlot()

    def toggleLog(self, b):        
        self.plot.left().setLogarithmic(b)
        self.plot.bottom().setLogarithmic(b)
        self.isLog = b
        if self.isLog:
            self.plot.bottom().range = QCPRange(0.001, self.x_max*5)
            self.plot.left().range = QCPRange(1, self.y_max*5)
        else:
            self.plot.bottom().range = QCPRange(0, self.x_max*1.1)
            self.plot.left().range = QCPRange(0, self.y_max*1.1)
        self.plot.replot()        

    def updatePlot(self):
        self.plot.clearGraphs()
        self.plot.clearItems()

        channel = self.settingsWidget.selectedChannelNames[0]
        if len(self.settingsWidget.calibration.blocks) < 1:
            print('No blocks to update BlockPlot with...')
            return

        extStd = data.timeSeries(channel).property('External standard')
        if type(extStd) != str or (extStd.split(',')[0] not in data.selectionGroupNames() and extStd != 'Model'):
            print('Invalid external standard so not going to update plot.')
            self.plot.rescaleAxes()
            self.plot.replot()
            return

        # First, make sure we have RM values for this channel:
        this_df = self.settingsWidget.calibration.blocks[0].dataFrame()[f'{channel}_RMppm']
        if not this_df.notna().values.any() and extStd != 'Model':
            print(f'No valid values for reference materials for {channel}')
            self.plot.rescaleAxes()
            self.plot.replot()
            return

        # Make sure that the axes have the same extents for all blocks
        # Makes it easier to see differences
        self.y_max = np.nanmax([np.nanmax(b.dataFrame()[f'{channel}']) for b in self.settingsWidget.calibration.blocks])
        self.y_min = np.nanmin([np.nanmin(b.dataFrame()[f'{channel}']) for b in self.settingsWidget.calibration.blocks])

        if np.all(np.isnan(np.array([b.slope(channel) for b in self.settingsWidget.calibration.blocks]))):
            print(f'No valid values for reference materials for {channel}')
            self.plot.rescaleAxes()
            self.plot.replot()
            return

        if extStd == 'Model':
            self.x_max = np.nanmax([self.y_max/b.slope(channel) for b in self.settingsWidget.calibration.blocks])
        else:
            self.x_max = np.nanmax([np.nanmax(b.dataFrame()[f'{channel}_RMppm']) for b in self.settingsWidget.calibration.blocks])


        if np.isnan(self.x_max):
            print(f'No valid values for reference materials for {channel}')
            self.plot.rescaleAxes()
            self.plot.replot()
            return

        block = self.settingsWidget.calibration.blocks[self.bn]
        df = block.dataFrameForChannel(channel)
        slope = block.slope(channel)
        intercept = block.intercept(channel)

        ss_res = np.sum( (df[channel] - (slope*df[f'{channel}_RMppm'] + intercept))**2 )
        ss_tot = np.sum( (df[channel] - np.mean(df[channel]))**2 )
        #Note, for modeled fits, the r_sq value will always be 1. This is a little misleading,
        # so will manually change this below, when the annotation is added to the plot
        r_sq = 1 - (ss_res / ss_tot) if len(df) > 1 else 1

        try:
            x_vals = np.logspace(-3, ceil(log(self.x_max)), 200)
        except:
            x_vals = np.logspace(-3, 3, 200)
            self.x_max = 1000.

        y_vals = slope * x_vals + intercept

        #Add error envelope
#        res = block.fit(channel)['sm_res']
#        try:
#            if data.timeSeries(channel).property('FitThroughZero'):
#                pred = res.get_prediction(x_vals)
#            else:
#                pred = res.get_prediction(sm.add_constant(x_vals))
#            frame = pred.summary_frame(alpha=0.05)

#            print(frame.obs_ci_upper)
#            pg = self.plot.addGraph()
#            pg.setData(x_vals, frame.obs_ci_upper)

#            mg = self.plot.addGraph()
#            mg.setData(x_vals, frame.obs_ci_lower)
#        except Exception as e:
#            print(f'Could not use wls_prediction_std... {e}')

        # Use actual graph with lots of data points so it looks right in log scale, straight line DOES NOT
        g_cal = self.plot.addGraph()
        g_cal.setData(x_vals, y_vals)
        g_cal.setColor(Qt.red)
        g_cal.removeFromLegend()

        grad = QCPColorGradient('gpJet')
        symbols = [
            'ssCross', 'ssPlus', 'ssCircle', 'ssSquare', 'ssDiamond', 'ssStar', 'ssTriangle',
            'ssTriangleInverted', 'ssCrossSquare', 'ssPlusSquare', 'ssCrossCircle', 'ssPlusCircle'
            ]

        for gi, groupName in enumerate(df['group'].unique()):
            gdf = df[df['group'] == groupName]
            g_data = self.plot.addGraph()
            g_data.setData(gdf[f'{channel}_RMppm'], gdf[channel])
            color = grad.color(gi, 0, len(df['group'].unique()))
            g_data.setScatterStyle(symbols[gi%len(symbols)], 8, color, color)
            g_data.setLineStyle('lsNone')
            g_data.setColor(color)
            g_data.setName(groupName)

            eby = QCPErrorBars(self.plot.bottom(), self.plot.left())
            eby.setData(gdf[f'{channel}_Uncert'])
            eby.setDataPlottable(g_data)
            eby.errorType = QCPErrorBars.etValueError
            eby.removeFromLegend()

            ebx = QCPErrorBars(self.plot.bottom(), self.plot.left())
            ebx.setData(gdf[f'{channel}_RMppm_Uncert'])
            ebx.setDataPlottable(g_data)
            ebx.errorType = QCPErrorBars.etKeyError
            ebx.removeFromLegend()

        # Add slope and intercept annotation here:
        ann = self.plot.annotate('', 0.01, 0.01, 'ptAxisRectRatio', Qt.AlignLeft | Qt.AlignTop)
        if extStd == 'Model':
            ann.text = '''
                    <p style="color:black;">
                    <b>Slope</b>: %s<br>
                    <b>Intercept</b>: %s<br>
                    <b><i>R²</i></b> : N/A</p>'''%(formatResult(slope, block.fit(channel)['slope_uncert'])[0], formatResult(intercept, block.fit(channel)['intercept_uncert'])[0])

            self.plot.annotate('*Modeled', 0.9, 0.9, 'ptAxisRectRatio', Qt.AlignCenter, Qt.AlignCenter, False)

        else:
            ann.text = '''
                    <p style="color:black;">
                    <b>Slope</b>: %s<br>
                    <b>Intercept</b>: %s<br>
                    <b><i>R²</i></b> : %.3f</p>'''%(formatResult(slope, block.fit(channel)['slope_uncert'])[0], formatResult(intercept, block.fit(channel)['intercept_uncert'])[0], r_sq)

        if self.isLog:
            self.plot.bottom().range = QCPRange(0.001, self.x_max*5)
            self.plot.left().range = QCPRange(1, self.y_max*5)
        else:    
            self.plot.bottom().range = QCPRange(0, self.x_max*1.1)
            self.plot.left().range = QCPRange(0, self.y_max * 1.1)

        self.plot.setLegendVisible(self.showLegend)
        self.plot.setLegendAlignment(Qt.AlignBottom | Qt.AlignRight)
        self.plot.replot()

        self.blockLabel.setText(f"Block {self.bn + 1}/{len(self.settingsWidget.calibration.blocks)}")


class FitsPlot(Plot):

    def __init__(self, parent):
        super().__init__(parent)
        self.settingsWidget = parent
        self.setToolsVisible(False)
        self.setBackground(CUI().tabBackgroundColor())
        self.left().label = 'Intensity'
        self.bottom().label = 'Concentration'
        self.grad = QCPColorGradient('gpViridis')

        self.setLegendAlignment(Qt.AlignTop | Qt.AlignLeft)
        self.showLegend = False
        self.showColorScale = False
        self.blockCount = 0

        menu = self.contextMenu()
        menu.addSeparator()
        la = menu.addAction('Legend')
        la.setCheckable(True)
        la.setChecked(self.showLegend)
        la.triggered.connect(self.toggleLegend)

        ca = menu.addAction('Color scale')
        ca.setCheckable(True)
        ca.setChecked(self.showColorScale)
        ca.triggered.connect(self.toggleColorScale)

    def toggleLegend(self, b):
        self.showLegend = b
        self.setLegendVisible(b)
        self.replot()

    def toggleColorScale(self, b):
        self.showColorScale = b
        if not b:
            self.removeColorScale()
        else:
            self.addColorScale('Block number', self.grad, 1, self.blockCount)

        self.replot()


    def updatePlot(self):
        self.clearGraphs()
        self.clearItems()

        block_times = [block.midTime() for block in self.settingsWidget.calibration.blocks]
        channel = self.settingsWidget.selectedChannelNames[0]
        self.blockCount = len(self.settingsWidget.calibration.blocks)
        extStd = data.timeSeries(channel).property('External standard')
        if type(extStd) != str or (extStd.split(',')[0] not in data.selectionGroupNames() and extStd != 'Model'):
            self.replot()
            return

        for i, block in enumerate(self.settingsWidget.calibration.blocks):
            df = self.settingsWidget.calibration.blocks[i].dataFrameForChannel(channel)

            if extStd == 'Model':
                x_max = 100.
                self.annotate('*Modeled', 0.9, 0.9, 'ptAxisRectRatio', Qt.AlignCenter, Qt.AlignCenter, False)

            else:
                if not df[f'{channel}_RMppm'].notna().values.any():
                    self.replot()
                    return

                x_max = df[f'{channel}_RMppm'].dropna().max()
                x_max += x_max * 0.1

            slope = block.slope(channel)
            intercept = block.intercept(channel)

            x_vals = np.linspace(0, x_max)
            y_vals = slope * x_vals + intercept
            g = self.addGraph()
            g.setName(f'Block {block.label}')
            g.setColor(self.grad.color(i+1, 1, self.blockCount))
            g.setData(x_vals, y_vals)

        self.rescaleAxes()
        self.replot()


class FitParamsPlot(Plot):

    def __init__(self, parent):
        super().__init__(parent)
        self.settingsWidget = parent
        self.setToolsVisible(False)
        self.setBackground(CUI().tabBackgroundColor())
        self.left().label = 'Slope'
        self.right().label = 'Intercept'
        self.bottom().label = 'Time'
        self.bottom().setDateTime(True)
        self.top().label = 'Block'
        self.right().visible = True
        self.top().visible = True

        self.right().labelColor = QColor(Qt.red)
        self.left().labelColor = QColor(Qt.blue)

    def updatePlot(self):
        self.clearGraphs()
        self.clearItems()

        channel = self.settingsWidget.selectedChannelNames[0]
        extStd = data.timeSeries(channel).property('External standard')
        if type(extStd) != str or (extStd.split(',')[0] not in data.selectionGroupNames() and extStd != 'Model'):
            self.replot()
            return

        block_times = [block.midTime() for block in self.settingsWidget.calibration.blocks]

        if not self.settingsWidget.calibration.blocks[0].dataFrame()[f'{channel}_RMppm'].notna().values.any() and extStd != 'Model':
            self.replot()
            return

        slopes = [block.slope(channel) for block in self.settingsWidget.calibration.blocks]
        slopes_err = [block.slopeUncert(channel) for block in self.settingsWidget.calibration.blocks]
        inters = [block.intercept(channel) for block in self.settingsWidget.calibration.blocks]
        inters_err = [block.interceptUncert(channel) for block in self.settingsWidget.calibration.blocks]

        slopes_graph = self.addGraph(self.bottom(), self.left())
        slopes_eb = QCPErrorBars(self.bottom(), self.left())
        slopes_eb.setDataPlottable(slopes_graph)
        slopes_graph.setScatterStyle('ssCircle')

        inters_graph = self.addGraph(self.bottom(), self.right())
        inters_graph.setScatterStyle('ssDiamond', 6, Qt.red, Qt.red)
        inters_graph.setColor(Qt.red)
        inters_eb = QCPErrorBars(self.bottom(), self.right())
        inters_eb.setDataPlottable(inters_graph)

        slopes_graph.setData(block_times, slopes)
        slopes_eb.setData(slopes_err)

        inters_graph.setData(block_times, inters)
        inters_eb.setData(inters_err)

        self.rescaleAxes()
        self.top().range = QCPRange(1, len(block_times))
        self.bottom().scaleRange(1.1)
        self.left().scaleRange(1.1)
        self.top().scaleRange(1.1)

        if extStd == 'Model':
            self.left().label = 'Slope (Modeled)'
            self.right().label = 'Intercept (Modeled)'
        else:
            self.left().label = 'Slope'
            self.right().label = 'Intercept'

        self.replot()


class FractionationPlot(Plot):

    def __init__(self, parent):
        super().__init__(parent)
        self.settingsWidget = parent
        self.setToolsVisible(False)
        self.setBackground(CUI().tabBackgroundColor())
        self.bottom().label = 'Beam seconds'
        self.setLegendVisible(True)

        self.msg = OverlayMessage(
                    self, 'Warning',
                    'No internal standards set.',
                    0, 0.8, 40, 10
                    )
        self.msg.setCloseButtonVisible(False)
        self.msg.hide()

    def updatePlot(self):
        self.clearGraphs()  
        self.clearItems()
        try:
            channel = self.settingsWidget.selectedChannelNames[0]
        except:
            return

        extStd = data.timeSeries(channel).property('External standard')

        if extStd == 'Model':
            text  = self.annotate('Fractionation correction not available\nfor channels with modelled sensitivity',
                          0.5, 0.5, 'ptAxisRectRatio', Qt.AlignCenter, Qt.AlignCenter, False)
            text.padding = QMargins(10,10,10,10)
            self.replot()
            return

        if type(extStd) != str or extStd.split(',')[0] not in data.selectionGroupNames():
            self.replot()
            return

        fdf = self.settingsWidget.calibration.fractionation(channel)

        if len(fdf) == 0:
            self.msg.show()
            self.replot()
            return

        self.msg.hide()
        grad = QCPColorGradient('gpViridis')
        externalsInUse = list(set(list(itertools.chain(*[c.property('External standard').split(',') for c in data.timeSeriesList(data.Input) if 'TotalBeam' not in c.name]))))

        for i, intStd in enumerate(fdf['IS'].unique()):
            for ei, ext in enumerate(externalsInUse):
                t, r, rsd, fit = self.settingsWidget.calibration.fitFractionation(channel, intStd, group=ext)
                if t is None:
                    continue

                color = grad.color(i+ei+1, 1, len(fdf['IS'].unique())+len(externalsInUse))
                g = self.addGraph()
                g.setData(t, r)
                g.setLineStyle('lsNone')
                g.setScatterStyle('ssDisc', 6, color, color)
                g.setColor(color)
                name = (intStd[:20]+'...') if len(intStd) > 20 else intStd
                g.setName(f'{name} - {ext}')

                self.eb = QCPErrorBars(self.bottom(), self.left())
                self.eb.setDataPlottable(g)
                self.eb.setData(rsd)

                self.eb.removeFromLegend()

                if fit:
                    sg = self.addGraph()
                    sx = np.linspace(t.min(), t.max(), 100)
                    sg.setData(sx, fit(sx))
                    sg.setColor(color)
                    sg.removeFromLegend()

#        try:
#            sel = self.settingsWidget.sels[0]
#            num = data.timeSeries(channel+'_CPS').dataForSelection(sel)
#            den = data.timeSeries(sel.property('Internal element')+'_CPS').dataForSelection(sel)
#            r = savgol_filter(num, 11, 3)/savgol_filter(den, 11, 3)
#            g2 = self.addGraph(self.bottom(), self.right())
#            t = data.timeSeries(channel+'_CPS').timeForSelection(sel)
#            t -= t[0]
#            g2.setData(t, r)
#            print(r)

#        except:
#            print('No selection for fractionation plot')


        self.left().label = f'Normalized {data.timeSeries(channel).property("Element")}/Internal standard'
        self.rescaleAxes()
        self.replot()


class JacksonPlot(Plot):

    def __init__(self, parent):
        super().__init__(parent)
        self.settingsWidget = parent
        self.setToolsVisible(False)
        self.setBackground(CUI().tabBackgroundColor())
        self.bottom().label = 'Condensation temperature (K)'
        self.grad = QCPColorGradient('gpViridis')
        self.setLegendVisible(True)

    def updatePlot(self):
        self.clearGraphs()
        self.clearItems()

        if not hasattr(self.settingsWidget, 'sels') or not self.settingsWidget.sels:
            print('Did not update JacksonPlot because no selections...')
            return

        try:
            sel = self.settingsWidget.sels[0]
            isElements = sel.property('Internal element').split(',')

            channelName = self.settingsWidget.selectedChannelNames[0]
            self.left().label = f'{channelName} concentration (ppm)'
        except:
            return

        if not isElements or not isElements[0]:
            return

        thistc = data.elements[data.timeSeries(channelName).property('Element')]['Tcond_Lodders']
        tc = np.array([data.elements[data.timeSeries(ise).property('Element')]['Tcond_Lodders'] for ise in isElements])
        #tc = (thistc - tc)/tc

        cps = data.timeSeries(f'{channelName}_CPS')
        channelSurface = self.settingsWidget.calibration.surface(channelName, inv=True)

        conc = np.empty(len(tc))

        for i, ise in enumerate(isElements):
            try:
                isChannel = data.timeSeries(f'{ise}_CPS')
                isSurface = self.settingsWidget.calibration.surface(ise, inv=True)
                sqPPM = channelSurface(sel.midTimeInSec, data.result(sel, cps).value())
                rmPPM = data.referenceMaterialData(sel.group().name)[data.timeSeries(ise).property('Element')].valueInPPM()
                isPPM = isSurface(sel.midTimeInSec, data.result(sel, isChannel).value())
                conc[i] = rmPPM*sqPPM/isPPM
            except Exception as e:
                print(e)

        for i, ise in enumerate(isElements):
            g = self.addGraph()
            g.setData([tc[i]], [conc[i]])
            g.setLineStyle('lsNone')
            color = self.grad.color(i, 0, len(conc))
            g.setScatterStyle('ssDisc', 6, color, color)
            g.name = isElements[i]

        self.addStraightLine([thistc, conc.min()], [thistc, conc.max()])

        self.rescaleAxes()
        self.bottom().range = QCPRange(500, 1800)
        self.left().scaleRange(1.5)
        self.replot()

class ExtModelDialog(QDialog):

    def __init__(self, calibration, selectedChannels, parent=None):
        QDialog.__init__(self, parent)
        self.calibration = calibration
        self.calibration.updateBlocks()
        self.selectedChannels = selectedChannels
        self.blockSplines = {}
        try:
            self.setWindowTitle(f"Sensitivity model for {', '.join([ch.name for ch in self.selectedChannels])}")
        except:
            pass

        self.setLayout(QVBoxLayout())
        mainLayout = QHBoxLayout()
        leftWidget = QWidget(self)
        leftWidget.setLayout(QFormLayout())
        leftWidget.setFixedWidth(225)
        leftWidget.layout().setContentsMargins(0, 0, 0, 0)
        self.layout().addLayout(mainLayout)

        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, self)
        bb.accepted.connect(self.checkAndAccept)  #Can't exit Dialog without setting abundance.
        bb.rejected.connect(self.reject)

        self.layout().addWidget(bb)
        self.channelsTable = QTableWidget(self)
        self.allChannels = [c.name for c in data.timeSeriesList(data.Input) if 'Total' not in c.name]
        self.channelsTable.setRowCount(len(self.allChannels))
        self.channelsTable.setColumnCount(2)
        self.channelsTable.setHorizontalHeaderLabels(['Channel', 'Abundance'])
        self.channelsTable.horizontalHeader().setStretchLastSection(True)
        self.channelsTable.verticalHeader().setVisible(False)
        self.channelsTable.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.channelsTable.setSelectionMode(QAbstractItemView.ExtendedSelection)

        for ri, name in enumerate(self.allChannels):
            self.channelsTable.setItem(ri, 0, QTableWidgetItem(name))
            self.channelsTable.item(ri, 0).setFlags(self.channelsTable.item(ri, 0).flags() & ~Qt.ItemIsEditable)
            element = data.timeSeries(name).property('Element')
            try:
                mass = int(data.timeSeries(name).property('Mass'))
                ab = [i for i in data.elements[element]['isotopes'] if i['massNumber'] == mass][0]['abundance']
            except Exception as e:
                ab = 1
            self.channelsTable.setItem(ri, 1, QTableWidgetItem(f'{ab:.6f}'))
            self.channelsTable.item(ri, 0).setData(Qt.UserRole, True)
            if data.timeSeries(name).property('External standard') == 'Model':
                # If this is a channel using a model, highlight the items differently
                self.channelsTable.item(ri, 0).setData(Qt.UserRole, False)
                self.channelsTable.item(ri, 0).setBackground(Qt.black)
                self.channelsTable.item(ri, 1).setBackground(Qt.black)
                self.channelsTable.item(ri, 0).setFlags(self.channelsTable.item(ri, 0).flags() & ~Qt.ItemIsSelectable)
                self.channelsTable.item(ri, 1).setFlags(self.channelsTable.item(ri, 1).flags() & ~Qt.ItemIsSelectable)

            if 'ModelChannels' in self.selectedChannels[0].properties():
                if name in self.selectedChannels[0].property('ModelChannels'):
                    index = self.channelsTable.model().index(ri, 0)
                    self.channelsTable.selectionModel().select(index, QItemSelectionModel.Rows | QItemSelectionModel.Select)

        self.plot = Plot(self)
        self.plot.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        leftWidget.layout().addRow(self.channelsTable)
        self.fitComboBox = QComboBox(self)
        self.fitComboBox.addItems([
            'MeanMean',
            'MeanMedian',
            'LinearFit',
            'WeightedLinearFit',
            'StepLinear',
            'StepForward',
            'StepBackward',
            'StepAverage',
            'Nearest',
            'Akima',
            'Spline_NoSmoothing',
            'Spline_Smooth1',
            'Spline_Smooth2',
            'Spline_Smooth3',
            'Spline_Smooth4',
            'Spline_Smooth5',
            'Spline_Smooth6',
            'Spline_Smooth7',
            'Spline_Smooth8',
            'Spline_Smooth9',
            'Spline_Smooth10',
            'Spline_AutoSmooth'
        ])
        self.fitComboBox.setCurrentText(self.fitType())
        self.fitComboBox.textActivated.connect(self.updateFitType)
        leftWidget.layout().addRow('Fit type', self.fitComboBox)
        mainLayout.addWidget(leftWidget)
        mainLayout.addWidget(self.plot)
        self.resize(1200, 600)

        self.plot.left().label = 'Abundance corrected CPS/ppm'
        self.plot.bottom().label = 'Channel (m/z)'
        self.colorGrad = QCPColorGradient('gpSpectrum')

        self.channelsTable.selectionModel().selectionChanged.connect(self.processSelectionChanged)
        self.channelsTable.cellChanged.connect(self.processCellChanged)

        self.showLegend = False

        menu = self.plot.contextMenu()
        menu.addSeparator()
        la = menu.addAction('Legend')
        la.setCheckable(True)
        la.setChecked(self.showLegend)
        la.triggered.connect(self.toggleLegend)

        self.updatePlot()


    def checkAndAccept(self):
        #Check that selected channels have an abudance set:
        for ch in self.selectedChannels:
            item = self.channelsTable.item(
                self.channelsTable.findItems(ch.name, Qt.MatchExactly)[0].row(), 1)
            if float(item.text()) < 0.000001:
                self.abunNotice = OverlayMessage(self.plot,
                                                 'Error',
                                                 f'Please set the isotopic abundance of {ch.name} before closing this dialog',
                                                 3, 0.8, 40, 50)
                self.abunNotice.setCloseButtonVisible(False)
                self.abunNotice.show()
                return

        self.accept()

    def toggleLegend(self, b):
        self.showLegend = b
        # Get axis ranges so that toggling the legend doesn't reset zoom
        l_range = self.plot.left().range
        b_range = self.plot.bottom().range
        self.updatePlot()
        self.plot.left().setRange(l_range)
        self.plot.bottom().setRange(b_range)
        self.plot.replot()

    def fitType(self):
        try:
            fitType = [ch.property('ModelFitType') for ch in self.selectedChannels][0]
            if not fitType:
                fitType = 'LinearFit'
        except:
            fitType = 'LinearFit'

        return fitType

    def updateFitType(self, fit_name):
        for ch in self.selectedChannels:
            ch.setProperty('ModelFitType', fit_name)

        self.updatePlot()

    def processCellChanged(self, row, col):
        # Only process changes in col == 1 == abundance
        if col != 1:
            return

        channel_name = self.channelsTable.item(row, 0).text()
        print(f'Abundance changed for {channel_name}')

    def processSelectionChanged(self):
        channel_names = [index.data() for index in self.channelsTable.selectionModel().selectedRows()]
        for ch in self.selectedChannels:
            ch.setProperty('ModelChannels', channel_names)
        self.updatePlot()

    def updatePlot(self):
        print('Updating plot...')
        self.plot.clearGraphs()
        self.plot.clearItems()

        channelsForFit = [index.data() for index in self.channelsTable.selectionModel().selectedRows()]
        iForFit = [index.row() for index in self.channelsTable.selectionModel().selectedRows() if data.timeSeries(index.data()).property('External standard') != 'Model']
        iForFit.sort()

        channelMasses = np.array([int(data.timeSeries(n).property('Mass')) for n in self.allChannels], dtype=np.float)
        xOut = np.linspace(channelMasses.min(), channelMasses.max(), 1000)
        channelAbundances = np.array([float(self.channelsTable.item(ri, 1).text()) for ri in range(self.channelsTable.rowCount)])

        def sensForChannel(block, channel_name):
            if data.timeSeries(channel_name).property('External standard') == 'Model':
                return np.nan
            return block.slope(channel_name)

        # Quick note for future selves: can't plot channel sensitivities in one go because it will have to update if the user
        # changes any of the isotopic abundances in the table. So, can't put main plotting in __init__() for example
        for bi, block in enumerate(self.calibration.blocks):
            block_graph = self.plot.addGraph()
            block_graph.setName(f"Block {bi + 1}")
            channelSens = np.array([sensForChannel(block, cn) for cn in self.allChannels], dtype=np.float) / channelAbundances
            block_graph.setData(channelMasses[np.isfinite(channelSens)], channelSens[np.isfinite(channelSens)])
            color = self.colorGrad.color(bi, 0, len(self.calibration.blocks))
            block_graph.setScatterStyle('ssDisc', 6, color, color)
            block_pen = QPen(QColor(color))
            block_pen.setStyle(Qt.DotLine)
            block_graph.setPen(block_pen)

            xForFit = channelMasses[iForFit]
            yForFit = channelSens[iForFit]

            if len(xForFit) < 1:
                continue

            if len(xForFit) >= 1:
                yOut = data.spline(xForFit, yForFit, np.ones(len(yForFit)), self.fitType(), xOut)
                self.blockSplines[block.label] = (xOut, yOut)
                fitGraph = self.plot.addGraph()
                fitGraph.pen = QPen(color, 2)
                fitGraph.setData(xOut, yOut)
                fitGraph.removeFromLegend()

            # Also, highlight the masses used in the fit on the plot
            used_graph = self.plot.addGraph()
            used_graph.setData(xForFit, yForFit)
            used_graph.setScatterStyle('ssCircle', 10, Qt.red, Qt.transparent)
            used_graph.setLineStyle('lsNone')
            used_graph.removeFromLegend()


        for ch in self.selectedChannels:
            try:
                this_mass = float(ch.property('Mass'))
                y_Max = self.plot.left().range.upper()
                mass_line = self.plot.addStraightLine([this_mass,0], [this_mass, y_Max])
                mass_line.pen = QPen(Qt.DashLine)
            except Exception:
                # TODO: get mass from channel name???
                pass

        self.plot.rescaleAxes()
        # Going to get x maximum from channels because if there are a lot of channels in the actinides that are being
        # fitted, they won't be in the table. So get it from the dataManager
        try:
            last_ch_mass = channelMasses.max()
            self.plot.bottom().range = QCPRange(0, last_ch_mass + last_ch_mass * 0.1)
        except Exception:
            self.plot.bottom().range = QCPRange(0, self.plot.bottom().range.upper() + self.plot.bottom().range.upper() * 0.1)

        #self.plot.left().range = QCPRange(0, y_currentMax + y_currentMax * 0.1)

        self.plot.setLegendVisible(self.showLegend)
        self.plot.setLegendAlignment(Qt.AlignTop | Qt.AlignLeft)

        self.plot.replot()

    def blockSensitivities(self, channel_name):
        from scipy.interpolate import interp1d
        # Need to get the abundance for this
        try:
            mass = float(data.timeSeries(channel_name).property('Mass'))
            table_row = self.channelsTable.findItems(channel_name, Qt.MatchFixedString)[0].row()
            abund = float(self.channelsTable.item(table_row, 1).text())
        except:
            print(f'Could not calculate block sensitivities for {channel_name}')
            mass = np.nan
            abund = np.nan

        bs = {}
        for bi, block in enumerate(self.calibration.blocks):
            x, y = self.blockSplines[block.label]
            f = interp1d(x, y, 'linear', fill_value='extrapolate')
            S = f(mass)*abund
            bs[block.label] = S

        return bs


class DataStatus(Flag):
    Ready = 0
    NoBaselineSelections = auto()
    NoRMSelections = auto()
    NoInputs = auto()
    NoCPS = auto()
    NoInternalStandard = auto()
    NoExternalStandard = auto()
    MissingChannelMetadata = auto()


class SettingsWidget(QWidget):

    def __init__(self, parent=None):
        super().__init__(parent)
        self.calibration = Calibration()

        settings = QSettings()
        self.ui_path = settings.value("Paths/DataReductionSchemesPath")
        self.ui_file = QFile(self.ui_path + "/3d_trace_elements.ui")

        self.setLayout(QVBoxLayout())
        if not self.ui_file.open(QIODevice.ReadOnly):
            raise RuntimeError('Could not load settings ui')

        ui = QUiLoader().load(self.ui_file, self)
        self.layout().addWidget(ui)
        self.layout().setContentsMargins(0, 0, 0, 0)

        self.extTable = ui.findChild(QTableView, 'externalTable')
        self.extFilter = ui.findChild(QLineEdit, 'filterES')
        self.intGroupBox = ui.findChild(QGroupBox, 'ISGroupBox')
        self.intTable = ui.findChild(QTableView, 'internalTable')
        self.intFilter = ui.findChild(QLineEdit, 'filterIS')
        self.tabWidget = ui.findChild(QTabWidget, 'tabWidget')
        self.setExtButton = ui.findChild(QToolButton, 'setExternalButton')
        self.setIntButton = ui.findChild(QToolButton, 'internalElementButton')
        self.setIntValueButton = ui.findChild(QToolButton, 'internalValueButton')
        self.setIntUnitsButton = ui.findChild(QToolButton, 'internalUnitsButton')
        self.importValuesButton = ui.findChild(QToolButton, 'importValuesButton')
        self.setIndexChannelCBox = ui.findChild(QComboBox, 'indexChannel')
        self.maskGroupBox = ui.findChild(QGroupBox, 'maskGroupBox')
        self.maskChannelCBox = ui.findChild(QComboBox, "maskChannel")
        self.maskMethodCBox = ui.findChild(QComboBox, "maskMethodComboBox")
        self.maskLineEdit = ui.findChild(QLineEdit, "maskCutoff")
        self.maskTrimSB = ui.findChild(QDoubleSpinBox, "maskTrim")
        self.maskFrame = ui.findChild(QFrame, 'maskFrame')
        self.normalizeExtCheckBox = ui.findChild(QCheckBox, 'normalizeExtCheckBox')
        self.fitMethodComboBox = ui.findChild(QComboBox, 'fitMethodComboBox')
        self.externalSplitter = ui.findChild(QSplitter, 'externalSplitter')
        self.normToComboBox = ui.findChild(QComboBox, 'normToComboBox')
        self.splineTypeComboBox = ui.findChild(QComboBox, 'splineTypeComboBox')
        self.splineTypeComboBox.addItems([
            'MeanMean',
            'MeanMedian',
            'LinearFit',
            'WeightedLinearFit',
            'StepLinear',
            'StepForward',
            'StepBackward',
            'StepAverage',
            'Nearest',
            'Akima',
            'Spline_NoSmoothing',
            'Spline_Smooth1',
            'Spline_Smooth2',
            'Spline_Smooth3',
            'Spline_Smooth4',
            'Spline_Smooth5',
            'Spline_Smooth6',
            'Spline_Smooth7',
            'Spline_Smooth8',
            'Spline_Smooth9',
            'Spline_Smooth10',
            'Spline_AutoSmooth'
        ])
        self.criteriaButton = ui.findChild(QToolButton, 'criteriaButton')
        self.throughZeroButton = ui.findChild(QToolButton, 'throughZeroButton')
        self.fracButton = ui.findChild(QToolButton, 'fracToolButton')
        self.modelComboBox = ui.findChild(QComboBox, 'modelComboBox')
        self.bsMethodComboBox = ui.findChild(QComboBox, 'bsMethodComboBox')
        self.bsChannelComboBox = ui.findChild(QComboBox, 'bsChannelComboBox')
        self.bsLineEdit = ui.findChild(QLineEdit, 'bsLineEdit')
        self.bsFrame = ui.findChild(QFrame, 'bsFrame')
        self.affinityButton = ui.findChild(QToolButton, 'affinityButton')
        self.affinityCorrectionCheckBox = ui.findChild(QCheckBox, 'affinityCorrectionCheckBox')

        self.throughZeroButton.setIcon(CUI().icon('bullseye'))
        self.throughZeroButton.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        self.throughZeroButton.setToolTip("Force calibration curve through 0")

        self.fracButton.setIcon(CUI().icon('magic'))
        self.fracMenu = QMenu(self)
        fracActions = [self.fracMenu.addAction(t) for t in ['None', 'Linear', 'Spline']]
        self.fracButton.setMenu(self.fracMenu)
        self.fracButton.setPopupMode(QToolButton.DelayedPopup)

        self.popoutButton = QToolButton(self.tabWidget)
        self.popoutButton.setIcon(CUI().icon('externallink'))
        self.popoutButton.setCheckable(True)
        self.popoutButton.setChecked(False)
        self.popoutButton.clicked.connect(self.togglePlotEmbedding)
        self.tabWidget.installEventFilter(self)
        self.tabWidget.setCornerWidget(self.popoutButton)
        self.tabWidget.setWindowTitle('3D Trace Elements')
        self.tabWidget.setAttribute(Qt.WA_DeleteOnClose, False)

        timeSeriesNames = data.timeSeriesNames(data.Input)
        defaultChannelName = ""
        if timeSeriesNames:
            defaultChannelName = timeSeriesNames[0]

        drs.setDefaultSetting("IndexChannel", defaultChannelName)
        drs.setDefaultSetting("Mask", False)
        drs.setDefaultSetting("MaskChannel", defaultChannelName)
        drs.setDefaultSetting('MaskMethod', 'Laser log')
        drs.setDefaultSetting("MaskCutoff", 0.1)
        drs.setDefaultSetting("MaskTrim", 0.0)
        drs.setDefaultSetting("NormalizeExternals", True)
        drs.setDefaultSetting("UseIntStds", True)
        drs.setDefaultSetting('SplineType', 'Spline_AutoSmooth')
        drs.setDefaultSetting('MasterExternal', '')
        bsDefault = 'Laser log' if len(data.laserLogSamples()) > 0 else 'Cutoff threshold'
        drs.setDefaultSetting('BeamSecondsMethod', bsDefault)
        drs.setDefaultSetting('BeamSecondsChannel', defaultChannelName)
        drs.setDefaultSetting('BeamSecondsValue', 1000.)
        drs.setDefaultSetting('AffinityCorrection', False)
        drs.setDefaultSetting('BlockFindingMethod', 'Simple')
        drs.setDefaultSetting('NClusters', -1)

        drs.setSetting('AffinityCorrection', False)

        self.maskTrimSB.setMinimum(-1E5)
        self.maskTrimSB.setMaximum(1E5)

        self.setupChannelCBs()
        self.setupExt()
        self.setupInt()
        self.setup3dPlot()
        self.setupBlockPlot()
        self.setupFitsPlot()
        self.setupFitParamsPlot()
        self.setupFractionationPlot()

        self.blockRMs = []

        # Restore settings
        print(f'Restoring settings {drs.settings()}')
        self.normalizeExtCheckBox.setChecked(drs.setting('NormalizeExternals'))
        self.setIndexChannelCBox.currentText = drs.setting('IndexChannel')
        self.maskLineEdit.setText(str(drs.setting('MaskCutoff')))
        self.maskTrimSB.setValue(float(drs.setting('MaskTrim')))
        self.maskChannelCBox.currentText = drs.setting('MaskChannel')
        self.maskMethodCBox.currentText = drs.setting('MaskMethod')
        self.maskFrame.setVisible('log' not in drs.setting('MaskMethod'))
        self.intGroupBox.setChecked(bool(drs.setting('UseIntStds')))
        self.splineTypeComboBox.currentText = drs.setting('SplineType')
        self.maskGroupBox.setChecked(drs.setting('Mask'))
        self.normToComboBox.currentText = drs.setting('MasterExternal')
        self.bsMethodComboBox.currentText = drs.setting('BeamSecondsMethod')
        self.bsChannelComboBox.currentText = drs.setting('BeamSecondsChannel')
        self.bsLineEdit.setText(drs.setting('BeamSecondsValue'))
        self.bsFrame.setVisible('log' not in drs.setting('BeamSecondsMethod').lower())
        self.affinityCorrectionCheckBox.setChecked(drs.setting('AffinityCorrection'))

        self.sels = []

        # Connections
        self.setIndexChannelCBox.activated.connect(lambda t: drs.setSetting("IndexChannel", self.setIndexChannelCBox.currentText))
        self.maskGroupBox.toggled.connect(lambda b: drs.setSetting('Mask', b))
        self.maskMethodCBox.activated.connect(lambda t: drs.setSetting('MaskMethod', self.maskMethodCBox.currentText))
        self.maskMethodCBox.activated.connect(lambda t: self.maskFrame.setVisible(t))
        self.maskChannelCBox.activated.connect(lambda t: drs.setSetting("MaskChannel", self.maskChannelCBox.currentText))
        self.maskLineEdit.textEdited.connect(lambda t: drs.setSetting("MaskCutoff", float(t)))
        self.maskTrimSB.valueChanged.connect(lambda t: drs.setSetting("MaskTrim", float(t)))
        self.normalizeExtCheckBox.toggled.connect(lambda b: self.updateDRSSetting('NormalizeExternals', b))
        self.normToComboBox.activated.connect(lambda t: self.updateDRSSetting('MasterExternal', self.normToComboBox.currentText))
        self.splineTypeComboBox.activated.connect(lambda t: self.updateDRSSetting('SplineType', self.splineTypeComboBox.currentText))
        self.intGroupBox.toggled.connect(lambda b: drs.setSetting('UseIntStds', b))
        self.criteriaButton.clicked.connect(self.editCriteria)
        self.throughZeroButton.clicked.connect(self.toggleThroughZero)
        self.importValuesButton.clicked.connect(self.importValues)
        self.fracButton.clicked.connect(self.toggleFractionation)
        self.fracMenu.triggered.connect(self.toggleFractionation)
        data.dataChanged.connect(self.setupChannelCBs)
        self.modelComboBox.activated.connect(self.changeModel)
        self.bsMethodComboBox.activated.connect(lambda t: drs.setSetting('BeamSecondsMethod', self.bsMethodComboBox.currentText))
        self.bsMethodComboBox.activated.connect(lambda t: self.bsFrame.setVisible(t in [1, 2]))
        self.bsChannelComboBox.activated.connect(lambda t: drs.setSetting('BeamSecondsChannel', self.bsChannelComboBox.currentText))
        self.bsLineEdit.textEdited.connect(lambda t: drs.setSetting('BeamSecondsValue', float(t)))
        drs.finished.connect(lambda: drs.setProperty('isRunning', False))

        drs.finished.connect(lambda: self.intModel.refreshSelections())

        self.intTable.installEventFilter(self)

        # Create notices that may or may not be hidden
        self.extNotice = OverlayMessage(self.extTable, 'Warning', 'Baseline subtracted channels are missing.', 0, 0.8, 40, 10)
        self.extNotice.setCloseButtonVisible(False)
        self.blsButton = QToolButton(self.extNotice)
        self.blsButton.setText('Calculate now!')
        self.blsButton.clicked.connect(self.extModel.refreshChannels)
        self.extNotice.addButton(self.blsButton)
        self.extNotice.hide()

        self.noSelsNotice = OverlayMessage(self.extTable, 'Warning', 'There are no selections for one of the External Standards',
                                           0, 0.8, 40, 10)
        self.noSelsNotice.setCloseButtonVisible(False)
        self.noSelsNotice.hide()

        self.intNotice = OverlayMessage(self.intTable, 'Warning', 'Elements and values must be specified for > 0 selections.', 0, 0.8, 40, 10)
        self.intNotice.setCloseButtonVisible(False)
        self.intNotice.hide()

        self.mdNotice = OverlayMessage(self.extTable, 'Warning', 'There is missing metadata for some channels.', 0, 0.8, 40, 10)
        self.mdNotice.setCloseButtonVisible(False)
        self.mdNotice.hide()

        self.updateStatus()

    def updateStatus(self):
        self.status = DataStatus.Ready
        if len(data.timeSeriesNames(data.Input)) == 0:
            self.status = self.status | DataStatus.NoInputs
        if len(data.selectionGroupList(data.Baseline)) == 0 or len(data.selectionGroupList(data.Baseline)[0].selections()) == 0:
            self.status = self.status | DataStatus.NoBaselineSelections
        if len(data.timeSeriesNames(data.Intermediate, {'DRSType': 'BaselineSubtracted'})) == 0:
            self.status = self.status | DataStatus.NoCPS
        if len(data.selectionGroupList(data.ReferenceMaterial)) == 0 or len(data.selectionGroupList(data.ReferenceMaterial)[0].selections()) == 0:
            self.status = self.status | DataStatus.NoRMSelections
        if len(self.externalsInUse()) == 0:
            self.status = self.status | DataStatus.NoExternalStandard
        if len(self.internalsInUse()) == 0:
            self.status = self.status | DataStatus.NoInternalStandard
        #print(f'Updated status to {self.status}')
        self.maybeShowExtNotices()


    def eventFilter(self, obj, event):
        if obj == self.tabWidget and event.type()== QEvent.Close:
            print(f'Caught close event on tabWidget')
            event.ignore()
            return True

        if obj == self.intTable and event.type() == QEvent.KeyPress and event.matches(QKeySequence.Paste):
            text = QApplication.clipboard().text()
            values = re.split('\s+|,', text)
            values = [float(v) for v in values]

            if len(values) > 1: # If multiple values pasted, use them up
                indCount = 0
                for ind in self.intTable.selectionModel().selectedRows(3):
                    if indCount >= len(values):
                        break

                    self.intModel.setData(ind, values[indCount])
                    indCount += 1
            else: # If only 1 value pasted, paste it to all selected rows
                for ind in self.intTable.selectionModel().selectedRows(3):
                    self.intModel.setData(ind, values[0])


        return False

    def setupChannelCBs(self):        
        cbs = [self.setIndexChannelCBox, self.maskChannelCBox, self.bsChannelComboBox]

        for cb in cbs:
            ct = cb.currentText
            cb.clear()
            cb.addItems(data.timeSeriesNames(data.Input))

            if cb.findText(ct) >= 0:
                cb.currentText = ct

    def externalsInUse(self):
        groupNames = list(itertools.chain(*[str(self.extModel.index(r, 1).data()).split(',') for r in range(self.extModel.rowCount())]))
        groupNames = list(set(groupNames))
        try:
            groupNames.remove('Model')
            groupNames.remove('None')
        except:
            pass
        return groupNames

    def internalsInUse(self):
        isElements = [str(s.property('Internal element')) for s in self.intModel.selections]
        isElements = list(set(isElements))
        if 'None' in isElements:
            isElements.remove('None')
        if '' in isElements:
            isElements.remove('')
        return isElements


    def updateNormOptions(self):
        current = self.normToComboBox.currentText
        self.normToComboBox.clear()
        self.normToComboBox.addItems(self.externalsInUse())
        index = self.normToComboBox.findText(current)
        if index >= 0:
            self.normToComboBox.currentIndex = index

    def setupExt(self):
        self.extTable.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.extTable.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.extFilterModel = QSortFilterProxyModel(self)
        self.extModel = ExternalsModel(self)
        self.extFilterModel.setSourceModel(self.extModel)
        self.extTable.setModel(self.extFilterModel)
        self.extFilter.textEdited.connect(self.extFilterModel.setFilterFixedString)
        self.extTable.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self.extTable.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        self.extTable.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeToContents)
        self.extTable.horizontalHeader().setSectionResizeMode(3, QHeaderView.ResizeToContents)

        self.extTable.setContextMenuPolicy(Qt.CustomContextMenu)
        self.extTable.customContextMenuRequested.connect(self.showExtContextMenu)

        self.rmMenu = ReferenceMaterialsMenu(self)
        self.setExtButton.setMenu(self.rmMenu)
        self.setExtButton.setPopupMode(QToolButton.InstantPopup)
        self.setExtButton.setIcon(CUI().icon('trophy'))
        self.setExtButton.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)

        self.extModel.dataChanged.connect(lambda a: self.updateAffected())

        self.extModel.dataChanged.connect(lambda a: self.updateBlocks())
        self.extTable.selectionModel().selectionChanged.connect(self.processExtSelection)
        self.rmMenu.rmsChanged.connect(self.extModel.updateData)
        self.extModel.dataChanged.connect(self.updateNormOptions)
        self.updateNormOptions()
        if 'MasterExternal' in drs.settings():
            self.normToComboBox.currentText = drs.setting('MasterExternal')

        self.extModel.dataChanged.connect(self.processExtSelection)
        self.extTable.setItemDelegate(ExternalsDelegate(self.extTable))
        #data.dataChanged.connect(lambda: self.extModel.refreshChannels())
        data.dataChanged.connect(self.updateStatus)

    def showExtContextMenu(self, point):
        index = self.extTable.indexAt(point)
        rmsIndex = self.extModel.index(index.row(), 1)
        rms = str(rmsIndex.data())
        if rms != 'Model':
            print('This channel does not use a model, so not showing menu...')
            return
        menu = QMenu()
        action = menu.addAction('Configure Model')
        action.triggered.connect(self.configureExtModel)
        menu.exec(self.extTable.mapToGlobal(point))

    def configureExtModel(self):
        try:
            w = ExtModelDialog(self.calibration, self.selectedChannels)
        except MissingRMGroupError:
            self.status = DataStatus.NoRMSelections
            self.maybeShowExtNotices()
            return

        if w.exec() != QDialog.Accepted:
            return

        # Apply dialog settings to selected channels...
        for channel in self.selectedChannels:
            bs = w.blockSensitivities(channel.name)
            channel.setProperty('ModelSensitivities', bs)

        self.processExtSelection()

    def maybeShowExtNotices(self):
        if drs.property('isRunning'):
            self.extNotice.hide()
            self.intNotice.hide()
            return

        if self.status & DataStatus.NoCPS:
            self.extNotice.show()
        else:
            self.extNotice.hide()

        if self.status & DataStatus.NoRMSelections:
            self.noSelsNotice.show()
        else:
            self.noSelsNotice.hide()

        if self.status & DataStatus.NoInternalStandard and drs.setting('UseIntStds'):
            self.intNotice.show()
        else:
            self.intNotice.hide()

        if self.status & DataStatus.MissingChannelMetadata:
            self.mdNotice.show()
        else:
            self.mdNotice.hide()

    def setupInt(self):
        self.intTable.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.intTable.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.intFilterModel = QSortFilterProxyModel(self)
        self.intModel = InternalsModel(self)
        self.intFilterModel.setSourceModel(self.intModel)
        self.intTable.setModel(self.intFilterModel)
        self.intFilter.textEdited.connect(self.intFilterModel.setFilterFixedString)
        self.intTable.selectionModel().selectionChanged.connect(self.processIntSelection)

        self.elementMenu = ChannelsMenu(self)
        self.setIntButton.setMenu(self.elementMenu)
        self.setIntButton.setPopupMode(QToolButton.InstantPopup)

        self.setIntValueButton.clicked.connect(self.getInternalValue)
        self.elementMenu.channelsChanged.connect(self.intModel.updateData)
        self.elementMenu.channelsChanged.connect(lambda: self.calibration.clearFractionationCache())
        self.elementMenu.channelsChanged.connect(lambda: self.fracPlot.updatePlot())
        self.elementMenu.channelsChanged.connect(self.updateStatus)

        self.affinityMenu = ChannelsMenu(self, propName = 'Affinity elements')
        self.affinityMenu.removeAction(self.affinityMenu.actions()[-1])
        majorsAction = self.affinityMenu.addAction('Majors')
        majorsAction.triggered.connect(self.setAffinitiesToMajors)
        self.affinityButton.setMenu(self.affinityMenu)
        self.affinityMenu.channelsChanged.connect(self.intModel.updateData)
        self.affinityButton.clicked.connect(lambda: self.intModel.refreshSelections())
        self.affinityCorrectionCheckBox.clicked.connect(lambda b: drs.setSetting('AffinityCorrection', b))

        affRMMenu = QMenu('Reference materials', self)
        self.affinityMenu.addMenu(affRMMenu)
        for rm in data.referenceMaterialNames():
            affRMMenu.addAction(rm)
        affRMMenu.triggered.connect(lambda a: self.setAffinitiesToMaterial(a.text))

        self.unitsMenu = QMenu(self)
        for u in ['ppm', 'ppb', 'wtpc', 'wtpc_oxide']:
            self.unitsMenu.addAction(u)
        self.setIntUnitsButton.setMenu(self.unitsMenu)
        self.setIntUnitsButton.setPopupMode(QToolButton.InstantPopup)
        self.unitsMenu.triggered.connect(self.processUnits)
        self.intTable.setItemDelegate(InternalsDelegate(self.intTable))
        data.selectionGroupsChanged.connect(lambda: self.intModel.refreshSelections())
        data.selectionGroupsChanged.connect(self.updateStatus)

    def setAffinitiesToMajors(self):
        allMajors = ['Si', 'Al', 'Ca', 'Mg', 'Na', 'K', 'Ti', 'Fe', 'Mn', 'P']
        majors = [ch.name for ch in data.timeSeriesList(data.Input) if ch.property('Element') in allMajors]
        sels_to_update = self.sels

        if len(sels_to_update) < 1:     #If no selections selected in table, apply to all
            sels_to_update = self.intModel.selections

        for s in sels_to_update:
            s.setProperty('Affinity elements', ','.join(majors))

        self.intModel.updateData(self.sels)

    def setAffinitiesToMaterial(self, material):
        sels_to_update = self.sels

        if len(sels_to_update) < 1:     #If no selections selected in table, apply to all
            sels_to_update = self.intModel.selections

        for s in sels_to_update:
            s.setProperty('Affinity elements', material)

        self.intModel.updateData(self.sels)

    def setup3dPlot(self):
        self.plot3d = Plot3d(self)
        self.plot3d.setLabelFont(QFont('Helvetica', 16, QFont.Bold))
        self.tabWidget.addTab(self.plot3d, '3D')        
        self.resetViewButton = OverlayButton(self.plot3d, 'TopLeft', 10, 10)
        self.resetViewButton.setIcon(CUI().icon('cube'))
        self.resetViewButton.setFixedSize(25, 25)
        self.resetViewButton.clicked.connect(self.reset3dView)
        self.reset3dView()

    def reset3dView(self):
        self.plot3d.setView(1, [1, 1, 1], [30, 0, 45])

    def togglePlotEmbedding(self, b):
        if b:
            # Out
            label = QLabel('Visualization in external window')
            label.setAlignment(Qt.AlignCenter)
            CUI().splitterReplace(self.externalSplitter, 1, label)
            self.tabWidget.show()
            self.tabWidget.resize(600, 600)
            try:
                self.update3d(self.selectedChannelNames[0])
            except:
                pass
        else:
            # In
            CUI().splitterReplace(self.externalSplitter, 1, self.tabWidget)
            try:
                self.update3d(self.selectedChannelNames[0])
            except:
                pass

    def setupBlockPlot(self):
        self.blockPlot = BlockPlot(self)
        self.tabWidget.addTab(self.blockPlot, 'Blocks')

    def setupFitsPlot(self):
        self.fitsPlot = FitsPlot(self)
        self.tabWidget.addTab(self.fitsPlot, 'Fits')

    def setupFitParamsPlot(self):
        self.fitParamsPlot = FitParamsPlot(self)
        self.tabWidget.addTab(self.fitParamsPlot, 'Fit Params')

    def setupFractionationPlot(self):
        self.fracPlot = FractionationPlot(self)
        self.tabWidget.addTab(self.fracPlot, 'Fractionation')

        #self.jackPlot = JacksonPlot(self)
        #self.tabWidget.addTab(self.jackPlot, 'Condensation T')

    def updateBlocks(self):        
        if not drs.setting('MasterExternal') or drs.setting('MasterExternal') not in self.externalsInUse():
            print(f'Master external is invalid...')
            if self.externalsInUse() and len(self.externalsInUse()) > 0:
                me = self.externalsInUse()[0]
                drs.setSetting('MasterExternal', me)
                self.normToComboBox.currentText = me
                print(f'   ... setting it to {me}.')
            else:
                print('   ... because there are no externals in use.')

        if self.blockRMs != self.externalsInUse():            
            try:
                self.calibration.updateBlocks()
            except MissingRMGroupError:
                return

            self.blockRMs = self.externalsInUse()

    def compileMeasurementData(self, channelName, ex, ey, ez):
        mx = np.array([]); my = np.array([]); mz = np.array([])
        msx = np.array([]); msy = np.array([]); msz = np.array([])
        grad = QCPColorGradient('gpViridis')
        colors = []

        for bi, block in enumerate(self.calibration.blocks):
            df = block.dataFrameForChannel(channelName)
            mx = np.append(mx, df['sel_mid_time'])
            msx = np.append(msx, df['sel_duration'])
            my = np.append(my, df[f'{channelName}_RMppm'])
            msy = np.append(msy, df[f'{channelName}_RMppm_Uncert'])
            mz = np.append(mz, df[f'{channelName}'])
            msz = np.append(msz, df[f'{channelName}_Uncert'])
            c = grad.color(bi, 0, len(self.calibration.blocks))
            c.setAlpha(120)
            colors += [c]*len(df)

        # Normalize to the range being sent for surface
        msx = (msx/mx); msy = msy/my; msz = msz/mz
        mx = (mx-ex[0])/(ex[1] - ex[0])
        my = (my-ey[0])/(ey[1] - ey[0])
        mz = (mz-ez[0])/(ez[1] - ez[0])
        msx = msx*mx; msy = msy*my; msz = msz*mz
        return mx, msx, my, msy, mz, msz, colors

    def update3d(self, channelName, n=100):
        try:
            cpsChannel = data.timeSeries(f'{channelName}_CPS')
            minSlope = np.min([abs(block.slope(channelName)) for block in self.calibration.blocks])
            minInt = np.min([block.intercept(channelName) for block in self.calibration.blocks])
            # cps = slope*ppm + intercept
            af, _ = calculateRelativeYields()
            aff = 1./min(af.values()) # To compensate for low yield
            maxppm = aff*np.nanmax( (cpsChannel.data() - minInt)/(minSlope) )
            x = np.linspace(cpsChannel.time().min(), cpsChannel.time().max(), n)
            y = np.linspace(0, maxppm, n)
            X,Y = np.meshgrid(x,y)
            Z = self.surface(X,Y)
            Z = Z.astype(np.float32)
            self.plot3d.setSurfaceData(Z)
            self.plot3d.setXTitle('Time')
            self.plot3d.setYTitle('Concentration')
            self.plot3d.setZTitle('Intensity')
            for ax in ['X1', 'X2', 'X3', 'X4', 'Y1', 'Y2', 'Y3', 'Y4', 'Z1', 'Z2', 'Z3', 'Z4']:
                self.plot3d.setTickLabels(ax, [0, 1], ['Min', 'Max'])

            self.plot3d.setOrtho(True)
            self.plot3d.setMeasurementData(*self.compileMeasurementData(channelName, (x.min(), x.max()), (y.min(), y.max()), (Z.min(), Z.max())))
        except Exception as e:
            # if there is any problem getting the above data, clear the 3d plot
            # probably they have not set an external for the selected group
            self.plot3d.clearData()

    def updateAffected(self):
        # If one of the channels modified is part of an IS we also need to
        # clear all the stored fractionation and surface caches
        selectedNames = [sr.data(Qt.UserRole).name for sr in self.extTable.selectionModel().selectedRows()]
        print(f'Update Affected {selectedNames}')
        if any(name in ''.join(self.internalsInUse()) for name in selectedNames):
            print('... clearing caches!')
            self.calibration.clearFractionationCache()
            self.calibration.clearSurfaceCache()

        [self.calibration.updateSurface(name) for name in selectedNames]

        try:
            [self.calibration.updateFractionation(name) for name in selectedNames]
        except RuntimeError:
            return


    def processExtSelection(self):
        self.selectedChannels = [sr.data(Qt.UserRole) for sr in self.extTable.selectionModel().selectedRows()]
        if not self.selectedChannels:
            self.selectedChannels = [c for c in data.timeSeriesList(data.Input) if 'TotalBeam' not in c.name]

        self.selectedChannelNames = [c.name for c in self.selectedChannels]

        # Check that all channels have an element set so that it doesn't cause errors below
        for ch in data.timeSeriesList(data.Input):
            if ch.name == 'TotalBeam':
                continue
            if not ch.property('Element') or len(ch.property('Element')) < 1:
                print(f"There was an issue with {ch.name}")
                self.status = self.status | DataStatus.MissingChannelMetadata
                self.maybeShowExtNotices()
                return

        self.rmMenu.rmsForActiveChannels = []
        for c in self.selectedChannels:
            if c.property('External standard') is not None:
                self.rmMenu.rmsForActiveChannels.append(c.property('External standard'))

        self.rmMenu.activeChannels = self.selectedChannelNames
        if len(self.selectedChannels) < 1:
            return

        channel = self.selectedChannels[0]
        fitThroughZero = bool(channel.property('FitThroughZero'))
        frac = bool(channel.property('FractionationCorrection'))
        self.throughZeroButton.setChecked(fitThroughZero)
        self.fracButton.setChecked(frac)
        model = channel.property('Model')
        if model:
            self.modelComboBox.setCurrentText(model)

        if len(self.calibration.blocks) == 0:
           try:
               self.updateBlocks()
           except RuntimeError:
               self.status = DataStatus.NoRMSelections
               return

        if len(self.calibration.blocks) == 0:
            print('No blocks. Aborting update!')
            return

        self.blockPlot.updatePlot()
        self.fitsPlot.updatePlot()
        self.fitParamsPlot.updatePlot()
        self.fracPlot.updatePlot()
        #self.jackPlot.updatePlot()

        # Check to make sure we have RM values for this channel
        extStd = channel.property('External standard')
        if not self.calibration.blocks[0].dataFrame()[f'{channel.name}_RMppm'].notna().values.any() and extStd != 'Model':
            self.surface = None
        elif type(extStd) == str and (extStd.split(',')[0] in data.selectionGroupNames() or extStd == 'Model'):
            self.surface, _ = fitSurface(self.calibration.blocks, channel.name)
            # fitSurface may have changed the spline type, so update UI here:
            self.splineTypeComboBox.currentText = drs.setting('SplineType')
        else:
            self.surface = None

        self.update3d(channel.name)

    def processIntSelection(self):
        self.sels = [sr.data(Qt.UserRole) for sr in self.intTable.selectionModel().selectedRows()]
        if not self.sels:
            self.sels = self.intModel.selections

        self.elementMenu.setSelections(self.sels)
        self.affinityMenu.setSelections(self.sels)
        #self.jackPlot.updatePlot()

    def getInternalValue(self):
        d = QDialog()
        d.setWindowFlags(Qt.Popup)
        l = QHBoxLayout()
        d.setLayout(l)
        l.setContentsMargins(3,3,3,3)
        l.setSpacing(0)
        valueLineEdit = QLineEdit(d)
        l.addWidget(valueLineEdit)
        okButton = QToolButton(d)
        okButton.setIcon(CUI().icon('check'))
        okButton.clicked.connect(lambda b: d.accept())
        cancelButton = QToolButton(d)
        cancelButton.setIcon(CUI().icon('remove'))
        cancelButton.clicked.connect(lambda b: d.reject())
        l.addWidget(okButton)
        l.addWidget(cancelButton)

        okShortcut = QShortcut(QKeySequence(Qt.Key_Return), okButton)
        okShortcut.activated.connect(lambda: okButton.click())
        okShortcut2 = QShortcut(QKeySequence(Qt.Key_Enter), okButton)
        okShortcut2.activated.connect(lambda: okButton.click())
        cancelShortcut = QShortcut(QKeySequence(Qt.Key_Escape), cancelButton)
        cancelShortcut.activated.connect(lambda: cancelButton.click())

        valueLineEdit.setFocus()
        d.move(self.setIntValueButton.mapToGlobal(QPoint(0, self.setIntValueButton.height)))

        if d.exec_() != QDialog.Accepted:
            return

        sels = [sr.data(Qt.UserRole) for sr in self.intTable.selectionModel().selectedRows()]
        if not sels:
            sels = self.intModel.selections

        for sel in sels:
            if not sel.property('Internal element') or sel.property('Internal element') != 'Criteria':
                sel.setProperty('Internal value', float(valueLineEdit.text))
            else:
                sel.setProperty('Internal value', valueLineEdit.text)

        self.intModel.updateData()

    def processUnits(self, action):
        u = action.text
        sels = [sr.data(Qt.UserRole) for sr in self.intTable.selectionModel().selectedRows()]
        if not sels:
            sels = self.intModel.selections

        for sel in sels:
            sel.setProperty('Internal units', u)

        self.intModel.updateData()

    def editCriteria(self):
        d = CriteriaDialog()
        d.table.model().setCriteria(drs.setting('ISCriteria'))

        if d.exec_() == QDialog.Rejected:
            return

        drs.setSetting('ISCriteria', d.criteria())

    def updateDRSSetting(self, name, value):
        print(f'Updating DRS setting {name} to {value}')
        if drs.setting(name) == value:
            return
        drs.setSetting(name, value)
        self.processExtSelection()

    def toggleThroughZero(self, b):
        channels = [sr.data(Qt.UserRole) for sr in self.extTable.selectionModel().selectedRows()]
        if not channels:
            channels = [c for c in data.timeSeriesList(data.Input) if 'TotalBeam' not in c]

        for channel in channels:
            channel.setProperty('FitThroughZero', b)

        self.updateAffected()

        self.processExtSelection()
        self.extFilterModel.invalidate()

    def toggleFractionation(self, b):
        on = (isinstance(b, bool) and b) or (isinstance(b, QAction) and b.text != 'None')
        fit = 'Linear' if on else 'None'
        if on and isinstance(b, QAction):
            fit = b.text

        channels = [sr.data(Qt.UserRole) for sr in self.extTable.selectionModel().selectedRows()]
        if not channels:
            channels = [c for c in data.timeSeriesList(data.Input) if 'TotalBeam' not in c.name]

        for channel in channels:
            channel.setProperty('FractionationCorrection', on)
            channel.setProperty('FractionationFitType', fit)

        self.processExtSelection()
        self.extFilterModel.invalidate()

    def changeModel(self, act):
        channels = [sr.data(Qt.UserRole) for sr in self.extTable.selectionModel().selectedRows()]
        if not channels:
            channels = [c for c in data.timeSeriesList(data.Input) if 'TotalBeam' not in c.name]

        for channel in channels:
            channel.setProperty('Model', self.modelComboBox.currentText)

        self.processExtSelection()
        self.extFilterModel.invalidate()

    def importValues(self):
        lastRawDataPath = QSettings().value('paths/lastrawdatapath', QDir().home().absolutePath())
        fileName = QFileDialog.getOpenFileName(self, 'Open IS file', lastRawDataPath, 'IS file (*.txt *.csv)')

        if not fileName:
            return

        res = data.importISValues(fileName)
        print(res)



def settingsWidget():
    drs.setSettingsWidget(SettingsWidget())
