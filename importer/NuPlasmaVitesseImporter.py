# A python-based importer for iolite 4 starts with some metadata
#/ Type: Importer
#/ Name: NuPlasma Vitesse importer
#/ Authors: Joe Petrus and Bence Paul
#/ Description: An importer for loading Nu Plasma Vitesse files
#/ References: None
#/ Version: 0.1
#/ Contact: support@iolite-software.com

import time
import numpy as np
import pandas as pd
import re
from datetime import datetime
from iolite.QtCore import QDateTime


def accepted_files():
    return ".run"


def correct_format():
    IoLog.debug("correct_format for Vitesse file importer called on file = %s"%(importer.fileName))

    if importer.fileName.endswith('csv'):
        nrf_regex = r"Cycle time \(ms\),x"
        with open(importer.fileName) as f:
            match = re.search(nrf_regex, f.read(), re.MULTILINE)
            if match:
                return True
            else:
                IoLog.debug(".csv file did not match the Vitesse format.")

    return False


def import_data():

    importer.message('Starting Vitesse import')
    importer.progress(0)

    mdf = None
    log_data = []

    with open(importer.fileName, 'r') as f:
        header = [next(f).strip() for x in range(15)]

    importer.message('reading metadata')
    importer.progress(5)

    metadata = {}
    for line in header:
        try:
            key, value = re.search(r'(.+):\s?,\s?(.+)', line).groups()
            metadata[key] = value
        except:
            pass

    startt = datetime.fromisoformat(metadata['Timestamp']).timestamp()
    df = pd.read_csv(importer.fileName, skiprows=len(metadata))
    df = df.rename(columns=lambda x: x.strip())
    df['Cycle time (ms)'] += startt*1000.
    df = df.append(pd.Series([np.nan]*len(df.columns)), ignore_index=True)
    df.iloc[-1]['Cycle time (ms)'] = df.iloc[-2]['Cycle time (ms)'] + 0.1

    if mdf is None:
        mdf = df
    else:
        mdf = pd.concat([mdf, df])

    startDT = datetime.fromtimestamp(df.iloc[0]['Cycle time (ms)']/1000.)
    endDT = datetime.fromtimestamp(df.iloc[-1]['Cycle time (ms)']/1000.)
    data.createFileSampleMetadata(metadata['Laser line name'], startDT, endDT, importer.fileName)
    log_data.append(
        {
            'startTime': QDateTime.fromMSecsSinceEpoch(startDT.timestamp()*1000.),
            'endTime': QDateTime.fromMSecsSinceEpoch(endDT.timestamp()*1000.),
            'name': metadata['Laser line name'],
            'repRate': np.nan,
            'width': float(metadata['Spot size']),
            'height': float(metadata['Spot size']),
            'shape': 'square',
            'scanSpeed': np.nan,
            'angle': 0.,
            'startX': df.iloc[0]['x [um]'],
            'startY': df.iloc[0]['y [um]'],
            'endX': df.iloc[-2]['x [um]'],
            'endY': df.iloc[-2]['y [um]']
        }
    )

    mdf = mdf.sort_values(by=['Cycle time (ms)'])

    importer.message('Adding channel data')
    importer.progress(30)

    cps_cols = [c for c in mdf.columns if 'cps' or 'ppm' in str(c)]
    channel_list = []  # List of channels loaded to be reported in UI
    for col in cps_cols:
        try:
            m, el = re.search(r'(\d+)(\w+) cps', col).groups()
            data.addDataToInput(f'{el}{m}', mdf['Cycle time (ms)']/1000., mdf[col], {'machineName': 'Vitesse'})
            channel_list.append(f'{el}{m}')
        except:
            try:
                m, el = re.search(r'(\d+)(\w+) ppm', col).groups()
                data.addDataToInput(f'{el}{m}', mdf['Cycle time (ms)']/1000., mdf[col], {'machineName': 'Vitesse'})
                channel_list.append(f'{el}{m}')
            except:
                pass

    importer.message('Adding file objects')
    importer.progress(90)

    data.createLaserLog(importer.fileName, log_data)

    data.createImportedFileMetadata(
        datetime.fromtimestamp(mdf.iloc[0]['Cycle time (ms)']/1000.),
        datetime.fromtimestamp(mdf.iloc[-1]['Cycle time (ms)']/1000.),
        importer.fileName,
        datetime.now(),
        len(mdf),
        channel_list
    )

    # Now calculate Total Beam:
    data.calculateTotalBeam()

    importer.message('Finished')
    importer.progress(100)
    importer.finished()
