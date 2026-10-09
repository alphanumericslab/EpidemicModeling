"""Offline data readers and explicit preprocessing for examples.

Author: Reza Sameni, Emory University.
"""
import numpy as np
import pandas as pd
from scipy.signal import lfilter


def read_covid19_data(confirmed_datafile, death_datafile, recovered_datafile, region_list, min_cases):
    """Aggregate Johns Hopkins wide time-series tables by country substring.

    Inputs are CSV paths with four metadata columns followed by aligned
    date columns. Returns (total, infected, recovered, deceased, first_index,
    threshold_index, num_days). Case arrays have shape (regions, days).
    Indices are one-based for MATLAB parity, with 0 when never reached.
    """
    tables = [pd.read_csv(p) for p in (confirmed_datafile,death_datafile,recovered_datafile)]
    dates = list(tables[0].columns[4:])
    if any(list(t.columns[4:]) != dates for t in tables): raise ValueError("date columns must align")
    results = []
    for table in tables:
        country = table.iloc[:,1].fillna('').astype(str)
        results.append(np.array([table.loc[country.str.contains(str(region),regex=False),dates].sum().to_numpy(float) for region in region_list]))
    total,deceased,recovered = results
    first = np.array([np.flatnonzero(row>0)[0]+1 if np.any(row>0) else 0 for row in total])
    threshold = np.array([np.flatnonzero(row>=min_cases)[0]+1 if np.any(row>=min_cases) else 0 for row in total])
    return total,total-deceased-recovered,recovered,deceased,first,threshold,len(dates)


def read_oxford_data(path):
    """Read an Oxford/XPRIZE CSV, normalize region blanks, and sort daily rows.

    Date accepts YYYYMMDD or ISO strings. CountryName, Date and ConfirmedCases
    are required. Duplicate country/region/date observations are rejected
    because they make daily fitting ambiguous. No online download is used.
    """
    frame = pd.read_csv(path)
    for column in ('CountryName','Date','ConfirmedCases'):
        if column not in frame: raise ValueError(f"missing column {column}")
    if 'RegionName' not in frame: frame['RegionName'] = ''
    frame['RegionName'] = frame.RegionName.fillna('')
    date = frame.Date.astype(str)
    frame['Date'] = pd.to_datetime(date,format='%Y%m%d' if date.str.fullmatch(r'\d{8}').all() else '%Y-%m-%d')
    if frame.duplicated(['CountryName','RegionName','Date']).any():
        raise ValueError("duplicate country/region/date rows")
    return frame.sort_values(['CountryName','RegionName','Date']).reset_index(drop=True)


def prepare_cases(cumulative, window=7):
    """Return cleaned daily cases and their causal zero-padded moving average.

    Daily cases begin with zero, negative revisions are clipped to zero,
    and nonfinite differences become zero. This preprocessing uses only
    present and past observations and is shared by both updated languages.
    """
    if int(window) != window or window < 1: raise ValueError("window must be a positive integer")
    daily = np.maximum(0,np.nan_to_num(np.diff(np.asarray(cumulative,float),prepend=np.nan),nan=0,posinf=0,neginf=0))
    return daily,lfilter(np.ones(window)/window,[1],daily)


def read_geo_table(path):
    """Read a geographic CSV and normalize missing RegionName values to empty strings.

    CountryName is required; absent RegionName is added. Original column
    names and non-geographic fields are preserved. Matches the MATLAB helper.
    """
    data = pd.read_csv(path)
    if 'CountryName' not in data: raise ValueError('missing CountryName')
    if 'RegionName' not in data: data['RegionName'] = ''
    data['RegionName'] = data.RegionName.fillna('')
    return data
