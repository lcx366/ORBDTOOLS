# Welcome to the ORBDTOOLS package

[![PyPI version shields.io](https://img.shields.io/pypi/v/orbdtools.svg)](https://pypi.python.org/pypi/orbdtools/) [![PyPI pyversions](https://img.shields.io/pypi/pyversions/orbdtools.svg)](https://pypi.python.org/pypi/orbdtools/) [![PyPI status](https://img.shields.io/pypi/status/orbdtools.svg)](https://pypi.python.org/pypi/orbdtools/) [![GitHub contributors](https://img.shields.io/github/contributors/lcx366/ORBDTOOLS.svg)](https://GitHub.com/lcx366/ORBDTOOLS/graphs/contributors/) [![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/lcx366/ORBDTOOLS/graphs/commit-activity) [![GitHub license](https://img.shields.io/github/license/lcx366/ORBDTOOLS.svg)](https://github.com/lcx366/ORBDTOOLS/blob/master/LICENSE) [![Documentation Status](https://readthedocs.org/projects/orbdtools/badge/?version=latest)](http://orbdtools.readthedocs.io/?badge=latest) [![Build Status](https://travis-ci.org/lcx366/orbdtools.svg?branch=master)](https://travis-ci.org/lcx366/orbdtools)

This package is on its way to become an archive of scientific routines for data processing related to Arc Matching, Arc Associating, Initial Orbit Determination(IOD), Cataloging OD, and Precise OD. Currently, the package only implements a small number of functional modules,  and subsequent new modules will be added and updated one after another. So far,  operations on Arc Matching include:

1. Matching of observation arcs based on optical angle measurement data to space objects in TLE file; 
2. Matching of observation arcs based on radar measurement data(range+angle) to space objects in TLE file;
3. Data processing related to TLE file:
   - Reading and parsing of TLE file
   - Calculation of the mean orbital elements(only long-term items are considered) in a certain epoch
   - Orbital Propagation

## How to Install

On Linux, macOS and Windows architectures, the binary wheels can be installed using `pip` by executing one of the following commands:

```
pip install orbdtools
pip install orbdtools --upgrade # to upgrade a pre-existing installation
```

## How to use

### Data processing related to TLE files

#### Download TLE files from [SPACETRACK](https://www.space-track.org)

```python
>>> from orbdtools import TLE
>>> tle_file = TLE.download([52132,51454,37637,26758,44691])
>>> # tle_file = TLE.download('satno.txt')
```

A sample of  the *satno.txt* file is as follows:

```
#satno
12345
23469
45678
```

#### Read and parse a TLE file

```python
>>> import numpy as np
>>>
>>> tle_file = 'test/tle_20220524.txt'
>>> epoch_obs = '2022-05-24T08:38:34.000Z' # Approximate epoch of the observed arc, which is optional
>>> tle = TLE.from_file(tle_file,epoch_obs) # Only TLEs within a week before and after the observation epoch are considered
>>> # tle = TLE.from_file(tle_file) # All TLEs are considered
>>> print(tle.range_epoch) # Epoch coverage for TLE files in format of [min, median, max]
>>> print(tle.df)
```

    ['2022-05-18T23:07:15.444Z', '2022-05-23T22:00:02.000Z', '2022-05-24T06:32:39.874Z']

<p align="middle">
  <img src="readme_figs/df.png" width="600" />
</p>

#### Calculate the mean orbital elements at a certain epoch

```python
>>> tle_epoch = tle.atEpoch(epoch_obs)
>>> print(tle_epoch.df)
```

<p align="middle">
  <img src="readme_figs/df_epoch.png" width="600" />
</p>

#### Orbital Propagation

##### Calculate the cartesian coordinates of space objects in GCRF(Geocentric Celetial Reference Frame) over a period of time

```python
>>> t_list = ['2022-05-23T16:22:49.408Z', '2022-05-23T18:48:34.488Z',
           '2022-05-23T18:09:35.640Z', '2022-05-23T20:54:20.228Z',
           '2022-05-23T20:29:03.621Z', '2022-05-23T23:24:11.831Z',
           '2022-05-24T00:38:08.803Z', '2022-05-24T02:33:32.466Z',
           '2022-05-23T21:10:37.703Z', '2022-05-23T17:48:24.865Z']
>>> xyz_gcrf = tle.predict(t_list)
>>> # sats_index = [10,20,30,40] # Index of space objects
>>> # xyz_gcrf = tle.predict(t_list,sats_index)
>>> print(xyz_gcrf.shape)
```

    (4904, 10, 3)

### Arc Matching

Load the observation file and extract the necessary data.

```python
>>> import numpy as np
>>> 
>>> obs_data = np.loadtxt('test/obs.dat',dtype=str,skiprows=1) # Load the observation file
>>> # extract the necessary data
>>> t = obs_data[:,0] # Obsevation time in UTC
>>> xyz_site = obs_data[:,1:4].astype(float) # Cartesian coordinates of the site in GCRF, [km]
>>> radec = obs_data[:,4:6].astype(float) # Ra and Dec of space object, [deg]
>>> r = obs_data[:,6].astype(float) # Slant distance of the space object relative to the site, [km]
```

#### Match the observation arc based on optical angle measurement data to space objects in TLE file

```python
>>> from orbdtools import ArcObs
>>> arc_optical = ArcObs({'t':t,'radec':radec,'xyz_site':xyz_site}) # Load the necessary data
>>> 
>>> # Use the LOWESS(LOcally Weighted Scatterplot Smoothing) method to remove outliers, optional
>>> arc_optical.lowess_smooth() 
>>> arc_optical.arc_match(tle) # Match the observation arc to TLE
>>> 
>>> print(arc_optical.code_match)
>>> print(arc_optical.satnum)
>>> print(arc_optical.disp_match)
```

    1
    1616
    Target ID: 1616

Three types of matching results are summaried as follows

| Match Code | Object ID         | Solution Case      | Status  | What to Do Next    |
|:----------:|:-----------------:|:------------------:|:-------:|:------------------:|
| 1          | NORAD ID          | Unique solution    | Success |                    |
| 0          | None              | No solution        | Failure | increase threshold |
| -1         | list of NORAD IDs | Multiple solutions | Failure | decrease threshold |

#### Match the observation arc based on radar measurement data(range+angle) to space objects in TLE file

```python
>>> from orbdtools import ArcObs
>>> arc_radar = ArcObs({'t':t,'radec':radec,'r':r,'xyz_site':xyz_site}) # Added 'r' to optical angle measurement data 
>>> arc_radar.lowess_smooth()
>>> arc_radar.arc_match(tle)
>>> 
>>> print(arc_radar.code_match)
>>> print(arc_radar.satnum)
>>> print(arc_radar.disp_match)
```

    1
    1616
    Target ID: 1616

## Change log

- **0.0.1 — Jul 18, 2023**
  
  - The ***orbdtools*** package was released.

## Reference

- [Skyfield](https://rhodesmill.org/skyfield/)
- [sgp4](https://pypi.org/project/sgp4/)
- [lowess](https://pypi.org/project/loess/)