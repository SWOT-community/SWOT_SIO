#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This file contains routines for manipulating SWOT data.

Author: David Sandwell and Yao Yu
Date: First version: Oct. 6, 2023

Dependencies:
numpy 
xarray 
os
cartopy.crs 
matplotlib.pyplot 
pandas 
pygmt
scipy
"""

import numpy as np
import xarray as xr
import os
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import pandas as pd
import pygmt
import scipy.signal

def swath_unit_vectors(lon, lat):
    """
    compute the unit direction vectors for each cell in a swath
    
    Parameters
    ----------
    lat : array_like
        2D array of latitude values with shape (na, nc).
    lon : array_like
        2D array of longitude values with shape (na, nc).

    Returns
    -------
    uea : array_like
	2D array of the east component of the init vector in the along-track direction
    una : array_like
	2D array of the north component of the init vector in the along-track direction
    uec : array_like
	2D array of the east component of the init vector in the cross-track
    unc : array_like
	2D array of the north component of the init vector in the cross-track
    """
    d2r = np.pi/180.0
    cost = np.cos(lat*d2r)
    dna, dnc = np.gradient(lat,1,1)
    dea, dec = np.gradient(lon,1,1)*cost
    dsa = (dea*dea + dna*dna)**0.5
    uea = dea/dsa
    una = dna/dsa
    dsc = (dec*dec + dnc*dnc)**0.5
    uec = dec/dsc
    unc = dnc/dsc

    return uea, una, uec, unc
