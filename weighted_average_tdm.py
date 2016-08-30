#!/usr/bin/env python
"""
  Libraries and routines for calculation of the
  weighted averages of the TDM orientation
  from a histogram

  As most often in my coding to be simple and fast,
  I assume quite a bunch of stuff here, e.g.:
    Filenaming is of particular type
    the format is specific,
    data range and the numbers fit into some criteria
     (or just different data was not tested)
"""
# -*- coding: utf-8 -*-

#import pytram
import numpy as np
import math
import matplotlib.pyplot as plt
#import cPickle

# constants
kT = 2.5775  # kJ/mol*310K
degtorad = math.pi/180.0


def weighted_avg_and_std(values, weights):
    """
    Parameters
    ----------
    values : np.array
        data to be averaged
    weights : np.array
        the weight factors corresponding to the data
        -- both Numpy ndarrays with the same shape.

    Returns
    -------
    Returns the weighted average and its variance as a tuple
    """
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise
    return (average, math.sqrt(variance))



#%%

if __name__ == '__main__':
    dfnm = "tdm_dTRAM-FEP.dat"

    print("Loading-in the histogram")
    fep = np.loadtxt(dfnm)
    # weight factors, assuming fep in kT
    wfs = np.exp(-fep[:,1])
    x = fep[:,0]

#%%
    print "average weighted value and std: \n", weighted_avg_and_std(values=x, weights=wfs)
    cos2 = np.cos(x*degtorad)**2
    cos2avgstd = weighted_avg_and_std(values=cos2, weights=wfs)
    print "average weighted cos2 and std: \n", cos2avgstd
    sin2 = np.sin(x*degtorad)**2
    sin2avgstd = weighted_avg_and_std(values=sin2, weights=wfs)
    print "average weighted sin2 and std: \n", sin2avgstd
    print "dichroic ratio r: \n", sin2avgstd[0]/cos2avgstd[0]

    hist  = np.histogram(x, weights=wfs, bins=90, range=(0.0, 90.0), density=True)

    plt.xlabel("tdm orientation [deg]")
    plt.ylabel("probability density")
    binwidth = hist[1][1] - hist[1][0]
    plt.plot(hist[1][:-1]+binwidth*0.5, hist[0], lw=2.0, color='black')
    plt.savefig("tdm_WHAM-distrib.png", papertype="letter", dpi=300)
