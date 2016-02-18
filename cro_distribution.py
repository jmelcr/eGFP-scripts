#!/usr/bin/env python
"""
  Libraries and routines for calculation of the
  re-weighted distribution of the TDM orientation
  in a WHAM-like fashion.
   This assumes that each individual simulation (window)
  was sampled from its global equilibrium.
   This is usually a good assumption for Umbrella sampling simulations.

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
import cPickle

# constants
kT = 2.5775  # kJ/mol*310K


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


def reorder_pytram_trajs_from_PITCH(trajs):
    """
    PITCH files from Plumed contain data in
    a different order than assumed in pytram.Reader

    This is a simple routine, that changes this glitch and
    adds 'time' key to the traj-dictionary

    Keys
    ----
    time : float
        time
    m    : integer
        Markov state
    b    : float
        bias
    t    : integer
        thermodynamic state

    all above quantities are numpy arrays of the defined type
    """
    for i,traj in enumerate(trajs.trajs):
        if traj.has_key('time'):
            print("this traj was already processed")
        else:
            traj['time'] = traj.pop('m')
            traj['m'] = traj.pop('t')
            traj['t'] = np.ones(traj['m'].shape, dtype=int) * i
    return trajs


def read_trajs(dfnm, NWINS=None):
    """
    Reads in Plumed's PitchTdmBias files

    Parameters
    ----------
    dfnm : string
        default filenaming of the PitchTdmBias files
        these files typically contain several columns
        with values time, pitch, TDM-orientation, bias
    NWINS : int
        number of simulation windows
        this number is added at the end of the filename
         If not given (is None), no number is added at the end
        and only one filed is read

    Returns
    -------
    trajs : np.array
        array of everything in the files provided
    """
    files = []
    for i in range(NWINS):
        files.append(dfnm+"."+str(i))

    # load the PITCH.* files using pytram's reader
    print("Loading-in the PitchTdm trajectories...")
    trajs = []
    for filename in files:
        trajs.append(np.loadtxt(filename))

    return trajs


def get_weight_factors(trajs, winbiases):
    """
    Gives weight factors as read from Plumed's Pitch files

    Parameters
    ----------
    trajs : list of arrays
        obtained from read_trajs
    winbiases : np.array
        of windows biasing weight factors (in kT units)
        read-in from previous dTRAM or WHAM estimation of FEP

    Returns
    -------
    wfs : list of np.arrays
        list of arrays of weight factors (floats)
        with dTRAM or WHAM estimated window weights
    """
    print("Obtaining weight factors and reshaping biases ...")
    wfs = []
    for i, traj in enumerate(trajs):
        # the last column should contain biases in kJ/mol - divide by kT,
        # they are with +sign as they are used for reweighting
        # back to Boltzmann distribution.
        # window weights are already in kT units,
        # they are with -sign as they give
        # the Boltzmann-weights for individual windows
        wfs.append(np.exp(traj[:, -1]/kT - winbiases[i]))
    return wfs


    def weight_switch(vals, thres):
        """
        """
        j=0
        for i in xrange(vals.shape[0]):
            if vals[i] > thres:
                delta = abs(vals[i]-thres)
                vals[i] = thres + delta*math.exp(-delta/thres)
                j += 1
        print j, "values out of", vals.shape[0], "(", j/vals.shape[0]*100.0, "%) were above threshold", thres, ".\n Their weights were reduced exponentially. "
        return vals



#%%

if __name__ == '__main__':
    dfnm = "PitchTdmBias.traj"        # Plumed PITCH filename convention
    max_weight_thres = math.exp(3.0)  # threshold for maximal weight factor (~3kT here)
    # weight factors of the windows
    print("Loading-in the simulation windows' weight factors...")
    with open("pitch_dTRAM-FEP_windows.pickle","r") as f: winbiases = cPickle.load(f)
    NWINS = winbiases.shape[0]
    winbiases -= winbiases.min()

    trajs = read_trajs(dfnm, NWINS)

    wfs = get_weight_factors(trajs, winbiases)
#%%

    wfs_flat = wfs[0]
    for wf in wfs[1:]:
        wfs_flat = np.append(wfs_flat, wf)

    tdms_flat = trajs[0][:,-2]
    for traj in trajs[1:]:
        tdms_flat = np.append(tdms_flat, traj[:,-2])


    # get rid of over-weighted values (that might destroy my statistics)
    wfs_flat = weight_switch(wfs_flat, max_weight_thres)


#%%
    print "average weighted value and std: \n", weighted_avg_and_std(values=tdms_flat, weights=wfs_flat)

    hist  = np.histogram(tdms_flat, weights=wfs_flat, bins=90, range=(0.0, 90.0), density=True)
    plt.plot(hist[1][:-1], hist[0])
    plt.show()
    #plt.plot(tdms_flat)
    #plt.show()

    wfs_flat.max()

