#!/usr/bin/env python
"""
  Calculate the weighted distribution of the TDM orientation.
  Several thing assumed:
     - PITCH recorded at 2x higher pace than XTC trajs
     - for some reason, the time stamps do not agree
     - everything fits into the memory at once
     - simulations can be of different length
  orientation: theta in degrees (pitch)
  high specificity to the problem assumed
  -- not easily ransferrable to other applications
"""
# -*- coding: utf-8 -*-

import pytram
import numpy as np
import math

import cPickle

# constants
kT = 2.5775  # kJ/mol*310K
dfnm = "PITCH."      # Plumed PITCH filename convention
dfname = "tdm_theta" # default filename for the TDM orientation files
dfext = ".dat"       #  -- and their extention



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


def get_weight_factors(dfnm, NWINS, winwfs, tdms, stride=2):
    """
    Gives weight factors as read from Plumed's Pitch files

    Parameters
    ----------
    dfnm : string
        default filenaming of the PITCH files
    NWINS : int
        number of simulation windows
    winwfs : np.array
        weight factors of the individual simulation windows (e.g. from dTRAM)
    stride : int
        multiplication factor for the difference
        in the pace of the PITCH records and TDM records

    Returns
    -------
    wfs : np.array
        flat array of weight factors (floats)
    """
    files = []
    for i in range(NWINS):
        files.append(dfnm+str(i))

    # load the PITCH.* files using pytram's reader
    print("Loading-in the Pitch trajectories...")
    trajs = pytram.Reader(files=files)
    print("Reordering the Pitch trajectories ... ")
    trajs = reorder_pytram_trajs_from_PITCH(trajs)
    print("Obtaining biases...")
    # working in units of kT here!
    biases = []
    for win in range(NWINS):
        # assign bias for each TDM record - prevents the case, where there are less TDMs than biases
        #for i,rec in enumerate(tdms[win]):  #this was the old way -1
        ## THIS IS A SOURCE OF ERRORS -- NRECORDS OF TDMS AND BIASES DON'T AGREE -- INVESTIGATE!!!
        for i in xrange(len(trajs.trajs[win]['b'])/stride):
            # BEWARE!!
            # PITCH is recorded at 2x pace than xtc-traj frames
            biases.append( np.asscalar( trajs.trajs[win]['b'][stride*i] )/kT + winwfs[win] )

    print("Obtaining weight factors and reshaping biases ...")
    wfs = np.exp( -np.asarray(biases) )   #.reshape( (NWINS,-1) )
    return wfs


def get_tdms(dfname, NWINS):
    """
    Gives weight factors as read from Plumed's Pitch files
    Specific file-format assumed/expected (not checked!)

    Parameters
    ----------
    dfname : string
        default filenaming of the TDM-orientation files
    NWINS : int
        number of simulation windows

    Returns
    -------
    tdms : np.array
        flat array of TDM orientations (floats)
    """
    print("Loading-in TDM orientations from .dat files (specific input assumed)")
    # load-in tdm-orientations
    tdms = []
    for i in range(NWINS):
        tdm = np.loadtxt(dfname+str(i)+dfext)
        tdms.append(tdm)

    # get the number of records 
    nrecords=0
    for tdm in tdms:
        nrecords += tdm.shape[0]

    # turn the  TDM-records into a flat np.array    
    tdms_flat_arr = np.zeros( (nrecords) )
    i=0
    for tdm in tdms:
        for j in tdm[:,1]:
            tdms_flat_arr[i] = j
            i+=1
    
    return tdms_flat_arr


if __name__ == '__main__':
    # weight factors of the windows
    print("Loading-in the simulation windows' weight factors...")
    with open("pitch_dTRAM-FEP_windows.pickle","r") as f: winwfs = cPickle.load(f)
    NWINS = winwfs.shape[0]

    tdms = get_tdms(dfname, NWINS)

    wfs = get_weight_factors(dfnm, NWINS, tdms, winwfs, stride=2)

    print( "average weighted value and std: \n", weighted_avg_and_std(values=tdms, weights=wfs) )    