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
    y = wfs = np.exp(-fep[:,1])
    x = fep[:,0]
    fep_err_estim = 0.9999999999  # error-bar estimate of FEP profile in kT

#%%
    ### 1PPM  ###
    # intenzities projected onto z & y axes with 1PPM (I~|dip*E|^2)
    rho_z = np.trapz(y*np.power(np.cos(np.deg2rad(x)), 2), x) *2.0 *np.pi
    rho_y = np.trapz(y*np.power(np.sin(np.deg2rad(x)), 2), x)      *np.pi
    print "\n1PPM\ndichr. ratio r = Rho_z / Rho_y = {} / {} = {} \nlog_2(r) = {}".format(rho_z, rho_y, rho_z/rho_y, math.log(rho_z/rho_y, 2))
    # intenzities projected onto z & y axes with 1PPM (I~|dip*E|^2)
    rho_z = np.trapz((y*np.exp(-fep_err_estim*np.cos(np.deg2rad(2.0*x))))*np.power(np.cos(np.deg2rad(x)), 2), x) *2.0 *np.pi
    rho_y = np.trapz((y*np.exp(-fep_err_estim*np.cos(np.deg2rad(2.0*x))))*np.power(np.sin(np.deg2rad(x)), 2), x)      *np.pi
    print "\nError estimate using FEP+{}*cos(2x)\n1PPM\ndichr. ratio r = Rho_z / Rho_y = {} / {} = {} \nlog_2(r) = {}".format(fep_err_estim, rho_z, rho_y, rho_z/rho_y, math.log(rho_z/rho_y, 2))
    # intenzities projected onto z & y axes with 1PPM (I~|dip*E|^2)
    rho_z = np.trapz((y*np.exp(+fep_err_estim*np.cos(np.deg2rad(2.0*x))))*np.power(np.cos(np.deg2rad(x)), 2), x) *2.0 *np.pi
    rho_y = np.trapz((y*np.exp(+fep_err_estim*np.cos(np.deg2rad(2.0*x))))*np.power(np.sin(np.deg2rad(x)), 2), x)      *np.pi
    print "\nError estimate using FEP-{}*cos(2x)\n1PPM\ndichr. ratio r = Rho_z / Rho_y = {} / {} = {} \nlog_2(r) = {}".format(fep_err_estim, rho_z, rho_y, rho_z/rho_y, math.log(rho_z/rho_y, 2))


    ### 2PPM  ###
    # intenzities projected onto z & y axes with 2PPM (I~|dip*E|^4)
    rho_z = np.trapz(y*np.power(np.cos(np.deg2rad(x)), 4), x) * 2.0     *np.pi
    rho_y = np.trapz(y*np.power(np.sin(np.deg2rad(x)), 4), x) * 3.0/4.0 *np.pi
    print "\n2PPM\ndichr. ratio r = Rho_z / Rho_y = {} / {} = {} \nlog_2(r) = {}".format(rho_z, rho_y, rho_z/rho_y, math.log(rho_z/rho_y, 2))

    # intenzities projected onto z & y axes with 2PPM (I~|dip*E|^4)
    rho_z = np.trapz((y*np.exp(-fep_err_estim*np.cos(np.deg2rad(2.0*x))))*np.power(np.cos(np.deg2rad(x)), 4), x) * 2.0     *np.pi
    rho_y = np.trapz((y*np.exp(-fep_err_estim*np.cos(np.deg2rad(2.0*x))))*np.power(np.sin(np.deg2rad(x)), 4), x) * 3.0/4.0 *np.pi
    print "\nError estimate using FEP-{}*cos(2x)\n2PPM\ndichr. ratio r = Rho_z / Rho_y = {} / {} = {} \nlog_2(r) = {}".format(fep_err_estim, rho_z, rho_y, rho_z/rho_y, math.log(rho_z/rho_y, 2))
    # intenzities projected onto z & y axes with 2PPM (I~|dip*E|^4)
    rho_z = np.trapz((y*np.exp(+fep_err_estim*np.cos(np.deg2rad(2.0*x))))*np.power(np.cos(np.deg2rad(x)), 4), x) * 2.0     *np.pi
    rho_y = np.trapz((y*np.exp(+fep_err_estim*np.cos(np.deg2rad(2.0*x))))*np.power(np.sin(np.deg2rad(x)), 4), x) * 3.0/4.0 *np.pi
    print "\nError estimate using FEP-{}*cos(2x)\n2PPM\ndichr. ratio r = Rho_z / Rho_y = {} / {} = {} \nlog_2(r) = {}".format(fep_err_estim, rho_z, rho_y, rho_z/rho_y, math.log(rho_z/rho_y, 2))


    hist  = np.histogram(x, weights=wfs, bins=90, range=(0.0, 90.0), density=True)

    plt.xlabel("tdm orientation [deg]")
    plt.ylabel("probability density")
    binwidth = hist[1][1] - hist[1][0]
    plt.plot(hist[1][:-1]+binwidth, hist[0], lw=2.0, color='black')
    plt.savefig("tdm_WHAM-distrib_a.png", papertype="letter", dpi=300)
