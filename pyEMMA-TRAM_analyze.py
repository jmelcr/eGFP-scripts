#!/usr/bin/env python
"""
 TRAM/dTRAM analysis of a very SPECIFIC REUS/MSM simulations of xxeGFPs
 Beware: this script might still be a bit too specific to this task.

 ------------------------------------------------------------
 Made by Joe,  Last edit 2016/08/31
------------------------------------------------------------
 input: Plumed's pitch files,
        MetaData in Grossfield's-WHAM format
 output: FEP (pickled, textfile and plot)

--------------------------------------------------------
"""


# coding: utf-8

# import basic classes and routines developed for dTRAM analysis of eGFPs:
from eGFP_FEP_dTRAM_REUS_analysis import *
# import other public packages:
import pyemma
import msmtools
import numpy as np
import scipy, math
import matplotlib.pyplot as plt
import cPickle
from optparse import OptionParser
#import IPython

try:
    import read_xvg_calc_mean as rxvg
except:
    print "Couldn't read XVG-reading library, can't process Gromacs pullx.xvg files!"


#%%

if __name__ == "__main__":
    # help message is automatically provided
    # type=string, action=store is default
    parser = OptionParser()
    parser.add_option('-m', '--metafile', dest='wham_metadata_fname', help='wham metadata file name', default="wham_meta.data")
    parser.add_option('-w', '--whamfep',  dest='wham_fep_fname', help='wham FEP file name', default="wham_fep_histo.dat")
    parser.add_option('-o', '--outfname', dest='out_file_name',  help='output file name (prefix)', default="pitch")
    parser.add_option('-f', '--trajform', dest='traj_file_format', help='file format of the trajectories (def. plain, xvg)', default="plain")
    parser.add_option('-b', '--binw', dest='bin_width', help='bin width', default=1.0, type=float)
    parser.add_option('-i', '--niter', dest='nruns_max', help='no. iter*1k (=1run)', default=100, type=int)
    parser.add_option('-l', '--lag', dest='lag_time', help='lag time (in units of frame-rec-rate) for dTRAM c_k_i_j kinetic mtx', default=1000, type=int)
    parser.add_option('-t', '--ftol', dest='f_toler_dtram', help='convergence tolerance for dTRAM iterative sc-procedure', default=1.0E-10, type=float)
    parser.add_option('-e', '--errestim', dest='err_estim', help='number of bootstrapping iterations (randomized trajs) for error estimation', default=10, type=int)
    opts, args = parser.parse_args()

    niter_max = opts.nruns_max*1000
    #read-in the simulations' metadata in the Grossfield's WHAM format
    # and the trajectories in the PLUMED's format
    sim_data = SimData(opts.wham_metadata_fname, opts.bin_width, traj_file_format=opts.traj_file_format)

    # just an alias, gridpoints are saved as an attribute in the SimData obj
    gridpoints = sim_data.gridpoints

#%%

    #IPython.embed()
    #us_centres = [ sim['x0'] for sim in sim_data.sims[:] ]
    #us_force_constants = [ sim['kappa'] for sim in sim_data.sims[:] ]
    #md_trajs  = [ sim['x'] for sim in sim_data.trajs[:] ]
    md_ttrajs = [ np.ravel(sim['t']) for sim in sim_data.trajs[:] ]
    md_dtrajs = [ np.ravel(sim['m']) for sim in sim_data.trajs[:] ]
    if 'b' in sim_data.trajs[0].keys():
        md_bias   = [ sim['b'] for sim in sim_data.trajs[:] ]
    else:
        md_bias   = [ np.zeros(shape=sim['m'].shape) for sim in sim_data.trajs[:] ]
    #us_dtrajs = pyemma.coordinates.assign_to_centers(data=us_trajs, centers=us_centres) # biased discrete trajectories
    #md_dtrajs = pyemma.coordinates.assign_to_centers(data=md_trajs, centers=centers) # unbiased discrete trajectories

    # TRAM approach
    print " \r\n TRAM running ..."
    tram = pyemma.thermo.tram(
    ttrajs=md_ttrajs, dtrajs=md_dtrajs,
    bias=md_bias,
    lag=opts.lag_time,
    maxiter=niter_max, maxerr=opts.f_toler_dtram, save_convergence_info=10,
    init='mbar', init_maxiter=2000, init_maxerr=opts.f_toler_dtram)
    #tram_obj = tram[4]  # same abbreviation as in Jupyter notebook - prbably works with the estimate_... thing
    #estimate_umbrella_sampling(
    #us_trajs, us_dtrajs, us_centres, us_force_constants,

#%%
    # MBAR approach
    print " \r\n MBAR running ..."
    mbar = pyemma.thermo.mbar(
    ttrajs=md_ttrajs, dtrajs=md_dtrajs,
    bias=md_bias,
    maxiter=niter_max, maxerr=opts.f_toler_dtram, save_convergence_info=10)

#%%

    # dTRAM approach
    print " \r\n dTRAM running ..."
    dtramdata = TRAMData( sim_data.trajs,
    b_K_i=sim_data.gen_harmonic_bias_mtx(gridpoints=sim_data.gridpoints) )

    try:
        dtram_obj = dtram(dtramdata, lag=opts.lag_time, maxiter=1,
                          ftol=opts.f_toler_dtram, verbose=False )
    except pytram.ExpressionError:
        raise pytram.ExpressionError, "Input was faulty, ending..."
    except pytram.NotConvergedWarning:
        pass


#%%

    converged = False
    nruns=0

    while not converged:
       try:
           print "Running another 1k steps ..."
           nruns+=1
           if nruns > opts.nruns_max:
               print """BEWARE: Did not converge within 50k iterations, this is weird, ending. \nBEWARE: generated results are based on unconverged data! """
               break
           dtram_obj.sc_iteration(maxiter=1000)
       except:
           print "Not converged in "+str(nruns)+"k steps. ",
       else:
           print "Converged within "+str(nruns)+"k steps."
           converged = True



#%%

    print "Bootstrapping trajectories"
    # Trajectory Bootstrap error estimate
    md_dtrajs_bootstrapped = []
    for i in range(opts.err_estim):
        md_dtrajs_bootstrapped.append(msmtools.estimation.bootstrap_trajectories(md_dtrajs, 0.5))  #opts.lag_time))

#%%

    print "Estimating errors ..."
    tram_bootstrapped_objs = []
    for i,dtraj in enumerate(md_dtrajs_bootstrapped):
        # should be fine at least for pure MSM (not sure about biased simulations and sims that change thermo state)
        print " \r\n Bootstrapping traj no.", i
        ttraj = [ np.zeros(shape=part.shape, dtype=int)   for part in dtraj ]
        bias  = [ np.zeros(shape=part.shape+(1,), dtype=float) for part in dtraj ]  # the "+(1,) is a workaround that made it working, but I don't uderstand it completely
        tram_obj = pyemma.thermo.tram(  #mbar much faster than .. tram(
        dtrajs=dtraj,
        bias=bias, ttrajs=ttraj,
        lag=opts.lag_time,
        maxiter=niter_max, maxerr=opts.f_toler_dtram)
        tram_bootstrapped_objs.append(tram_obj)

#%%
    # get the standard error
    mean = tram.free_energies.mean()
    std_sq_sum = np.zeros(shape=tram.free_energies.shape)
    for tram_obj in tram_bootstrapped_objs:
        diff = tram_obj.free_energies.min() - tram.free_energies.min()  # shift to minimum  --or--  shift to mean: tram_obj.free_energies.mean() - mean
        active_set = tram_obj.active_set
        std_sq = np.square( tram_obj.free_energies[active_set-active_set.min()] - tram.free_energies[active_set] - diff )
        std_sq_sum[active_set] += std_sq[active_set-active_set.min()]
    # make average and sqrt to get standard error
    tram_std_err = np.sqrt(std_sq_sum/len(tram_bootstrapped_objs))

#%%
    # quickly plot the free energies from the original and bootstrap-randomized trajs
#    plt.plot(tram.active_set, tram.free_energies, lw=2, color='black')
#    for tram_obj in tram_bootstrapped_objs:
#        plt.plot(tram_obj.active_set, tram_obj.free_energies)
#    plt.show()

#%%

    print "Shifting and saving the profile into a pickle-object and text file..."

    # "aliases" for plotting and saving
    #dTRAM data
    fep_dtram = dtram_obj.f_i
    fep_dtram_shift = np.copy(fep_dtram)
    fep_dtram_shift -= fep_dtram_shift.min()
    #TRAM data
    fep_tram = tram.free_energies
    fep_tram_shift = np.copy(fep_tram)
    fep_tram_shift -= fep_tram_shift.min()
    #mbar data
    fep_mbar = mbar.free_energies
    fep_mbar_shift = np.copy(fep_mbar)
    fep_mbar_shift -= fep_mbar_shift.min()

    with open(opts.out_file_name+"_dTRAM-FEP_windows.pickle","w") as f : cPickle.dump(dtram_obj.f_K,f)
    with open(opts.out_file_name+"_dTRAM-FEP.pickle","w") as f : cPickle.dump(fep_dtram_shift,f)
    with open(opts.out_file_name+"_dTRAM-FEP.dat","w") as f:
        for i,c in enumerate(gridpoints):
           line = str(c)+"   "+str(fep_dtram_shift[i])+"\n"
           f.write( line )
    with open(opts.out_file_name+"_TRAM-FEP.dat","w") as f:
        for i,c in enumerate(gridpoints):
           line = str(c)+"   "+str(fep_tram_shift[i])+"\n"
           f.write( line )



#%%

    print "Plotting free energies ... "

    plt.plot( gridpoints, fep_dtram_shift, '-', color='black', lw=2.0, label="dTRAM" )
    plt.plot( gridpoints , fep_mbar_shift, '--', color='grey', lw=2.0, label="MBAR pyEMMA", alpha=0.5)
    plt.plot( gridpoints , fep_tram_shift, '--', color='red',  lw=1.0, label="TRAM pyEMMA")

    fep_tram_err1 =  tram_std_err + fep_tram_shift
    fep_tram_err2 = -tram_std_err + fep_tram_shift
    plt.fill_between( gridpoints, fep_tram_err1, fep_tram_err2, color='orange', alpha=0.25, label="TRAM pyEMMA bootstrapped err. estim.")

    plt.legend( loc='upper center', fontsize=10 )
    plt.xlabel( r"$pitch$ [deg]", fontsize=12 )
    plt.ylabel( r"$G$ [kT]", fontsize=12 )

    plt.savefig(opts.out_file_name+"_pyEMMA-TRAM_dTRAM.png", papertype="letter", dpi=300)

#%%
    # plotting dTRAM - Jacobi corrected (assuming Pitch-polar coordinates)
    kT = k_b*sim_data.sims[0]['temperature']
    plt.plot( gridpoints, fep_dtram_shift + kT*np.log(np.sin(gridpoints/180.0*math.pi)), '-', color='blue', lw=2.0, label="dTRAM Jacobi corrected" )
    plt.plot( gridpoints,  fep_tram_shift + kT*np.log(np.sin(gridpoints/180.0*math.pi)), '--', color='red' , lw=1.0, label="pyEMMA TRAM Jacobi corrected" )

    plt.legend( loc='upper left', fontsize=10 )
    plt.xlabel( r"$pitch$ [deg]", fontsize=12 )
    plt.ylabel( r"$G$ [kT]", fontsize=12 )

    plt.savefig(opts.out_file_name+"_pyEMMA-TRAM_dTRAM_JacobiCorr.png", papertype="letter", dpi=300)
