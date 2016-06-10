#!/usr/bin/env python
"""
  Simple adiabacity check of dAFED simulations
  assuming only 1CV and a particular file format
"""

# coding: utf-8

import numpy as np
import math
from optparse import OptionParser
import IPython

k_b = 0.0083144621  #kJ/Mol*K

if __name__ == "__main__":
    # help message is automatically provided
    # type=string, action=store is default
    parser = OptionParser()
    parser.add_option('-f', '--file', dest='data_fname', help='data file name', default="ValMetDist_s-v-x-b")
    parser.add_option('-k', '--kappa', dest='kappa', help='kappa -- spring constant', default=666000, type=float)
    parser.add_option('-t', '--tau', dest='tau', help='tau -- time coupling constant', default=0.628, type=float)
    opts, args = parser.parse_args()

    data = np.loadtxt(opts.data_fname)
    time, s, v, x, b = np.hsplit(data,5)
    kappa = opts.kappa
    mass = kappa*(opts.tau/(2*math.pi))**2
    n_rec = data.shape[0]
    sim_len_time = float((time[-1] - time[0]))  # ps
    print "\nsimulation length = ", sim_len_time/1000.0, "ns \n"
    
    # M. Cuendet: http://lausanne.isb-sib.ch/~mcuendet/downloads/dafed0.8_manual.pdf
    # Configuration temperature
    # T_conf= kappa/k_b * <(x-S)**2>   

    print "T_C of parts of the trajectory:"
    traj_parts = 5
    T_C_sum = 0.0
    for i in range(traj_parts):
        lbound = (n_rec-1)/traj_parts * i
        ubound = (n_rec-1)/traj_parts * i+1
        T_C_part = kappa/k_b * np.mean((x-s).ravel()[lbound:ubound]**2)
        T_C_sum += T_C_part
        print i, "/", traj_parts, " ... ", T_C_part, "K"

    print "average .. ", T_C_sum/float(traj_parts), "K"

    T_C     = kappa/k_b * np.mean((x-s).ravel()**2)
    T_C_std = kappa/k_b *  np.std((x-s).ravel()**2)
    print "trj tot .. ", T_C, "+/-", T_C_std, "K \n"

    # W_S = int_0^t dt*kappa*(s-x)*v_s = int_0^s_final ds*kappa*(s-x)
    W_Sv = np.trapz( (kappa*(s-x)*v).ravel(), time.ravel() )
    print "W_Sv Transferred heat = ", W_Sv, "kJ/mol"
    W_Ss = np.trapz( kappa*(s-x).ravel(), s.ravel())
    print "W_Ss Transferred heat = ", W_Ss, "kJ/mol"
    print "P_S  Heat production  = ", float(1000*W_Ss / sim_len_time), "kJ/(mol*ns) \n"

    H_S     = np.mean(0.5*(mass*v**2 + kappa*(x-s)**2))
    H_S_std = np.std( 0.5*(mass*v**2 + kappa*(x-s)**2))
    print "H_S Energy of extended system = ", H_S, "+/-", H_S_std, "kJ/mol, rel. std = ", 100*H_S_std/H_S, "% \n"

    print "delta-H_S, end-begin, fraction of traj:"
    for traj_frac in (0.05, 0.1, 0.2, 0.3): # i.e. 10,20,30,40%
        n_rec_frac = int(traj_frac*n_rec)
        H_S_beg = np.mean(0.5*(mass*v**2 + kappa*(x-s)**2)[:n_rec_frac])
        H_S_end = np.mean(0.5*(mass*v**2 + kappa*(x-s)**2)[-n_rec_frac:-1])
        print traj_frac*100, "% ... ", H_S_end-H_S_beg, "kJ/mol"
    #IPython.embed()

    # Diffusion: <dx**2> = 1*D*t ; 1 for 1-Dimension
    D = np.mean((s-s[0])**2)/sim_len_time *1000.0
    print "\nDiffusion coef = ", D, "len_unit/ns \n"

