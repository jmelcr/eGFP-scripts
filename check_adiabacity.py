#!python 

import numpy as np
from optparse import OptionParser
#import IPython

k_b = 0.0083144621  #kJ/Mol*K

if __name__ == "__main__":
    # help message is automatically provided
    # type=string, action=store is default
    parser = OptionParser()
    parser.add_option('-f', '--file', dest='data_fname', help='data file name', default="ValMetDist_fict_bias")
    parser.add_option('-k', '--kappa', dest='kappa', help='kappa -- spring constant', default=666000, type=float)
    opts, args = parser.parse_args()

    data = np.loadtxt(opts.data_fname)
    time, x, s, b = np.hsplit(data,4)
    kappa = opts.kappa
    
    print "T_config = ", np.mean(np.square(x-s))*kappa/k_b, "K"

    W_S = np.trapz( kappa*(s-x).ravel(), s.ravel())
    print "Transferred heat = ", W_S, "kJ/mol"
    print "Heat production = ", float(W_S / (time[-1] - time[0])), "kJ/(mol*ps)"
    #IPython.embed()
