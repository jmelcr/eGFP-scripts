#!/usr/bin/env python
"""
 dTRAM analysis on the SPECIFIC REUS simulations of xxeGFPs
 this script expects a lot of problem-specific namings and
 conventions, definitely not flexible, nor general.

 TODO: can guess gridpoints from the visited states
       should be non-visited state proof (program crashes if there's one)
"""


# coding: utf-8

from pytram import TRAMData, dtram # this is the dTRAM API function
import pytram
import numpy as np
import scipy
import matplotlib.pyplot as plt
import cPickle

k_b = 0.0083144621  #kJ/Mol*K
kT = 2.5775  # kJ/mol*310K
NBINS = 90
lbound = 0.5
ubound = 89.5
nruns_max = 50
lag_time = 100
f_toler_dtram = 1.0E-6
wham_fep_fname = "wham_fep_histo.dat"
wham_metadata_fname = "wham_meta.data"


#%%
class SimFiles:
    """
    Class for storing meta-data about simulation files
    e.g. filename, ubrella-position etc.
    and readin-in and manipulation simulation data

    the expected file-format is the same as for Grossfields WHAM code

    See the code for further details
    """
    def __init__(self, fname):
        if not isinstance(fname, str):
            raise RuntimeError, "provided filename is not a string!"
        self.fname = fname
        try:
            metadatafile = open(self.fname, "r")
            lines = metadatafile.readlines()
            metadatafile.close()
        except:
            raise IOError, "Can't open/close/read the provided file!"

        self.nsims = len(lines)

        # list of simulation meta-data dictionaries
        self.sims = []
        for line in lines:
            if not ( line.startswith("#") and line.strip() ) :
                split_line = line.split()
                fname  = split_line.pop(0)
                numpad = int(split_line.pop(2))
                x_0, kappa, temperature = [float(i) for i in split_line]
                # kappa is twice that large for Grossfields WHAM code than
                # in the definition here + convert it to kT units from kJ/mol
                kappa /= 2.0*k_b*temperature
                tmpdict = {'fname': fname, 'numpad':numpad, 'x0': x_0, 'kappa': kappa, 'temperature': temperature}
                self.sims.append(tmpdict)
            else:
                self.nsims -= 1


    def read_trajs(self):
        """
        Reads in Plumed's Pitch files

        Returns
        -------
        trajs : list of dictionaries containing arrays
            with the key convention of pytram

        TODO:
        -----
        dTRAM fails if the trajectories contain unvisited states
          -- not checked here
        """
        def discretize(data, binwidth=1.0, offset=0.0):
            """
            turns data into integers, i.e.
            discretizes data based on the offset and bin width

            Parameters
            ----------
            data : np.array
                data to be discretized
            binwidth : float
                width of one bin in the discretization scheme
            offset : float
                starting value for the discretization scheme
            """
            if offset != 0.0:
                data -= offset
            if binwidth != 1.0:
                data /= binwidth
            return data.astype(int)

        trajs = []
        for i, sim in enumerate(self.sims):
            filename = sim['fname']
            tmpdict = {}
            tmpdict['time'], tmpdict['pitch'], tmpdict['b'] = np.hsplit(np.loadtxt(filename), 3)  # should contain exactly 3 columns
            tmpdict['m'] = discretize(tmpdict['pitch'])
            tmpdict['t'] = np.ones(tmpdict['m'].shape, dtype=int) * i
            trajs.append(tmpdict)
        return trajs


    def gen_harmonic_bias_mtx(self, gridpoints):
        """ calculate the bias energies for all umbrellas and sttore them in the b_K_i array """
        def umb_pot( x, xk , k ):
            """ harmonic bias potential """
            if None == xk:
                return 0.0
            return 0.5*k * ( x - xk )**2

        b_K_i = np.zeros( shape=(self.nsims, gridpoints.shape[0]), dtype=np.float64 )
        for K in xrange( self.nsims ):
            b_K_i[K,:] = umb_pot( gridpoints, self.sims[K]['x0'], self.sims[K]['kappa'] )
        return b_K_i



#%%

if __name__ == "__main__":
    #read-in the simulations' metadata in the Grossfield's WHAM format
    metadata = SimFiles(wham_metadata_fname)

    # this is BAD! -- generalize!
    # discretization centres assuming bin0=(0:1), bin1=(1:2) ...
    #lbound = 80
    #ubound = 180
    #NBINS = ubound - lbound
    gridpoints = np.linspace( lbound, ubound, NBINS )

    trajs = metadata.read_trajs()

    dtramdata = TRAMData(trajs, b_K_i=metadata.gen_harmonic_bias_mtx(gridpoints=gridpoints))

    try:
        dtram_obj_joe = dtram( dtramdata, lag=lag_time, maxiter=1, ftol=f_toler_dtram, verbose=False )
    except pytram.ExpressionError:
        raise pytram.ExpressionError, "Input was faulty, ending..."
    except pytram.NotConvergedWarning:
        pass


    #%%

    # set a better initial guess of f_i from WHAM:
    try:
        wham_fep = np.loadtxt(wham_fep_fname)
    except:
        print """File "+wham_fep_fname+" not accesible, or something else...
        using implicit zeros as the initial values"""
    # this is very DIRTY! f_i is meant to be private!
    # FEP in kJ/mol is in the 2nd column (index 1)
    dtram_obj_joe._f_i = wham_fep[:, 1] / (k_b*metadata.sims[0]['temperature'])
    # apply the norm used in the .sc_iteration method of dtram_obj:
    dtram_obj_joe._f_i += scipy.misc.logsumexp(np.append(-dtram_obj_joe._f_i, dtramdata.n_markov_states))


    # In[45]:

    converged = False
    nruns=0

    while not converged:
       try:
           print "Running another 1k steps ..."
           nruns+=1
           if nruns > nruns_max:
               print """BEWARE: Did not converge within 50k iterations, this is weird, ending. \nBEWARE: generated results are based on unconverged data! """
               break
           dtram_obj_joe.sc_iteration(maxiter=1000)
       except:
           print "Not converged in "+str(nruns)+"k steps. ",
       else:
           print "Converged within "+str(nruns)+"k steps."
           converged = True



    # In[29]:

    print "Saving the profile into a pickle-object and text file..."

    # "aliases" for plotting and saving
    #dTRAM data
    cnz = (gridpoints) #[nz]
    fep_dtram = dtram_obj_joe.f_i #[nz]
    fep_dtram_shift = np.copy(fep_dtram)
    fep_dtram_shift -= fep_dtram_shift.min()

    with open("pitch_dTRAM-FEP_windows.pickle","w") as f : cPickle.dump(dtram_obj_joe.f_K,f)
    with open("pitch_dTRAM-FEP.pickle","w") as f : cPickle.dump(fep_dtram_shift,f)
    with open("pitch_dTRAM-FEP.dat","w") as f:
        for i,c in enumerate(cnz):
           line = str(c)+"   "+str(fep_dtram_shift[i])+"\n"
           f.write( line )



    # In[46]:

    print "Plotting dTRAM, and WHAM data if present ... "

    # Load the result from WHAM: (generated with the same set and with the same converg. criterion 1.0E-6, 90 bins ..)
    try:
        wham_data = np.loadtxt(wham_fep_fname)
        cnz_wham = wham_data[:,0]
        fep_wham = wham_data[:,1]/kT # scale by kT
        fep_wham_err1 =  wham_data[:,2]/kT + fep_wham
        fep_wham_err2 = -wham_data[:,2]/kT + fep_wham
        plt.plot( cnz_wham, fep_wham, '-', color='green', lw=2.0, label="WHAM" )
        plt.fill_between( cnz_wham, fep_wham_err1, fep_wham_err2, color='green', alpha=0.2)
    except:
        print "Wham data missing or just something wrong happened in the process. Look into the code, lad. "

    plt.plot( cnz, fep_dtram_shift, '-', color='black', lw=2.0, label="dTRAM" )

    plt.legend( loc='upper left', fontsize=10 )
    plt.xlabel( r"$pitch$ [deg]", fontsize=12 )
    plt.ylabel( r"$G$ [kT]", fontsize=12 )

    plt.savefig("pitch_fep_dTRAM_WHAM.png", papertype="letter", dpi=300)
