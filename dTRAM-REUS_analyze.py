#!/usr/bin/env python
"""
 dTRAM analysis on the SPECIFIC REUS simulations of xxeGFPs
 Beware: this script might still be a bit too specific to this task.

 ------------------------------------------------------------
 Made by Joe,  Last edit 2016/03/15
------------------------------------------------------------
 input: Plumed's pitch files,
        MetaData in Grossfield's-WHAM format
 output: FEP (pickled, textfile and plot)

--------------------------------------------------------
"""


# coding: utf-8

from pytram import TRAMData, dtram # this is the dTRAM API function
import pytram
import numpy as np
import scipy
import matplotlib.pyplot as plt
import cPickle
from optparse import OptionParser

k_b = 0.0083144621  #kJ/Mol*K


#%%
class SimFiles:
    """
    Class for storing meta-data about simulation files
    e.g. filename, ubrella-position etc.
    and reading-in and manipulating simulation data

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
                tmpdict = {'fname': fname, 'numpad': numpad,
                           'x0': x_0, 'kappa': kappa,
                           'temperature': temperature}
                self.sims.append(tmpdict)
            else:
                self.nsims -= 1


    def read_trajs(self, binwidth):
        """
        Reads in Plumed's Pitch files

        File format:
        ------------
        3 columns : floats
            Time -- PITCH -- bias
            crashes or unexpected behaviour if different format (e.g. more cols)

        Returns
        -------
        trajs : list of dictionaries containing arrays
            with the key convention of pytram
            dict-key 'm' -- Markov state -- 0-based indexing
        """
        def discretize(in_data, offset=0.0):
            """
            turns data into integers, i.e.
            discretizes data based on the offset and bin width
            Discretized data are 0-based indices
            (shifted with local var smallest)

            Parameters
            ----------
            data : np.array
                data to be discretized
            offset : float
                starting value for the discretization scheme
            """
            data = np.copy(in_data)
            if offset != 0.0:
                data -= offset
            # binwidth variable taken from the
            # enclosing scope (at definition time) -- closure
            if binwidth != 1.0:
                data /= binwidth
            return data.astype(int)

        smallest = None
        trajs = []
        for i, sim in enumerate(self.sims):
            filename = sim['fname']
            tmpdict = {}
            tmpdict['time'], tmpdict['x'], tmpdict['b'] = np.hsplit(np.loadtxt(filename), 3)  # should contain exactly 3 columns
            tmpdict['m'] = discretize(tmpdict['x'])
            # find the smallest state-no
            if  smallest == None or smallest > tmpdict['m'].min():
                smallest = tmpdict['m'].min()
            tmpdict['t'] = np.ones(tmpdict['m'].shape, dtype=int) * i
            trajs.append(tmpdict)

        # shift the trajs by the smallest state-no (smallest var)
        for traj in trajs:
            traj['m'] -= smallest

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


class SimData(SimFiles):
    """
    SubClass of SimFiles for storing trajectories-data and
    meta-data about simulation files
    e.g. filename, ubrella-position etc.
    and readsing-in and manipulating simulation data

    the expected file-format is the same as for Grossfields WHAM code

    See the code for further details
    """
    def __init__(self, fname, bin_width, verbose=False):
        SimFiles.__init__(self, fname)
        if isinstance(bin_width, float):
            self.verbose = verbose
            self.bin_width = bin_width
            self.trajs = self.read_trajs(binwidth=self.bin_width)
            self.remove_unvisited()
            self.gridpoints = self.get_gridpoints()
            self._n_therm_states = None
            self._n_markov_states = None
        else:
            print "Incorrect bin_width provided: ", bin_width
            raise RuntimeError, "Terminating with exception ..."

    def remove_unvisited(self):
        """
        Returns the same list of trajectories with
        reindexed Markov states,
        so that there are no unvisited states

        Parameters:
        -----------
        trajs : list
            list of trajectories (dictionaries in pytram-format)
        """
        # store the indices in an attribute
        self.unvisited_states = self.get_unvisited_states()
        # shift the Markov-states indices
        self.shift_m_indices()

    def shift_m_indices(self):
        """
        Shifts Markov-state indices IN-PLACE so that there are no
        unvisited states (stored in attribute for possible reversal)
        and the state-index series is continuous (0,1,2,3,4...)
        """
        if not isinstance(self.unvisited_states, list):
            print "Did you already look for unvisited states?"
        elif len(self.unvisited_states) > 0:
            # sort unvisited states descending -> req'd for subseq. shifting
            self.unvisited_states.sort(reverse=True)
            for i in self.unvisited_states:
                for traj in self.trajs:
                    # select all states with a higher index than i
                    sel = (traj['m'] > i)
                    # reduce their index-number by 1
                    traj['m'][sel] -= 1
        else:
            print "shift_m_indices: Nothing to shift/remove."

    def get_trajs_minmax(self, key='x'):
        """
        Returns the min/max values of key='x' (some coordinate, e.g. pitch)
        from list of trajs (each list-item is a dictionary -- pyTram format)
        """
        # initialize with the first trajectory
        traj = self.trajs[0]
        minim_all, maxim_all = get_minmax(traj[key])
        # go through the rest trajs
        for traj in self.trajs[1:]:
            minim, maxim = get_minmax(traj[key])
            if minim < minim_all:
                minim_all = minim
            if maxim > maxim_all:
                maxim_all = maxim

        return (minim_all, maxim_all)

    def get_unvisited_states(self):
        """
        Returns the Markov-indices of unvisited states (sorted) in a list
        Assuming 0-alignment of bins
        e.g. bin_width=1.0 gives bins (0.0:1.0), (1.0,2.0), ...
        """
        # find min/max of Markov state indices
        minim, maxim = self.get_trajs_minmax(key='m')

        # find the unvisited states and put them into list
        if minim > 0:
            unvisited_states = range(minim)
        else:
            unvisited_states = []

        for j in range(minim, maxim):
            visited = False
            for traj in self.trajs:
                if (traj['m'] == j ).any():
                    visited = True
                    break
            if not visited:
                unvisited_states.append(j)
        if len(unvisited_states) > 0:
            unvisited_states.sort()
        return unvisited_states

    def get_avg_bias_mtx(self):
        """
        Return the b_K_i matrix (required for pytram)
        obtained by averaging the biasfactors from simulation data
        (does not require any assumption on the potential shape and its params)
        """
        print "TODO!!!"

    def get_gridpoints(self):
        """
        Returns 0-aligned gridpoints according to the bin_width and
        min/max of all trajectories
        """
        # get lower/upper bounds of the trajectories (all data)
        bounds_minmax = self.get_trajs_minmax()
        # bounds will contain the min/max gridpoints/bin-centres
        # instead of absolute real min/max
        bounds = list(bounds_minmax)
        for i, bound in enumerate(bounds):
            bounds[i] = get_bin_centre(bound, self.bin_width)
        # how many bins do we need then?
        self.nbins = int(abs(bounds[1]-bounds[0]) / self.bin_width) + 1  # +1 for the right edge
        # generate gripoints aligned with "absolute 0"
        gridpoints = np.linspace( bounds[0], bounds[1], self.nbins )
        if self.unvisited_states > 0:
            sel = np.array(gridpoints, dtype=bool, copy=False)
            sel[:] = True
            sel[self.unvisited_states] = False
            gridpoints = gridpoints[sel]

        return gridpoints


##########################
#  Copied from dTRAM code
##########################

    def get_C_K_ij( self, lag=10, sliding_window=True ):
        r"""
        Parameters
        ----------
        lag : int
            lagtime tau, at which the countmatrix should be evaluated
        sliding_window : boolean (default=True)
            lag is applied by mean of a sliding window or skipping data entries.

        Returns
        -------
        C_K_ij : numpy.ndarray(shape=(T,M,M))
            count matrices C_ij at each termodynamic state K
        """
        C_K_ij = np.zeros(
            shape=(self.n_therm_states, self.n_markov_states, self.n_markov_states),
            dtype=np.intc)
        for traj in self.trajs:
            t = 0
            while t < traj['m'].shape[0]-lag:
                K = traj['t'][t]
                if np.all(traj['t'][t:t+lag+1] == K):
                    C_K_ij[K, traj['m'][t], traj['m'][t+lag]] += 1
                if sliding_window:
                    t += 1
                else:
                    t += lag
        return C_K_ij

    ############################################################################
    #
    #   n_markov_states / n_therm_states getters
    #
    ############################################################################

    @property
    def n_markov_states(self):
        if self._n_markov_states is None:
            if self.verbose:
                print "# Counting Markov states"
            self._n_markov_states = 0
            for traj in self.trajs:
                max_state = np.max(traj['m'])
                if max_state > self._n_markov_states:
                    self._n_markov_states = max_state
            self._n_markov_states += 1
            if self.verbose:
                print "# ... found %d Markov states" % self._n_markov_states
        return self._n_markov_states

    @property
    def n_therm_states(self):
        if self._n_therm_states is None:
            if self.verbose:
                print "# Counting thermodynamic states"
            self._n_therm_states = 0
            for traj in self.trajs:
                max_state = np.max(traj['t'])
                if max_state > self._n_therm_states:
                    self._n_therm_states = max_state
            self._n_therm_states += 1
            if self.verbose:
                print "# ... found %d thermodynamic states" % self._n_therm_states
        return self._n_therm_states



#%%

def get_minmax(data):
    """ Returns (minim, maxim) tuple of data (np.array) """
    minim = data.min()
    maxim = data.max()
    return (minim, maxim)

def get_bin_centre(value, bin_width):
    """
    Returns centre (i.e. gridpoint) of
    the bin the value belongs to
    Gridpoints/bins assumed 0-aligned
    """
    return (int(value/bin_width)+0.5)*bin_width


def get_init_f_i_from_wham(wham_fep_fname, sim_data):
    """
    Reads in the WHAM-guess of f_i factors
    and returns them as a new guess for assignment to dTRAM obj attribute ._f_i
    -- sets a better initial guess of f_i from WHAM,
    BUT: the gridpoints don't have to necessarily agree!
    """
    try:
        wham_fep = np.loadtxt(wham_fep_fname)
        # this is very DIRTY! f_i is meant to be private!
        # FEP in kJ/mol is in the 2nd column (index 1)
        new_f_i = wham_fep[:, 1] / (k_b*sim_data.sims[0]['temperature'])
        # apply the norm used in the .sc_iteration method of dtram_obj:
        new_f_i += scipy.misc.logsumexp(np.append(-new_f_i, sim_data. len(new_f_i)))
    except:
        print """File "+wham_fep_fname+" not accesible, or something else...
        using implicit zeros as the initial values"""

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

def weight_switch(vals, thres):
    """
    applies a switching function on an array element-wise for
    elements that are above given thershold

    For values that are delta above threshold,
    the switching function is:
      thres + delta*math.exp(-delta/thres)

    Parameters
    ----------
    vals : array
        input values for the switching function
    thres : float
        threshold value for the switching function

    Returns
    -------
    weight_switch : np.array
        of switched/corrected values (floats)
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

if __name__ == "__main__":
    # help message is automatically provided
    # type=string, action=store is default
    parser = OptionParser()
    parser.add_option('-m', '--metafile', dest='wham_metadata_fname', help='wham metadata file name', default="wham_meta.data")
    parser.add_option('-w', '--whamfep', dest='wham_fep_fname', help='wham FEP file name', default="wham_fep_histo.dat")
    parser.add_option('-b', '--binw', dest='bin_width', help='bin width', default=1.0, type=float)
    parser.add_option('-i', '--niter', dest='nruns_max', help='no. iter*1k (=1run)', default=100, type=int)
    parser.add_option('-l', '--lag', dest='lag_time', help='lag time (in units of frame-rec-rate) for dTRAM c_k_i_j kinetic mtx', default=10, type=int)
    parser.add_option('-t', '--ftol', dest='f_toler_dtram', help='convergence tolerance for dTRAM iterative sc-procedure', default=1.0E-8, type=float)
    opts, args = parser.parse_args()

    #read-in the simulations' metadata in the Grossfield's WHAM format
    # and the trajectories in the PLUMED's format
    sim_data = SimData(opts.wham_metadata_fname, opts.bin_width)

    # just an alias, gridpoints are saved as an attribute in the SimData obj
    gridpoints = sim_data.gridpoints

#%%

######################################
#  C_K_ij trans-mtx Reweighting part
######################################


    with open("pitch_dTRAM-FEP_windows.pickle","r") as f: winbiases = cPickle.load(f)

#%%

    ckij = sim_data.get_C_K_ij()

#%%
    bki = sim_data.gen_harmonic_bias_mtx(gridpoints)
    # inserting the windows-biases into the b_K_i matrix
#%%
    # the result of dTRAM shoud not (and does not) depend on this
    for k in range(ckij.shape[0]):
        bki[k,:] -= winbiases[k]

#%%

    dtrammtx = pytram.dtram_from_matrix(ckij,bki, maxiter=5000)
#%%
    dtrammtx.sc_iteration(maxiter=5000, ftol=1.0E-5)
    plt.plot(gridpoints, dtrammtx.f_i) #-fep_dtram)
    plt.savefig("figtest1.png",dpi=200)
    plt.show()

#%%



#%%

################################
#  OLD code continues...
################################

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



    # In[29]:

    print "Saving the profile into a pickle-object and text file..."

    # "aliases" for plotting and saving
    #dTRAM data
    fep_dtram = dtram_obj.f_i
    fep_dtram_shift = np.copy(fep_dtram)
    fep_dtram_shift -= fep_dtram_shift.min()

    with open("pitch_dTRAM-FEP_windows.pickle","w") as f : cPickle.dump(dtram_obj.f_K,f)
    with open("pitch_dTRAM-FEP.pickle","w") as f : cPickle.dump(fep_dtram_shift,f)
    with open("pitch_dTRAM-FEP.dat","w") as f:
        for i,c in enumerate(gridpoints):
           line = str(c)+"   "+str(fep_dtram_shift[i])+"\n"
           f.write( line )



    # In[46]:

    print "Plotting dTRAM, and WHAM data if present ... "

    # Load the result from WHAM: (generated with the same set and with the same converg. criterion 1.0E-6, 90 bins ..)
    try:
        wham_data = np.loadtxt(opts.wham_fep_fname)
        cnz_wham = wham_data[:,0]
        kT = k_b*sim_data.sims[0]['temperature']
        fep_wham = wham_data[:,1]/kT # scale by kT
        fep_wham_err1 =  wham_data[:,2]/kT + fep_wham
        fep_wham_err2 = -wham_data[:,2]/kT + fep_wham
        plt.plot( cnz_wham, fep_wham, '-', color='green', lw=2.0, label="WHAM" )
        plt.fill_between( cnz_wham, fep_wham_err1, fep_wham_err2, color='green', alpha=0.2)
    except:
        print "Wham data missing or just something wrong happened in the process. Look into the code, lad. "

    plt.plot( gridpoints, fep_dtram_shift, '-', color='black', lw=2.0, label="dTRAM" )

    plt.legend( loc='upper left', fontsize=10 )
    plt.xlabel( r"$pitch$ [deg]", fontsize=12 )
    plt.ylabel( r"$G$ [kT]", fontsize=12 )

    plt.savefig("pitch_fep_dTRAM_WHAM.png", papertype="letter", dpi=300)
