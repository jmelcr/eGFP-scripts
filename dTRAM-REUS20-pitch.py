#!/usr/bin/env python
"""
 dTRAM analysis on the SPECIFIC REUS simulations of xxeGFPs with 20 umbrella windows at particular points in pitch-space
 this script expects a lot of problem-specific namings and conventions, definitely not flexible, nor general.
"""


# coding: utf-8

from pytram import TRAMData, dtram # this is the dTRAM API function
import pytram
import numpy as np
import scipy
import matplotlib.pyplot as plt
import cPickle

kT = 2.5775  # kJ/mol*310K
NBINS = 90

print """dTRAM analysing script of eGFP-like REUS simulations with
20 umbrella windows at particular hard-coded positions.
the hard-coded tolerance for energy factors is e-6
      """

#%%
class SimFiles:
    """
    Class for storing meta-data about simulation files
    e.g. filename, ubrella-position etc.

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

        # list of simulation meta-data dictionaries
        self.sims = []
        for line in lines:
            split_line = line.split()
            fname  = split_line.pop(0)
            numpad = int(split_line.pop(2))
            x_0, kappa, temperature = [float(i) for i in split_line]
            tmpdict = {'fname': fname, 'numpad':numpad, 'x0': x_0, 'kappa': kappa, 'temperature': temperature}
            self.sims.append(tmpdict)

#%%

def read_trajs(dfnm, NWINS=None):
    """
    Reads in Plumed's Pitch files

    Parameters
    ----------
    dfnm : string
        default filenaming of the Pitch files
        these files typically contain several columns
        with values time, pitch, and bias
    NWINS : int
        number of simulation windows
        this number is added at the end of the filename
         If not given (is None), no number is added at the end
        and only one file is read

    Returns
    -------
    trajs : list of dictionaries containing arrays
        with the key convention of pytram
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

    files = []
    for i in range(NWINS):
        files.append(dfnm+"."+str(i))

    # load the PITCH.* files using pytram's reader
    print("Loading-in the Pitch trajectories...")
    trajs = []
    for i, filename in enumerate(files):
        tmpdict = {}
        tmpdict['time'], tmpdict['pitch'], tmpdict['b'] = np.hsplit(np.loadtxt(filename), 3)  # should contain exactly 3 columns
        tmpdict['m'] = discretize(tmpdict['pitch'])
        tmpdict['t'] = np.ones(tmpdict['m'].shape, dtype=int) * i
        trajs.append(tmpdict)
    return trajs


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
    for i, traj in enumerate(trajs.trajs):
        if 'time' in traj:
            print("this traj was already processed")
        else:
            traj['time'] = traj.pop('m')
            traj['m'] = traj.pop('t')
            traj['t'] = np.ones(traj['m'].shape, dtype=int) * i
    return trajs


# In[2]:
m = SimFiles("wham_meta.data")
m.sims

#%%

print "Preparing ... "

# discretization centres assuming bin0=(0:1), bin1=(1:2) ...
centres = np.linspace( 0.5, 89.5, NBINS )

# umbrella positions
umb_pos_str = "1.0 5.0 9.5 14.0 18.5 23.0 27.5 32.0 36.5 41.0 45.5 50.0 54.5 59.0 63.5 68.0 73.0 77.5 82.0 86.5"
umb_pos = [ float(x) for x in umb_pos_str.split() ]
NWINS = len(umb_pos)

# umbrella = k/2 * dx**2   ; k=kappa as used in PLUMED
KAPPAkjm = 0.1
KAPPAkT = KAPPAkjm/kT

# harmonic bias potential
def umb_pot( x, xk , k ):
    if None == xk:
        return 0.0
    return 0.5*k * ( x - xk )**2

# calculate the bias energies for all umbrellas and sttore them in the b_K_i_joe array
b_K_i_joe = np.zeros( shape=(NWINS, centres.shape[0]), dtype=np.float64 )
for K in xrange( NWINS ):
    b_K_i_joe[K,:] = umb_pot( centres, umb_pos[K], KAPPAkT )


# In[3]:
dfnm = "PITCH"
trajs = read_trajs(dfnm,NWINS)

print "Preparing data for dTRAM ... "

dtramdata = TRAMData(trajs, b_K_i=b_K_i_joe)



# In[43]:

# Tolerance and iterations are hardcoded

print "Running dTRAM for 1000 iters with ftol e-6"

nruns=0

try:
    #nruns+=1
    dtram_obj_joe = dtram( dtramdata, 1, maxiter=1, ftol=1.0E-6, verbose=False )
except:
    print "Not converged in "+str(nruns)+"k steps. (or some other error)"
#xtram_obj_joe = xtram( xtramdata, 1, maxiter=1, ftol=1.0E-2, verbose=True )


#%%
wham_fep_fname = "wham_fep_histo.dat"
try:
    wham_fep = np.loadtxt(wham_fep_fname)
except:
    print """File "+wham_fep_fname+" not accesible, or something else...
    using implicit zeros as the initial values"""
# set a better initial guess of f_i from WHAM:
# this is very DIRTY! f_i is meant to be private!
# FEP in kJ/mol is in the 2nd column (index 1)
dtram_obj_joe._f_i = wham_fep[:, 1]/kT
# apply the norm used in the .sc_iteration method of dtram_obj:
dtram_obj_joe._f_i += scipy.misc.logsumexp(np.append(-dtram_obj_joe._f_i, dtramdata.n_markov_states))

# run dTRAM/xTRAM iterative solver:


# In[45]:

converged = False
nruns_max = 50

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

#xtram_obj_joe.sc_iteration(maxiter=4, verbose=True)


# In[29]:

print "Saving the profile into a pickle-object and text file..."

# "aliases" for plotting and saving
#dTRAM data
cnz = (centres) #[nz]
fep_dtram = dtram_obj_joe.f_i #[nz]
fep_dtram_shift = np.copy(fep_dtram)
fep_dtram_shift -= fep_dtram_shift.min()

with open("pitch_dTRAM-FEP_windows.pickle","w") as f : cPickle.dump(dtram_obj_joe.f_K,f)

with open("pitch_dTRAM-FEP.pickle","w") as f : cPickle.dump(fep_dtram_shift,f)
#with open("pitch_dTRAM-FEP.pickle","r") as f : fep_dtram_shift = cPickle.load(f)

with open("pitch_dTRAM-FEP.dat","w") as f:
    for i,c in enumerate(cnz):
       line = str(c)+"   "+str(fep_dtram_shift[i])+"\n"
       f.write( line )



# In[46]:

print "Plotting dTRAM, and WHAM data if present ... "

# Load the result from WHAM: (generated with the same set and with the same converg. criterion 1.0E-6, 90 bins ..)
try:
    wham_data = np.loadtxt("wham_fep_histo.dat")
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
plt.xlabel( r"$x$ / deg", fontsize=12 )
plt.ylabel( r"$U(x)$ / kT", fontsize=12 )

plt.savefig("pitch_fep_dTRAM_WHAM.png", papertype="letter", dpi=300)
<<<<<<< HEAD

#%%

m = SimFiles("wham_meta.data")
=======
>>>>>>> 10f8fe10e1d1440527956a220c551413b3271657
