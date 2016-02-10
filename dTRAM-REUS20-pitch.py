#!/usr/bin/env python
"""
 dTRAM analysis on the SPECIFIC REUS simulations of xxeGFPs with 20 umbrella windows at particular points in pitch-space
 this script expects a lot of problem-specific namings and conventions, definitely not flexible, nor general. 
"""


# coding: utf-8

## Pitch dTRAM/xTRAM analysis

# import everything necessary and prepare constants, variables...

# In[1]:

#get_ipython().magic(u'pylab inline')
from pytram import TRAMData, dtram, xtram # this is the dTRAM API function
import pytram
import numpy as np
import matplotlib.pyplot as plt
import cPickle
#import os
#os.chdir('/home/pepa/eGFP_sensors/wilCleGFP/wil_opls_berger/v-0.1_reus')

kT = 2.5775  # kJ/mol*310K
NBINS = 90

print """dTRAM analysing script of eGFP-like REUS simulations with 
20 umbrella windows at particular hard-coded positions.
the hard-coded tolerance for energy factors is e-6
      """

# In[2]:

print "Preparing ... "

# discretization centres assuming bin0=(0:1), bin1=(1:2) ...
centres = np.linspace( 0.5, 89.5, NBINS )

# umbrella positions
umb_pos_str = "1.0 5.0 9.5 14.0 18.5 23.0 27.5 32.0 36.5 41.0 45.5 50.0 54.5 59.0 63.5 68.0 73.0 77.5 82.0 86.5"
umb_pos = [ float(x) for x in umb_pos_str.split() ]
NWINS = len(umb_pos)

# umbrella = k/2 * dx**2   ; k=kappa as used in PLUMED
KAPPAkjm = 0.5
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

files = []
dfnm = "PITCH."
for i in range(NWINS):
   files.append(dfnm+str(i))
    


# load the PITCH.* files using pytram's reader

# In[4]:

print "Reading-in trajectories ... "

trajs = pytram.Reader(files=files)


# PITCH files don't have the right format for pytram - corrections required ...

# Now, trajs should contain list of dictionaries with Markov state 'm' and thermodyn. state 't', 
# but they contain t,m,b.
# m - time, 
# t - Markov state, 
# b - bias factor in kJ/mol

# In[5]:

# rename and reorder and fill the trajectory-dictionaries
for i,traj in enumerate(trajs.trajs):
    if traj.has_key('time'):
        print "this traj was already processed"
    else:
       traj['time'] = traj.pop('m')
       traj['m'] = traj.pop('t')
       traj['t'] = np.ones(traj['m'].shape, dtype=int) * i


# Trajectories "trajs" should now be in the right format for TRAMdata processing and then to be thrown into the DTRAM procedure: m - Markov state, t - thermodynamic state (umbrella), b - bias (always was) and there's a new key 'time'

# prepare data for dTRAM procedure:

# In[6]:

print "Preparing data for dTRAM ... "

dtramdata = TRAMData(trajs.trajs, b_K_i=b_K_i_joe)
#xtramdata = TRAMData(trajs.trajs, b_K_i=b_K_i_joe, kT_K=np.ones(NWINS)*kT, kT_target=0)


# run dTRAM/xTRAM iterative solver:

# In[43]:

# Tolerance and iterations are hardcoded

print "Running dTRAM for 1000 iters with ftol e-6" 

nruns=0

try: 
    nruns+=1
    dtram_obj_joe = dtram( dtramdata, 1, maxiter=1000, ftol=1.0E-6, verbose=False )
except: 
    print "Not converged in "+str(nruns)+"k steps. (or some other error)"
#xtram_obj_joe = xtram( xtramdata, 1, maxiter=1, ftol=1.0E-2, verbose=True )


# In[45]:

converged = False

while not converged:
   try: 
       print "Running another 1k steps ..."
       nruns+=1
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


# In[ ]:



