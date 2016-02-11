#!/usr/bin/env python
"""
  Extract the orientation of the cromophore CRO from XTC-traj files
  orientation: theta in degrees (pitch)
  high specificity to the problem assumed -- not easily ransferrable to other applications
"""

# coding: utf-8

# In[12]:

import mdtraj as md
import numpy as np
import math

topfile = '../conf_heated-up_1POPCmol-removed_nowater.gro'
trajname = 'traj_comp'
trajext = '.xtc'

kT = 2.5775  # kJ/mol*310K
NBINS = 90
NWINS = 20


# In[6]:

# load the first frame to set everything up
try:
    t = md.load_frame( trajname+"0"+trajext , 0, top=topfile)
except IOError:
    print "Can't read the files.."
    raise IOError
except:
    print "Something went wrong, stopping..."
    raise Exception

# define TDM within the ChROmophore
cro_inds = t.topology.select('resname CRO')
cro_inds_amino = [ i for i in cro_inds if t.top.atom(i).name in ['C1', 'C2', 'N2', 'N3', 'CA2'] ]
cro_inds_tyr   = [ i for i in cro_inds if t.top.atom(i).name in ['CG2', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'] ]



# In[14]:

def get_tdm( frame, gr1, gr2 ):
 """ returns orientation of the TDM within the ChROmophore, dirty """
# frame -- mdtraj single-framey object
# grX = indices for group-X
 coms = {}  # dictionary with COMs of the selections
 for j,sel in enumerate( (gr1,gr2) ):
    coms[j] = np.zeros(frame.xyz[0,0,:].shape)
    for i in sel:
       #print t.top.atom(i).residue, t.top.atom(i).name
       coms[j] += frame.xyz[0,i,:]
    #print coms[j]/len(sel)

 tdm = coms[0]-coms[1]
 return tdm


# In[8]:

def theta( tdm ):
   """ get the angle theta of the TDM (vector) and the z-axis ( interval [0:90] in deg) """
   norm = math.sqrt( np.square(tdm).sum() )
   t = math.acos(abs(tdm[-1])/norm) /math.pi *180.0
   return t



# In[21]:

for i in xrange(0, NWINS):
 with open("tdm_theta"+str(i)+".dat","w") as f:
  traj_filename = trajname + str(i) + trajext
  for frames in md.iterload( traj_filename , top=topfile, chunk=100):
   for frame in frames:
      thet = theta( get_tdm(frame, cro_inds_amino, cro_inds_tyr) )
      time = frame.time[0]
      line = str(time) + "  " + str(thet) + "\n"
      f.write(line)
