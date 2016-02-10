#!/usr/bin/python
#  -*- mode: python; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*-
# 

__doc__="""
This module contains classes and routines for analysis of membrane proteins (and their GROMACS simulations)
It requires PMX, a module helps with manipulation of GMX files and contained models and 
it requires xdrfile library from GROMACS site to read in XTC trajectories. 
 (must be compiled with --enable-shared and libxdrfile.so must be in the LD_LIBRARY_PATH)

------------------------------------------------------------
 Made by Joe,  Last edit 2015/03/27
------------------------------------------------------------
 input: .xtc file specified as the first argument
 output:
        
------------------------------------------------------------
 The .xtc file should contain only the atoms of interest
 with PBC already REMOVED!!!
------------------------------------------------------------
        """

import pmx
from xdrfile import *
import sys
import numpy as np

# "constants" are capitalized
NDIM = 3


class residue_with_traj():
   ''' class containes an instance of pmx.molecule.Molecule (which serves for molecules/residues)
       which is accessible here as self.mol (the parent class instance this derived class is based upon)
       with additional parameters read form trajectory (position-time, histogram, average and std) 
       used to be a subclass of pmx.molecule.Molecule, but the handy methods from there work only with Molecule, not with this class anyway.
   '''
   number = 0 #number of residues (this class)
   list = []  #list of residues (this class)

   def __init__(self, parentMol, nframes=None):
       # increment class variables
       residue_with_traj.list.append(self)
       residue_with_traj.number += 1
       # init obj vars
       if isinstance(parentMol, pmx.molecule.Molecule):
          self.mol  = parentMol
       else:
          print "Given residue/molecule is not an instance of pmx.molecule.Molecule!"
          print "The residue class no. %5d will remain empty without its parent" % (residue_with_traj.number,)
          self.mol = None
       # create numpy array if optional argument nframes provided
       if isinstance(nframes,int) and nframes>=0:
          self.traj=np.ndarray([nframes,NDIM])
       else:
          # create empty list if not
          self.traj = []

       # these obj vars exist but have no value (better than raising exceptions)
       self.histogram = None
       self.meanpos   = None
       self.meanposz  = None
       self.std       = None

   def remove(self):
       """ remove the instance from Class list of instances """
       # decrement class variables
       residue_with_traj.list.remove(self)
       residue_with_traj.number -= 1
       del self






if __name__ == '__main__':
   import IPython
   import cPickle
   # just a shorthand for a long class name
   reswtraj = residue_with_traj
   nframes = 351194
   #nframes = 400 # testing
   membrane_resnames=["POPC"]

   if len(sys.argv)!=3:
     print "Missing file argument\nUsage: sample.py CONF TRAJ"
     sys.exit()

   #open trajectory
   traj_filename = sys.argv[2]
   try:
      trajectory=xdrfile(traj_filename) 
   except:
      print "Error: Trajectory filename '%s' not right" % (traj_filename,)
      sys.exit(1)
   print "Opened trajectory with %5d atoms. " % (trajectory.natoms,)

   #make neg = Model instance based on configuration in conf_filename
   conf_filename = sys.argv[1]
   try:
      neg = pmx.model.Model(filename=conf_filename)
   except:
      print "Error: Model filename '%s' not right" % (conf_filename,)
      sys.exit(1)

   membrane = reswtraj(pmx.molecule.Molecule(),nframes=nframes)
   membrane.mol.resname="POPC"
   for res in neg.fetch_residues(membrane_resnames):
       for atom in res.atoms:
           membrane.mol.append(atom)

   # pmx dictionary that contains masses of atoms
   massDict = pmx.library._atommass
   # special atom symbols for virtual sites
   # workaround for: zero masses are automatically reset to zero - use negligible values
   massDict['M']=0.0001
   massDict['D']=0.0001

   #assign masses as they are not read-in automatically
   for atom in neg.atoms:
       atom.m = massDict[atom.symbol]

   # create list of residues with trajectory
   protResnames = list(pmx.library._protein_residues) # make a list from the tuple so that we can...
   protResnames.extend(["CRO","LAZ","LAS"]) # add the specific residues CRO and LAZ, ...
   reslist = []
   for residue in neg.residues:
       if residue.resname in protResnames:
          reslist.append(reswtraj(residue,nframes=nframes))

   #add membrane to the list
   reslist.append(membrane)

   print "Class Reswtraj now contains %5d residues." % (reswtraj.number,)

   if nframes==None:
      nframes=0
      print "Counting trajectory frames ..."
      for frame in trajectory:
          nframes+=1
      print "Trajectory contains %9d frames \n" % (i,)

   # read some frames
   i=0
   for frame in trajectory:
       # this is SLOW...
       for res in reslist:  # can also use the class reswtraj.list variable instead of appending membrane
           for atom in res.mol.atoms:
               atom.x = frame.x[atom.id-1] # there should be exact 1-1 correspondence in ordering
                                           # xtc-id starts from 0, atom.id from 1
           if isinstance(res.traj,list):
              res.traj.append(res.mol.com(vector_only=True))
           elif isinstance(res.traj, np.ndarray):
              res.traj[i] = res.mol.com(vector_only=True)

       i+=1
       if i>=nframes :
          print "Read only %6d frames, but there might be more in the trajectory." % (nframes,)
          break



   for res in reslist:
       res.meanpos = np.mean(res.traj, axis=0)
       res.meanposz = res.meanpos[2]
       
   #sort the residues by their positions
   reslist.sort(key=lambda item: item.meanposz)

   for i,res in enumerate(reslist):
       if res.mol.resname!="SOL": print i+1, res.mol.resname, res.mol.id, res.meanposz

   #IPython.embed()

   print "Pickling reslist ..."
   with open("reslistwtraj.cpickle", "w") as file: cPickle.dump(reslist,file)
   print "All done \n"
   # depickle with reslist = cPickle.load(file)

