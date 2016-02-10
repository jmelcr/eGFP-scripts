#!/usr/bin/python
#  -*- mode: python; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*-
# 
#  $Id$
# 
#  Copyright (c) Erik Lindahl, David van der Spoel 2003-2007.
#  Coordinate compression (c) by Frans van Hoesel.
#  Python wrapper (c) by Roland Schulz
# 
#  IN contrast to the rest of Gromacs, XDRFILE is distributed under the
#  BSD license, so you can use it any way you wish, including closed source:
# 
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
# 
#  The above copyright notice and this permission notice shall be included in
#  all copies or substantial portions of the Software.
# 
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.
# 

#------------------------------------------------------------
# Modified by Joe,  Last edit 2014/05/16
#------------------------------------------------------------
# input: .xtc file specified as the first argument
# output: stdin dump of coordinates in a format suitable for
#        piping it into the pitch-n-roll analysis programme
#------------------------------------------------------------
# The .xtc file should contain only the atoms of interest
#  (e.g. CA carbons)
# with PBC already REMOVED!!!
#------------------------------------------------------------


from xdrfile import *
import sys

#you have to compile with --enable-shared
#and have libxdrfile.so in the LD_LIBRARY_PATH

if len(sys.argv)!=2:
  print "Missing file argument\nUsage: sample.py FILE"
  sys.exit()


x=xdrfile(sys.argv[1]) 
print x.natoms
i=0
for f in x:   #iterates frames
    #print "%8s %8s %8s %8s   Step: %8d "%("Atom","X","Y","Z",f.step) #print header
    print "%8d" % (f.time),
    for atcoor in f.x:  #iterate atoms
      print "%8.3f %8.3f %8.3f"%(atcoor[0],atcoor[1],atcoor[2]), #print coords one after another (note the comma at the end of the print command)
    print
    i+=1
    if i>=10: sys.exit()

