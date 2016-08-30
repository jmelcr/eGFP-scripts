#!/usr/bin/env python

# coding: utf-8

## lipidated eGFP sensors (plots)

import sys
sys.path.append('/home/joe/work/mempot/scripts')
sys.path.append('/home/pepa/work/mempot/scripts')

import cPickle
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

mymagenta2 = '#DD3377'
mymagenta  = '#991166'
mygreen2  = '#50BA33'
mygreen   = '#307720'
myblue    = '#1040AA'
myblack   = '#0A0832'
myyellow2 = '#EECD11'
myyellow  = '#EEAD22'
hlineColour = '#204020'
vlineColour = '#585100'

matplotlib.rc('font', family='TeX Gyre Adventor', size=14)

ener_unit = "kJ/mol"  #"kT" # or kJ/mol
ener_unit_factor = 2.5775  # 1.0 # or 2.5775 for kJmol


#%%

for cv in ["pitch", "tdm"]:

    # read-in the file-names of fep-profile containing files
    try:
        with open(cv+"_fep_files.list","r") as f: files_list = f.readlines()
    except:
        print "File list for CV "+cv+" probably does not exist!"
        raise RuntimeError


    # In[46]:

    # load-in fep data
    feps = {}
    tail_str = '\n'
    for fname in files_list:
        fn = fname.rstrip(tail_str)
        if "pops" in fn:
            key = "ps"
        else:
            key = "pc"
        if key in feps.keys():
            print "There already are data under the key '", key, "' in dictionary feps"
        else:
            print "Loadin data from:", fn
            feps[key] = np.loadtxt(fn)


    # Data loaded, now make abstraction aliases and the plots:

    # In[47]:

    pc = feps['pc']
    ps = feps['ps']

    crop = pc.shape[0]
    for fep in (pc,ps):
        while fep[:crop,0].max() > 90.0:
            crop -= 1


    ### pure (PO)PC .vs. (PO)PC/PS mix plots

    # In[48]:

    # plot limits
    ymax = 1.15
    ylims = (-ymax*0.5,ymax)

    #plt.xlim( (0.0, 90.0) )
    #plt.ylim( ylims )

    plt.xlabel(cv+' [deg]')
    plt.ylabel('$\Delta$ G [ '+ener_unit+' ]')

    # plotting

    plt.plot(pc[:crop,0], pc[:crop,1]*ener_unit_factor, label="PC",    lw=4, color=mygreen)
    plt.plot(ps[:crop,0], ps[:crop,1]*ener_unit_factor, label="PC/PS", lw=4, color=mymagenta)
    plt.legend(loc='upper center')

    plt.savefig(cv+"_PC-PS_"+ener_unit.replace("/","-")+".png", dpi=300)
    plt.show()

