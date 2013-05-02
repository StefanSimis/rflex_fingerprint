# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 17:38:37 2013
@author: simis
Rflex Fingerprint algorithm (Python)

 [Please retain the following traceback notice in your code]
 Version: 20130501.1
 Adaptation: Python port, Stefan Simis

 This code is the implementation of the 'fingerprint' method to derive Rrs from hyperspectral (ir)radiance measurements:
 Simis, S.G.H. and J. Olsson. Unattended processing of shipborne hyperspectral reflectance measurements. Remote Sensing of Environment, in press. DOI: 10.1016/j.rse.2013.04.001

 <a href="http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB"><img src="http://i.creativecommons.org/l/by-sa/3.0/88x31.png"></a>
 This code is licensed under a <a href="http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB">Creative Commons Attribution-ShareAlike 3.0 Unported License</a>

Rflex hardware/software are described at http://sourceforge.net/p/rflex/wiki/Home/
This code is maintained in a git repository at http://sourceforge.net/p/rflex/fingerprint/

"""
import numpy as np
from numpy import polyfit as npf
from numpy import polyval as npv

def optimize(rho,Lt,Ls,Ed,indices,bandwidth):
    bw = bandwidth
    ind = indices
    remainder = []
    Lw = Lt - (rho*Ls)
    ratio = Lw/Ed
    grid0 = range(-bw,bw+1,1)
    grid1 = range(-bw,bw+1,1)
    grid1.remove(0)
    for j in ind:
        # fit polynome through this part of the spectrum (around the selected drop or rise)
        grid_ratio = ratio[j+np.array(grid1)]
        Y = npf(grid1,grid_ratio,2)
        Z = npv(Y,grid0)
        #define minimization parameter = distance from peak to polynomial fit
        remainder.append(abs(ratio[j]-Z[bw]))
    
#    fval_features = remainder  # fit result (not passed, scipy.optimize.fminbound expects a single float returned
    remainder = sum(remainder)  # the optimization variable
#    wlind_features = ind       # passthrough, disabled as scipy.optimize.fmindbound expects a single float
    return remainder

def getfingerprint(Lt,Ls,featureseparator,edge_width):
    # nfeatures (default 20) can be optimized if none or too many (weak features) results are generated. 
    # For sensors with UV and NIR coverage (320-950 nm) try to get 10-15 indices returned.
    nfeatures = 40
    ew = edge_width
    #This mask covers edges of the spectrum as well as any areas affected by NaNs
    #If Ed has Additional NaNs these are not considered here.
    #1 left side of mask
    edgemask = np.array([1.0]*len(Lt))
    edgemask[0:ew+np.where(np.isfinite(Ls))[0][0]]=np.nan    #1a nan edges at start of Ls
    edgemask[0:ew+np.where(np.isfinite(Lt))[0][0]]=np.nan    #1b nan edges at start of Lt
    #2 right side of mask 
    edgemask[-ew-np.where(np.isfinite(Ls[range(len(Ls)-1,-1,-1)]))[0][0]::]=np.nan #2a nan edges at end of Ls 0(2)
    edgemask[-ew-np.where(np.isfinite(Lt[range(len(Lt)-1,-1,-1)]))[0][0]::]=np.nan #2b nan edges at end of Lt 2(4)
    
    Ltdiff = abs(np.diff(Lt))
    Ltind = Ltdiff.argsort()                    #indices sorted by value of Ltdiff
    Ltind = Ltind[range(len(Ltind)-1,-1,-1)]                    #reverse order
    Ltind = Ltind[np.isfinite(edgemask[Ltind])]   #remove indices where Ltdiff is not a number

    Lsdiff = abs(np.diff(Ls)) #same for Ls
    Lsind = Lsdiff.argsort()
    Lsind = Lsind[range(len(Lsind)-1,-1,-1)]                    #reverse order
    Lsind = Lsind[np.isfinite(edgemask[Lsind])]
 
    # identify matching feature locations from both ends of the sorted first derivative 
    indices = set(Ltind[0:nfeatures]).intersection(Lsind[0:nfeatures]) #Ls and Lt signals drop at these indices 
    indices = list(indices)
    indicesLtd = np.array([Ltdiff[i] for i in indices])   # Ltdiff values associated with drops
    indicesLtd = indicesLtd/max(indicesLtd)             # normalized
    indicesLsd = np.array([Lsdiff[i] for i in indices])   # Ltdiff values associated with drops
    indicesLsd = indicesLsd/max(indicesLsd)             # normalized
    scores = list(indicesLtd+indicesLsd)                # scores for combined Ls and Lt feature strength

    nout = []
    while len(indices)>0:
        hiscore = np.where(scores==sorted(scores)[-1])[0][0]
        nout.append(indices[hiscore])
        r = range(indices[hiscore]-featureseparator,indices[hiscore]+featureseparator+1,1)
        for i in r:
            if i in indices:
                scores.pop(indices.index(i))
                indices.remove(i)

    indices = list(nout)                                      # strongest indices, separated    

    #output
    return indices
