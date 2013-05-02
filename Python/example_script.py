# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 17:38:37 2013

@author: simis
 example script for implementation of the rho_fingerprint technique
 [Please retain the following traceback notice in your final code]
 Version: 20130410.1 (git: http://git.code.sf.net/p/rflex/fingerprint)

 Adaptation: Python port by Simis

 This code is the implementation of the 'fingerprint' method to derive Rrs from hyperspectral (ir)radiance measurements:
 Simis, S.G.H. and J. Olsson. Unattended processing of shipborne hyperspectral reflectance measurements. Remote Sensing of Environment, In press (Apr 2013).

 <a href="http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB"><img src="http://i.creativecommons.org/l/by-sa/3.0/88x31.png"></a>
 This work is licensed under a <a href="http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB">Creative Commons Attribution-ShareAlike 3.0 Unported License</a>

 Rflex hardware/software are described at http://sourceforge.net/p/rflex/wiki/Home/
 
 See also rho_fingerprint_optimize rho_fingerprint_getfingerprint fminbnd

"""

#import sample data (from Matlab binary format)
import scipy.io as spio
import scipy.optimize as spop
import rho_fingerprint as rp
import numpy as np

# sample data
mat = spio.loadmat('D:/DATA/Rflex/Fingerprint algorithm on git/Python/SampleData_matlab.mat')
LsSample, LtSample, EdSample, wl = mat['LsSample'], mat['LwSample'], mat['EdSample'], mat['wl'][:,0]
ID = range(1,13,1)
wl = np.around(wl,3)

#initialize
wlNativeRes = round(np.mean(np.diff(wl)),3)                #spectral interval of input data (should be monotonous). Rounded to a few digits.
fingerprint_res = 3.3                               #spectral interval used to lookup fingerprint bands. For narrow band spectrometers, set to 2 nm or higher. Ignored if equal to native_res.
wlFpRes = np.round(np.arange(wl[0],wl[-1],fingerprint_res),3)    #wavelength grid corresponding to fingerprint_res

bandwidth_opt=int(round(7.5/wlNativeRes))                    #half-width around each future considered in optimization (in units of native wavelength grid, equivalent to 7.5 nm)
featureseparator = int(round(10/np.mean(np.diff(wlFpRes))))       #number of channels separating subsequent features (in units of wavelength grid entered into fingerprint search, equivalent to 10 nm)
edge_width = int(np.ceil(bandwidth_opt*wlNativeRes/fingerprint_res)) #padding around the edges of the spectrum applied in fingerprint search (in units of grid entered into search, equivalent to bandwidth_opt)

range2exclude = range(np.where(wl>=750)[0][0], np.where(wl>=780)[0][0]+1, 1) #note: spectrum edges will be excluded automatically in getfingerprint
nonnegrange = np.intersect1d(np.where(wl>400)[0],np.where(wl<700)[0]) # spectral range where Rrs is not allowed to become negative
rho_Lo = 0.024

#loop through samples:
rho = [0.0]*len(ID)
EmptyFlag= [0]*len(ID)
HiFlag = [0]*len(ID)
LoFlag = [0]*len(ID)
SuspectFlag = [0]*len(ID)
wlind = [0]*len(ID)
nind = [0]*len(ID)

for j in range(0,len(ID),1):
    Ls,Lt,Ed = LsSample[:,j], LtSample[:,j], EdSample[:,j]
    
    #optional interpolation (only for feature identification)
    if wlNativeRes != fingerprint_res:
        LtFpRes = np.interp(wlFpRes,wl,Lt,left=np.nan,right=np.nan)
        LsFpRes = np.interp(wlFpRes,wl,Ls,left=np.nan,right=np.nan)
        EdFpRes = np.interp(wlFpRes,wl,Ed,left=np.nan,right=np.nan)
    else:
        LtFpRes = Lt
        LsFpRes = Ls
        EdFpRes = Ed
    
    #retrieve spectral indices of the Lt&Ls fingerprint
    indices = rp.getfingerprint(Lt,Ls,featureseparator,edge_width)

    #map the positions of the fingerprint bands onto the original wavelength grid
    indicesNativeRes =[]
    for i in indices:
        match = np.where((abs(wl-wlFpRes[i]))==min(abs(wl-wlFpRes[i])))[0][0]
        indicesNativeRes.append(match)

    indices = indicesNativeRes    
    
    for i in set(range2exclude).intersection(indices):
        indices.remove(i)
        
    wlind[j] = wl[indices]  #wavelength bands included in the fingerprint of this sample
    nind[j]= len(indices) #store number of fingerprint indices used in this sample (can use as quality filter)

    if nind[j] == 0:
        EmptyFlag[j] = 1
    else: # start minimization procedure
        #parameterize
        rho_Hi = min(Lt[nonnegrange]/Ls[nonnegrange])
        
        if rho_Hi<=rho_Lo:      #any solution will be negative
            rho_Hi=rho_Lo*1.1
            SuspectFlag[j] = 1 
        
        #run
        thisrho = spop.fminbound(rp.optimize,x1=rho_Lo, x2=rho_Hi,\
            args=(Lt,Ls,Ed,indices,bandwidth_opt),xtol=1e-12, maxfun=500, full_output=0, disp=0)
        rho[j]=thisrho
    
        if thisrho >= rho_Hi-0.0001:  # rho at rho_Hi bound, raise flag
            HiFlag[j] = 1
        if thisrho <= rho_Lo+0.0001:  # rho at rho_Lo bound, raise flag
            LoFlag[j] = 1
        
#        print wlind[j]
        
#print(ID,rho)
#print(LoFlag)
#print(HiFlag)
#print(SuspectFlag)
#print(EmptyFlag)
