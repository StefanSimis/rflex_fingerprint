#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calibration code: Created on Tue Apr  2 16:29:20 2013
Fingerprint code: Created on Wed May  1
@author: kathrin, stefan
"""
from sqlalchemy import *
from sqlalchemy.orm import *
from sqlalchemy.ext.declarative import declarative_base
import math
#import numpy
import scipy.optimize as spop
import rho_fingerprint as rp
import numpy as np

#set up database connection
db = create_engine('postgresql://wispweb:pw4wispweb@192.168.1.20:5432/wispweb_playground', echo=False)
metadata = MetaData(db)
Base = declarative_base(metadata=metadata)
session = sessionmaker(bind=db)()

#Reflect each database table we need to use, using metadata
class WispMeas(Base):
    __table__ = Table('wisp_measurement', metadata, autoload=True)

class Cals(Base):
    __table__ = Table('calibration', metadata, autoload=True)

def main(): 
    
    #get all calibrations
    calibrations = session.query(Cals)
    
    # parameters (calibration)
    a = np.array ([math.pi*(0.04/2)**2, math.pi*(0.04/2)**2,  math.pi*(0.39/2)**2]) #collection area for the three jazzes
    fov_gersh = 2 * math.pi * (1-math.cos(math.radians(3)/2))
    fov = np.array([fov_gersh, fov_gersh, 1])  #field of view for the three jazzes
    wl = np.arange(192,885.5,0.5)
    
    # parameters (fingerprint)
    wlNativeRes = 0.5       #spectral interval of input data (should be monotonous). Rounded to a few digits.
    fingerprint_res = 2.0   #spectral interval used to lookup fingerprint bands. For narrow band spectrometers, set to 2 nm or higher. Ignored if equal to native_res.
    wlFpRes = np.round(np.arange(wl[0],wl[-1],fingerprint_res),3)    #wavelength grid corresponding to fingerprint_res
    bandwidth_opt=int(round(7.5/wlNativeRes))                    #half-width around each future considered in optimization (in units of native wavelength grid, equivalent to 7.5 nm)
    featureseparator = int(round(10/np.mean(np.diff(wlFpRes))))       #number of channels separating subsequent features (in units of wavelength grid entered into fingerprint search, equivalent to 10 nm)
    edge_width = int(np.ceil(bandwidth_opt*wlNativeRes/fingerprint_res)) #padding around the edges of the spectrum applied in fingerprint search (in units of grid entered into search, equivalent to bandwidth_opt)
    rangeex1 = range(0,np.where(wl>=360)[0][0], 1)
    rangeex2 = range(np.where(wl>=750)[0][0], np.where(wl>=780)[0][0]+1, 1)
    range2exclude = rangeex1+rangeex2
    
    #loop through samples:
    ndata = session.query(WispMeas).count()
    ids = session.query(WispMeas.id).all()
    print '%s measurements in set' %(ndata)
    rho = [0.0]*ndata
    EmptyFlag= [0]*ndata
    HiFlag = [0]*ndata
    LoFlag = [0]*ndata
    SuspectFlag = [0]*ndata

    # loop through measurements, (1) calibrate, (2) fingerprint
    for sm in range(0,ndata,1):
        thisid = ids[sm][0]        
        wisp_meas = session.query(WispMeas).filter_by(id=thisid).first()
        print wisp_meas.id
        try:
            ID = wisp_meas.id
            # (1) Specrum calibration        
            it = np.array([wisp_meas.int_time_ld/1000000.0, wisp_meas.int_time_lu/1000000.0, wisp_meas.int_time_ed/1000000.0])
            # build count data and calibration data arrays
            # data is a 3-dimensional array holding measured spectra with 1st dim: jazes (3), 2nd dim: bands (2048), 3rd dim: parameters (0:wavelength, 1:count)
            # cals is a 2-dimensional array holding calibration data with 1st dim:  jazes (3), 2nd dim: calibration values (2048)
            cal = calibrations.filter(Cals.id == wisp_meas.calibration_ld_id)
            cal_ld = cal.one()
            cal = calibrations.filter(Cals.id == wisp_meas.calibration_lu_id)
            cal_lu = cal.one()
            cal = calibrations.filter(Cals.id == wisp_meas.calibration_ed_id)
            cal_ed = cal.one()
            
            data_ld = np.array([cal_ld.wvl, wisp_meas.counts_ld]).swapaxes(0,1)
            data_lu = np.array([cal_lu.wvl, wisp_meas.counts_lu]).swapaxes(0,1)
            data_ed = np.array([cal_ed.wvl, wisp_meas.counts_ed]).swapaxes(0,1)
            data = np.array([data_ld, data_lu, data_ed])
            cals = np.array([cal_ld.cal_values, cal_lu.cal_values, cal_ed.cal_values])
            
            radiance = np.empty((data.shape[0], data.shape[1]))
            
            # for first and last channel, wavelength spread is calculated as difference to second/next to last channel
            for i in range(data.shape[0]):
                radiance[i,0] = (0.01 * data[i,0,1] * cals[i,0]) / (it[i] * a[i] * (data[i,1,0] - data[i,0,0]) * fov[i])
                radiance[i,-1] = (0.01 * data[i,-1,1] * cals[i,-1]) / (it[i] * a[i] * (data[i,-1,0] - data[i,-2,0]) * fov[i])
                # for all other channels, spread is calculated as difference from previous to next band/2
                for j in range(1,data.shape[1]-1):
                    radiance[i,j] = (0.01 * data[i,j,1] * cals[i,j])/(it[i]*a[i]*((data[i,j+1,0] - data[i,j-1,0])/2) * fov[i])
        
            # interpolate radiance/irradiance spectra to equal intervals of 0.5nm step    
            radiance_interp05 = np.empty((data.shape[0], wl.shape[0]))
            
            for i in range(data.shape[0]):
                radiance_interp05[i,:] = np.interp(wl, data[i,:,0], radiance[i,:])
    #            wisp_meas.radiance_ld = radiance_interp05[0,:]
    #            wisp_meas.radiance_lu = radiance_interp05[1,:]
    #            wisp_meas.irradiance_ed = radiance_interp05[2,:]
            
            # (2) Fingerprint algorithm
            Ls,Lt,Ed = radiance_interp05[0,:], radiance_interp05[1,:], radiance_interp05[2,:] 
            # fingerprint search uses 2-nm resolution (use full res for optimization)        
            radiance_interp2 = np.empty((data.shape[0], wlFpRes.shape[0]))
            
            for i in range(data.shape[0]):
                radiance_interp2[i,:] = np.interp(wlFpRes, data[i,:,0], radiance[i,:],left=np.nan,right=np.nan)
            
            LsFpRes,LtFpRes = radiance_interp2[0,:], radiance_interp2[1,:]
    
            #retrieve spectral indices of the Lt&Ls fingerprint
            indices = rp.getfingerprint(LtFpRes,LsFpRes,featureseparator,edge_width)
        
            #map the positions of the fingerprint bands onto the original wavelength grid
            indicesNativeRes =[]
            for i in indices:
                match = np.where((abs(wl-wlFpRes[i]))==min(abs(wl-wlFpRes[i])))[0][0]
                indicesNativeRes.append(match)
        
            indices = indicesNativeRes    
            
            for i in set(range2exclude).intersection(indices):
                indices.remove(i)
            
            wlind = wl[indices]  #wavelength bands included in the fingerprint of this sample
            nind= len(indices)   #store number of fingerprint indices used in this sample (can use as quality filter)
        
            if nind == 0:
                EmptyFlag[sm] = 1
            else: # start minimization procedure
                #parameterize
                nonnegrange = np.intersect1d(np.where(wl>400)[0],np.where(wl<700)[0]) # spectral range where Rrs is not allowed to become negative
                rho_Hi = min(Lt[nonnegrange]/Ls[nonnegrange])
                rho_Lo = 0.024
                
                if rho_Hi<=rho_Lo:      #any solution will be negative
                    rho_Hi=rho_Lo*1.1
                    SuspectFlag[sm] = 1 
                
                #run
                thisrho = spop.fminbound(rp.optimize,x1=rho_Lo, x2=rho_Hi,\
                    args=(Lt,Ls,Ed,indices,bandwidth_opt),xtol=1e-12, maxfun=500, full_output=0, disp=0)
                
                rho[sm]=thisrho
            
                if thisrho >= rho_Hi-0.0001: # rho at rho_Hi bound, raise flag
                    HiFlag[sm] = 1
                if thisrho <= rho_Lo+0.0001:  # rho at rho_Lo bound, raise flag
                    LoFlag[sm] = 1
                
                print(sm,ID,thisrho)
        
        except:
            print 'Could not calculate new radiance for measurement %s' %(wisp_meas.id)

#    session.commit()

if __name__ == "__main__":
  main()