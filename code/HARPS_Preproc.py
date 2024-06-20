
# This notebook is for preprocessing the solar observations to make them ready for input to the autoencoder network.
# 
# The notebook performs the following preprocessing steps. 
#  1. Load Data
#  2. Connect observations with correct blaze and wave files
#  3. Blaze correct spectra
#  4. Filter out low flux and very high airmass observations 
#  5. Interpolate all observations to common wavelength grid 
#  6. Take natural logarithm and continuum normalise all observations
# 
# The notebook is set up to perform preprocessing of observations in the format of newer HARPS-N 2D solar observations. If one whises to train on solar spectra from another spectrograph this notebook will not work. You can either modify this notebook to the data format of the spectrograph or carry out your own preprocessing. The important part is to obtain blaze corrected continuum normalised log spectra interpolated to the same wavelength axis.
# ### Loading Data

import numpy as np
import pickle
from astropy.io import fits
import h5py
import time
import os
import matplotlib.pyplot as plt
from scipy import interpolate
import math
from spectrum_overload import Spectrum
from pathlib import Path

start = time.time()
print("Loading data...")

# Initializing 
data_matrix=[]
blaze_matrix=[]
wave_matrix=[]
airmass = []
berv = []
file_name = []
wave_name = []
blaze_name = []
wave_to_use = []
blaze_to_use = []

for path, subdirs, files in os.walk("/projects/data/HARPS/solar/2020"):
    for name in files:       
        # Loading Blaze files
        if name.endswith("blaze_A.fits"):
            fn = os.path.join(path, name) #filename
            blaze = fits.getdata(fn)
            blaze_matrix.append(blaze)
            blaze_name.append(fn)            

        # Loading Wave files  
        if name.endswith("wave_A.fits"):
            fn = os.path.join(path, name) #filename
            wave = fits.getdata(fn)            
            wave_matrix.append(wave)
            wave_name.append(fn)
            
        # Loading Spectrum files
        if name.endswith("e2ds_A.fits"):
            fn       = os.path.join(path, name) #filename
            tmp_data = fits.getdata(fn)       #read data entry 
            header = fits.getheader(fn)      # read header
            
            # Select only solar observations
            data_matrix.append(tmp_data)
            airmass.append(max(header['HIERARCH ESO TEL AIRM START'], header['HIERARCH ESO TEL AIRM END']))
            berv.append(header['HIERARCH ESO DRS BERV'])
            wave_to_use.append(header['HIERARCH ESO DRS CAL TH FILE'])
            blaze_to_use.append(header['HIERARCH ESO DRS BLAZE FILE'])
            file_name.append(fn)

print("Data was loaded in", time.time()-start, "seconds")
spectrum = np.asarray(data_matrix)
Airmass = np.asarray(airmass)
Berv = np.asarray(berv)
blaze = np.asarray(blaze_matrix)
wave = np.asarray(wave_matrix)
blaze_name = np.asarray(blaze_name)
wave_name = np.asarray(wave_name)
blaze_to_use = np.asarray(blaze_to_use)
wave_to_use = np.asarray(wave_to_use)

print(f"{spectrum.shape=}\n{Airmass.shape=}\n{Berv.shape=}\n{blaze.shape=}\n{wave.shape=}\n{blaze_name.shape=}\n{wave_name.shape=}\n{blaze_to_use.shape=}\n{wave_to_use.shape=}")


# ### Connecting spectra to their respective wave and blaze files

# HARPS -> 72 orders, from 378-691 nm (89-161)
# HARPS-N -> 69 orders, from 385-691 nm (92-161) (?)
# TAU is trained on HARPS-N solar observation, so we need to cut the HARPS orders to match the HARPS-N orders

spectrum, wave, blaze = spectrum[:, 3:, :], wave[:, 3:, :], blaze[:, 3:, :]

# Blaze correcting spectra
corrected_spectrum = spectrum / blaze

#Combining blaze corrected flux and wavelength into data array
Data = np.array([corrected_spectrum, wave])

# Checking if all wave and blaze files were found => no missing files
missing = []
for i in range(len(blaze)):
    if [i,0,0]==0:
        missing.append(i)
        
missing = np.array([]).astype(int)
missing = np.asarray(missing)  
if len(missing)==0:
    print('No missing blaze files')
if len(missing)>0:
    print(len(missing),'Missing blaze file')
    print(blaze_to_use[missing])
    print('For observation')
    print(missing)

missing = []
for i in range(len(wave)):
    if wave[i,0,0]==0:
        missing.append(i)
        
missing = np.array([]).astype(int)
missing = np.asarray(missing)   
if len(missing)==0:
    print('No missing wave files')
if len(missing)>0:
    print(len(missing),'Missing wave files')
    print(wave_to_use[missing])
    print('For observation')
    print(missing)

# ### Filtering Observations

# Filtering away observations with missing wave / blaze files:
if len(missing)>0:
    Data = np.delete(Data, missing ,axis=1)
    Airmass = np.delete(Airmass,missing,axis=0)
    Berv = np.delete(Berv,missing,axis=0)
    print('Number of removed spectra for missing blaze or wave files',len(missing))

# Filtering away observation with airmass > 2.0
index1 = []
for i in range(len(Airmass)):
    if  Airmass[i]  > 2: 
        index1.append(i)   

# Filtering the observations from index out 
if len(index1)>0:
    Data = np.delete(Data, index1 ,axis=1)
    Airmass = np.delete(Airmass,index1,axis=0)
    Berv = np.delete(Berv,index1,axis=0)
    file_name = np.delete(file_name,index1,axis=0)
    
    print('Number of removed spectra for high airmass:',len(index1))   

# Finding mean flux of spectra
means = np.zeros(len(Data[0]))
for i in range(len(Data[0])):
    means[i] = np.mean(Data[0,i,:])
    
#Some of the observations have very low flux values. 
# Finding index of low flux observations
index2 = []
for i in range(len(Data[0])):
    if  np.max(means) / np.mean(Data[0,i,:])  > 1.2: 
        index2.append(i)       

# Filtering the low flux observations from index out 
if len(index2)>0:
    filtered_Data = np.delete(Data, index2 ,axis=1)
    Airmass = np.delete(Airmass,index2,axis=0)
    Berv = np.delete(Berv,index2,axis=0)
    file_name = np.delete(file_name,index2,axis=0)
    
    print('Number of removed spectra for low flux:',len(index2)) 
else:
    filtered_Data = Data


# (flux/wave, observations, orders, pixels) -> (observations, orders, flux/wave, pixels)
D = filtered_Data.transpose(1,2,0,3) # Transposing data array to have observations as first index

"""
# For the specific training data we are using 
# orders 0,1,2,6 and 25 have negative flux values in their spectrum. 
# These are converted to a flux value of 1 for stability when taking the natural logarithm, 
D[:,0,0,:] = np.where(D[:,0,0,:]<1, 1, D[:,0,0,:])
D[:,1,0,:] = np.where(D[:,1,0,:]<1, 1, D[:,1,0,:])
D[:,2,0,:] = np.where(D[:,2,0,:]<1, 1, D[:,2,0,:])
D[:,6,0,:] = np.where(D[:,6,0,:]<1, 1, D[:,6,0,:])
D[:,25,0,:] = np.where(D[:,25,0,:]<1, 1, D[:,25,0,:])
"""


for i in range(len(D[0])):
    D[:, i, 0, :] = np.where(D[:, i, 0, :] < 1, 1, D[:, i, 0, :])


# Finding mean flux of spectra and saving for later scaling 
means = np.zeros(len(D))
for i in range(len(D)):
    means[i] = np.mean(D[i,:])
const = np.log(np.mean(means))

print('Data shape')
print(D.shape)

# ### Interpolating to common wavelength axis and continuum normalizing spectra

K = D.shape[1]       # number of apertures to combine (here 69)
resolution = D.shape[3] # number of pixels (here 4096)

combined_wave = np.zeros((K,resolution))
combined_flux = np.zeros((K,D.shape[0],resolution))
for k in range(K):

    # Interpolating to common restframe grid of wavepoints 
    aperture = k

    # Finding min and max for interpolation. 
    # This is needed as the observations do not cover the exact same wavelength regions
    MIN =  math.ceil(np.min(D[:,aperture,1,0]))   # Wave min
    MAX =  math.floor(np.max(D[:,aperture,1,-1])) # Wave max

    common_wave = np.linspace(MIN, MAX, num=resolution) 
    interpol_flux = []

    for i in range(len(D)):
        flux = D[i,aperture,0,:]  
        wave = D[i,aperture,1,:]
        f = interpolate.interp1d(wave, flux)

        int_flux = f(common_wave)   # use interpolation function returned by `interp1d`
        interpol_flux.append(int_flux)
    interpol_flux = np.array(interpol_flux)
    val_save = interpol_flux   
    
    # Continuum normalisation
    # Can be performed with your own choice of continuum normalization procedure
    # Here using spectrum_overload package
    interpol_flux=np.log(interpol_flux) # taking log of spectrum 
    normalized_flux = np.zeros((len(interpol_flux),len(interpol_flux[0])))
    for j in range(len(interpol_flux)):
        s = Spectrum(flux=interpol_flux[j], xaxis=common_wave)
        continuum = s.continuum(method="linear", nbins=10, ntop=5) # Optimal continuum normalisaton params can depend on spectrum size
        normalized_flux[j] = interpol_flux[j]/continuum.flux
    interpol_flux=normalized_flux
    
    # Combining the wave function and flux of each aperture
    combined_wave[k] = common_wave
    combined_flux[k] = interpol_flux
    
Preproc_wave=combined_wave
Preproc_flux=combined_flux.transpose(1,0,2) # (orders, observations, pixels) -> (observations, orders, pixels)
Preproc_airmass = Airmass
Preproc_berv = Berv

Path("/projects/data/HARPS/solar/2020/preproc/").mkdir(parents=True, exist_ok=True)

output = open('/projects/data/HARPS/solar/2020/preproc/preproc_wave.pkl', 'wb')
pickle.dump(Preproc_wave, output)
output.close()

output = open('/projects/data/HARPS/solar/2020/preproc/preproc_flux.pkl', 'wb')
pickle.dump(Preproc_flux, output)
output.close()

output = open('/projects/data/HARPS/solar/2020/preproc/preproc_airmass.pkl', 'wb')
pickle.dump(Preproc_airmass, output)
output.close()

output = open('/projects/data/HARPS/solar/2020/preproc/preproc_berv.pkl', 'wb')
pickle.dump(Preproc_berv, output)
output.close()

output = open('/projects/data/HARPS/solar/2020/preproc/preproc_const.pkl', 'wb')
pickle.dump(const, output)
output.close()

with open('/projects/data/HARPS/solar/2020/preproc/preproc_fnames.txt', 'w') as f:
    f.write("\n".join(file_name))
    