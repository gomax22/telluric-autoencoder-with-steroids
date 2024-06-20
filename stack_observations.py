import numpy as np
import os
from astropy.io import fits
import argparse

R = 69
ap = argparse.ArgumentParser()
ap.add_argument("-d", "--directory", required=True, help="Path to the directory containing the FITS files")
args = vars(ap.parse_args())

entries = sorted([os.path.join(args["directory"], file) \
    for file in os.listdir(args["directory"]) \
    if file.endswith('e2ds_A.fits') and not file.startswith('.')])


data = []
wave = []
airmass = []

# EXPSTART= '00:31:08' / Exposure start time (UT)
# EXPTIME =                600.0 / Effective exposure time (in secs)
# HIERARCH TNG METEO HUMIDITY = 35.0 / no comment                                 
# HIERARCH TNG METEO PRESSURE = 774.198768 / no comment                           
# HIERARCH TNG METEO WINDSPEED = 2.235 / no comment                               
# HIERARCH TNG METEO WINDDIR = 320.929269779 / no comment                         
# HIERARCH TNG METEO TEMP10M = 13.4 / no comment  
for entry in entries:
    with fits.open(entry) as hdu:
        flux = hdu[0].data
        am = hdu[0].header['AIRMASS']
                
        wavepar = np.arange(4 * R, dtype=float).reshape(R, 4) # HIERARCH TNG DRS CAL TH DEG X = 3 (N_pol + 1 = 4)
        for i in np.nditer(wavepar, op_flags=['readwrite']):
            i[...] = hdu[0].header['HIERARCH TNG DRS CAL TH COEFF LL{0}'.format(str(int(i)))]
        xx = np.arange(len(flux[0]))
        xs = [wavepar[r,0] + wavepar[r,1]*xx + wavepar[r,2]*xx**2 + wavepar[r,3]*xx**3 for r in range(R)]
        
        data.append(flux)
        airmass.append(am)
        wave.append(xs)

flux = np.array(data, dtype=np.float64)
airmass = np.array(airmass, dtype=np.float64)
wave = np.array(wave, dtype=np.float64)
fnames = np.array(entries)

out_fname = f'{args["directory"].split(os.sep)[-1]}.npz'
np.savez_compressed(os.path.join(args['directory'], out_fname), flux=flux, wave=wave, airmass=airmass, fname=fnames)