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
import argparse
import gc
import shutil

def check_ancillary(input_dir: str):
    fnames = os.listdir(input_dir)
    e2ds = wave = blaze = False
    for fname in fnames:
        if fname.endswith("e2ds_A.fits"):
            e2ds = True
        if fname.endswith("wave_A.fits"):
            wave = True
        if fname.endswith("blaze_A.fits"):
            blaze = True
    return e2ds and wave and blaze

ap = argparse.ArgumentParser()
ap.add_argument("-p", "--path", required=True, help="Path to the HARPS solar data")
ap.add_argument("-r", "--remove", required=False, default=False, help="Remove missing entries", type=bool)
args = vars(ap.parse_args())

start = time.time()
print("Loading data...")

# initialization
flux    = []
wave    = []
blaze   = []
fnames  = []
airmass = []
berv    = []
missing = []

# check all dirs of the path
entries = [os.path.join(args['path'], entry) for entry in os.listdir(args['path']) if os.path.isdir(os.path.join(args['path'], entry))]
for idx, entry in enumerate(entries):
    if idx % 100 == 0: print(f"Loading entry n. {idx} of {len(entries)}")
    
    # check ancillary files
    flag = check_ancillary(entry)
    if not flag:
        print(f"Skipping {entry} due to missing files")
        missing.append(entry.split('/')[-1])
        continue
    """
    # track spectrum, wave, blaze, airmass and berv
    for path, subdirs, files in os.walk(entry):
        for name in files:
            fn = os.path.join(path, name)
            fnames.append(fn)
            
            # loading Spectrum files
            if name.endswith("e2ds_A.fits"):
                with fits.open(fn) as hdul:
                    flux.append(hdul[0].data)
                    header = hdul[0].header
                    airmass.append(max(header['HIERARCH ESO TEL AIRM START'], header['HIERARCH ESO TEL AIRM END']))
                    berv.append(header['HIERARCH ESO DRS BERV'])
            
            # loading Blaze files
            if name.endswith("blaze_A.fits"):
                with fits.open(fn) as hdul:
                    blaze.append(hdul[0].data)
            
            # loading Wave files
            if name.endswith("wave_A.fits"):
                with fits.open(fn) as hdul:
                    wave.append(hdul[0].data)
                    
flux = np.stack(flux)
blaze = np.stack(blaze)
wave = np.stack(wave)
airmass = np.stack(airmass)
berv = np.stack(berv)
"""
print("Data was loaded in", time.time()-start, "seconds")
print(f"Missing {len(missing)} entries")

print(f"\nMissing entries: {missing}")
missing = [os.path.join(args['path'], entry) for entry in missing]

if args['remove']:
    for entry in missing:
        shutil.rmtree(entry, ignore_errors=True)
    print("Missing entries were removed successfully.")