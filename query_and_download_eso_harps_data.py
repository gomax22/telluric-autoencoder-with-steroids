import os 
import sys
from pprint import pprint
import pyvo as vo
from pyvo.dal import tap
from pathlib import Path
import argparse
import cgi
import requests
import tarfile
from astropy.io import fits
import glob
from concurrent.futures import ProcessPoolExecutor
import time

# Verify the version of pyvo 
from pkg_resources import parse_version
if parse_version(vo.__version__) < parse_version('1.4'):
    raise ImportError('pyvo version must be 1.4 or higher')
    
print('pyvo             version {version}'.format(version=vo.__version__))

# Defining the ESO tap service to query for phase 3 products:
tap = vo.dal.TAPService("https://archive.eso.org/tap_obs")


ap = argparse.ArgumentParser()
ap.add_argument("-f", "--file", required=True,
                help="path to file containing the archive ids")
ap.add_argument("-o", "--output", required=True,
                help="path to output files")
args = vars(ap.parse_args())


if not Path(args["output"]).exists():
    # create output path
    Path(args["output"]).mkdir(parents=True, exist_ok=True)

# select product_ids
with open(args['file'], 'r') as f:
    pids = f.read().splitlines()

# remove .fits
pids = [pid[:-5] for pid in pids]     

def getDispositionFilename( response ):
    """Get the filename from the Content-Disposition in the response's http header"""
    contentdisposition = response.headers.get('Content-Disposition')
    if contentdisposition == None:
        return None
    value, params = cgi.parse_header(contentdisposition)
    filename = params["filename"]
    return filename

def writeFile(response, output_dir):
    """Write on disk the retrieved file"""
    if response.status_code == 200:
        # The ESO filename can be found in the response header
        filename = getDispositionFilename( response )
        # create a folder for each archive
        # check if path exists and files already downloaded
        if not Path(output_dir).exists():
            # create output path
            Path(output_dir).mkdir(parents=True, exist_ok=True)
        

        # Let's write on disk the downloaded FITS spectrum using the ESO filename:
        with open(os.path.join(output_dir, filename), 'wb') as f:
            f.write(response.content)
        return filename 
    

def query(product_id):
    urls = {}
    # for product_id in pids:
    query = f"""
            SELECT archive_id, extension, eso_category, access_url, size_kb
            FROM phase3v2.product_files 
            where product_id='{product_id}'
        """

    res = tap.search(query)
    urls[product_id] = list(res['access_url'].data)
        
    # pprint(urls)
    return urls
    
def download(urls):

    for product_id, url in urls.items():
        for i, u in enumerate(url):
            # download file
            spec_response = requests.get(u)
            spec_fname = writeFile(spec_response, os.path.join(args['output'], product_id))
            
            # logging
            if spec_fname:
                print("Saved file: %s" % (spec_fname))
            else:
                print("Failed to download file")
    

def query_and_download_harps(pids):
    
    # get urls for the selected product_ids (spectrum + ancillary)
    urls = query(pids) 
    pprint(urls)

    download(urls)
    
def check_for_missing_files(entry):
    
    # missing files
    missing_files = []
    
    # get tar filename
    tar_file = [os.path.join(args['output'], entry, fname) for fname in os.listdir(os.path.join(args['output'], entry)) if fname.endswith('.tar')][0]

    # open tar file and extract only e2ds file
    with tarfile.open(tar_file, 'r:*') as tar:
        tar.extractall(path=os.path.join(args['output'], entry), members=[member for member in tar.getmembers() if member.name.endswith('e2ds_A.fits')])
    
    # open e2ds file
    e2ds_file = [os.path.join(args['output'], entry, fname) for fname in os.listdir(os.path.join(args['output'], entry)) if fname.endswith('e2ds_A.fits')][0]
    
    # open e2ds file and read 'HIERARCH ESO DRS BLAZE FILE' and 'HIERARCH ESO DRS CAL TH FILE' keywords
    with fits.open(e2ds_file) as hdul:
        # get blaze fname and wave fname
        blaze_file = hdul[0].header['HIERARCH ESO DRS BLAZE FILE']
        wave_file = hdul[0].header['HIERARCH ESO DRS CAL TH FILE']
    
    # add to missing files (fname.split("_")[0])
    # https://dataportal.eso.org/dataPortal/file/HARPS.2015-05-21T00:35:53.841.tar/HARPS.2015-05-21T00:35:53.841.tar
    missing_files.append(f"https://dataportal.eso.org/dataPortal/file/{blaze_file.split('_')[0]}.tar/{blaze_file.split('_')[0]}.tar")
    missing_files.append(f"https://dataportal.eso.org/dataPortal/file/{wave_file.split('_')[0]}.tar/{wave_file.split('_')[0]}.tar")
    
    # pprint({entry: missing_files})
    
    # download corresponding tar files from a different endpoint
    download({entry: missing_files})
    
    # extract only blaze file
    blaze_tar_file = os.path.join(args['output'], entry, f"{blaze_file.split('_')[0]}.tar")
    with tarfile.open(blaze_tar_file, 'r:*') as tar:
        tar.extractall(path=os.path.join(args['output'], entry), members=[member for member in tar.getmembers() if member.name.endswith('blaze_A.fits')])
    
    # extract only wave file
    wave_tar_file = os.path.join(args['output'], entry, f"{wave_file.split('_')[0]}.tar")
    with tarfile.open(wave_tar_file, 'r:*') as tar:
        tar.extractall(path=os.path.join(args['output'], entry), members=[member for member in tar.getmembers() if member.name.endswith('wave_A.fits')])
    
    # clear tar files (?)
    # for f in glob.glob(f"{os.path.join(args['output'], entry)}/*.tar"):
    #    os.remove(f)

if __name__  == "__main__":
    
    with ProcessPoolExecutor() as executor:

        print(f"started at {time.strftime('%X')}")
        
        # launch downloads over cpus
        futures = {executor.submit(query_and_download_harps, pid) for pid in pids[:10]}
    
    print(f"finished at {time.strftime('%X')}")
    print("Downloaded all files.")
    
    print("\nChecking for missing files...")
    entries = os.listdir(args['output'])
    
    with ProcessPoolExecutor() as executor:

        print(f"started at {time.strftime('%X')}")
        
        # launch downloads over cpus
        futures = {executor.submit(check_for_missing_files, entry) for entry in entries}
    
    print(f"finished at {time.strftime('%X')}")
    
    
    
    
    
    # 2. Strategy
    # open tar file and extract only e2ds file
    # open e2ds file and read 'HIERARCH ESO DRS BLAZE FILE' and 'HIERARCH ESO DRS CAL TH FILE' keywords
    # add to missing files (fname.split("_")[0])
    # repeat for all arcfiles
    # download all tar files
    # extract only wave or blaze files 
    