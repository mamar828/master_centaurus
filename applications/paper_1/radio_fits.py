from astropy.io import fits
from glob import glob

from src.hdu.map import Map


# Code to read radio files and convert them to FITS format
root_dir = "data/radio/greg_taylor"
files = glob("*.ICLNSH", root_dir=root_dir) + glob("*.ICLN", root_dir=root_dir) + [f"CENTAU-X12B.HIMAP"]
for file in files:
    f = fits.open(f"{root_dir}/{file}")
    data = f[0].data
    header = f[0].header
    Map(data=data, header=header).save(f"data/radio/{file.split(".")[0]}_{header["CRVAL3"]*1e-6:.0f}MHz.fits")

# Additionnal HST image
file = "U62G8401.HIMAP"
f = fits.open(f"{root_dir}/{file}")
data = f[0].data
header = f[0].header
Map(data=data, header=header).save(f"data/radio/{file.split(".")[0]}.fits")
