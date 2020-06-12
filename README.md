# pymosaics

[![PyPI version](https://badge.fury.io/py/pymosaic-fits.svg)](https://badge.fury.io/py/pymosaic-fits)
[![Build Status](https://travis-ci.org/volodymyrss/pymosaics.svg?branch=master)](https://travis-ci.org/volodymyrss/pymosaics)
[![codecov](https://codecov.io/gh/volodymyrss/pymosaics/branch/master/graph/badge.svg)](https://codecov.io/gh/volodymyrss/pymosaics)



__CAUTION__: *this is a simple functional project, although it should work as generally expected, do make sanity checks and make sure you understand what the mosaic does - as usual. In particular, the re-pixelization approach (mapping/interpolation) is suitable for maps with flux density, as long as the pixels are smaller than PSF*

Merge FITS file mosaic

* maps all images in the same pixels (taken from one of the images, or healpix)
* sums with a flexible alhorithm

* output mosaics format is compatible with input - allowing map-reduce operations

  note that pixilization used for mosaic and for the output is not neccessarily the same, i.e. for healpix it's useful to have "regular" projection output.

  clearly, heapix pixels are  useful primarily for all-sky.

## Installation

```bash
$ pip install pymosaic-fits
```

## Example:

```bash
$ mosaic tests/data/isgri_sky_ima_{1,2}.fits.gz out.fits # just a normal mosaic, pixels/output from first image
$ mosaic tests/data/isgri_sky_ima_{1,2}.fits.gz out.fits --mock # mock, to show assembly
$ mosaic tests/data/isgri_sky_ima_{1,2}.fits.gz out.fits --pixels healpix # healpix pixels, projected output
```

```python
import mosaic                                                            
mosaic.mosaic_fn_list(["tests/data/isgri_sky_ima_1.fits.gz", "tests/data/isgri_sky_ima_2.fits.gz"], "out.fits")
```
