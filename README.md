# pymosaics

Merge FITS file mosaic

* maps all images in the same pixels (taken from one of the images, or healpix)
* sums with a flexible alhorithm

* output mosaics format is compatible with input - allowing map-reduce operations

  note that pixilization used for mosaic and for the output is not neccessarily the same, i.e. for healpix it's useful to have "regular" projection output.

  clearly, heapix pixels are  useful primarily for all-sky.

## Installation

```bash
$ pip install pymosaics
```

## Example:

```bash
$ mosaic tests/data/isgri_sky_ima_{1,2}.fits.gz out.fits # just a normal mosaic, pixels/output from first image
$ mosaic tests/data/isgri_sky_ima_{1,2}.fits.gz out.fits --mock # mock, to show assembly
$ mosaic tests/data/isgri_sky_ima_{1,2}.fits.gz out.fits --pixels healpix # healpix pixels, projected output
```

```python
import mosaic$                                                               
mosaic.mosaic_fn_list(["tests/data/isgri_sky_ima_1.fits.gz", "tests/data/isgri_sky_ima_2.fits.gz"], "out.fits")
```
