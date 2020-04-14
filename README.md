# pymosaics

Merge FITS file mosaic

* maps all images in the same pixels (taken from one of the images, or healpix)
* sums with a flexible alhorithm

* output mosaics format is compatible with input - allowing map-reduce operations

Example:

```bash
$ mosaic tests/data/isgri_sky_ima_{1,2}.fits.gz  t.fits
```

```python
import mosaic$                                                               
mosaic.mosaic_fn_list(["tests/data/isgri_sky_ima_1.fits.gz", "tests/data/isgri_sky_ima_2.fits.gz"], "out.fits")
```

