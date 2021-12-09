import os
from astropy.io import fits
    
def test_healpix_single_mosaima():
    import mosaic
    out_fn = "out_p1.fits"
    mosaic.mosaic_list(["tests/data/isgri_mosa_ima_1.fits.gz"], out_fn, pixels="healpix", mock=False)


def test_healpix_mosaima():
    import mosaic
    out_fn = "out_p1.fits"
    mosaic.mosaic_list(["tests/data/isgri_mosa_ima_1.fits.gz", "tests/data/isgri_mosa_ima_2.fits.gz"], out_fn, pixels="healpix", mock=False)

def test_healpix_mosaima_mock():
    import mosaic
    out_fn = "out_p1.fits"
    mosaic.mosaic_list(["tests/data/isgri_mosa_ima_1.fits.gz", "tests/data/isgri_mosa_ima_2.fits.gz"], out_fn, pixels="healpix", mock=True)

def test_healpy():
    pass

def test_first_mosaima_mock():
    import mosaic
    out_fn = "out.fits"
    mosaic.mosaic_list(["tests/data/isgri_mosa_ima_1.fits.gz", "tests/data/isgri_mosa_ima_2.fits.gz"], out_fn, pixels="first", mock=True)

    import os
    assert os.path.exists(out_fn)

def test_first_skyima_mock():
    import mosaic
    out_fn = "out.fits"
    mosaic.mosaic_list(["tests/data/isgri_sky_ima_1.fits.gz", "tests/data/isgri_sky_ima_2.fits.gz"], out_fn, pixels="first", mock=True)

    import os
    assert os.path.exists(out_fn)

def test_first_skyima():
    import mosaic
    out_fn = "out.fits"
    mosaic.mosaic_list(["tests/data/isgri_sky_ima_1.fits.gz", "tests/data/isgri_sky_ima_2.fits.gz"], out_fn, pixels="first", mock=False)

    import os
    assert os.path.exists(out_fn)

def test_first_mosaimai_single():
    import mosaic
    out_fn = "out.fits"
    mosaic.mosaic_list(["tests/data/isgri_mosa_ima_1.fits.gz"], out_fn, pixels="first", mock=False)

    import os
    assert os.path.exists(out_fn)

def test_first_mosaima():
    import mosaic
    out_fn = "out.fits"
    mosaic.mosaic_list(["tests/data/isgri_mosa_ima_1.fits.gz", "tests/data/isgri_mosa_ima_2.fits.gz"], out_fn, pixels="first", mock=False)

    import os
    assert os.path.exists(out_fn)

def test_first_mosaima_double():
    import mosaic
    out_p1_fn = "out_p1.fits"
    mosaic.mosaic_list(["tests/data/isgri_mosa_ima_1.fits.gz", "tests/data/isgri_mosa_ima_2.fits.gz"], out_p1_fn, pixels="first", mock=False)
    import os
    assert os.path.exists(out_p1_fn)

    out_fn = "out.fits"
    mosaic.mosaic_list(["tests/data/isgri_mosa_ima_1.fits.gz", "tests/data/isgri_mosa_ima_2.fits.gz", out_p1_fn], out_fn, "first", mock=False)

def test_healpix_skyima():
    import mosaic
    out_fn = "out_p1.fits"
    mosaic.mosaic_list(["tests/data/isgri_sky_ima_1.fits.gz", "tests/data/isgri_sky_ima_2.fits.gz"], out_fn, pixels="healpix", mock=False)

def test_healpix_skyima_mock():
    import mosaic
    out_fn = "out_p1.fits"
    mosaic.mosaic_list(["tests/data/isgri_sky_ima_1.fits.gz", "tests/data/isgri_sky_ima_2.fits.gz"], out_fn, pixels="healpix", mock=True)

def test_first_mosaima_header():
    import mosaic

    out_fn = "out.fits"
    in_fn_1="tests/data/isgri_mosa_ima_1.fits.gz"
    in_fn_2 ="tests/data/isgri_mosa_ima_2.fits.gz"

    i1 = fits.open(in_fn_1)
    h1 = i1[2].header.copy()

    i2=fits.open(in_fn_2)
    h2 = i2[2].header.copy()

    out_header = mosaic.mosaic_list([i1, i2], pixels="first", mock=False)    
    out_header.writeto(out_fn)
    
    assert os.path.exists(out_fn)
    m = fits.open(out_fn)
    h = m[2].header
    
    assert abs(h['ONTIME'] - h1['ONTIME'] - h2['ONTIME']) < 0.01
    assert max(h1['TFIRST'], h2['TFIRST']) == h['TFIRST']
    assert max(h1['TLAST'], h2['TLAST']) == h['TLAST']
    assert max(h1['TELAPSE'], h2['TELAPSE']) <= h['TELAPSE']
    
    
