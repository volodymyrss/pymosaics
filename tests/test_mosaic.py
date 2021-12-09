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
    from astropy.io import fits
    import numpy

    out_fn = "out.fits"
    in_fn_1="tests/data/isgri_mosa_ima_1.fits.gz"
    in_fn_2 ="tests/data/isgri_mosa_ima_2.fits.gz"

    h1=fits.open(in_fn_1)
    ontime1=h1[2].header['ONTIME']
    h2=fits.open(in_fn_2)
    ontime2=h2[2].header['ONTIME']

    out_header = mosaic.mosaic_list([h1,h2], pixels="first", mock=False)
    sum_ontime = out_header.to_hdu_list()[2].header['ONTIME']

    out_header.writeto(out_fn)
    m = fits.open(out_fn)

    import os
    assert os.path.exists(out_fn)
    assert numpy.abs(sum_ontime - ontime1 - ontime2) < 0.01

    assert max(h1[2].header['TFIRST'], h2[2].header['TFIRST']) <= m[2].header['TFIRST']
    assert max(h1[2].header['TLAST'], h2[2].header['TLAST']) <= m[2].header['TLAST']

    assert max(h1[2].header['TELAPSE'], h2[2].header['TELAPSE']) <= m[2].header['TELAPSE']
    #TODO: test more kw
