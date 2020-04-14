
def test_healpy():
    pass

def test_first_mosaima_mock():
    import mosaic
    out_fn = "out.fits"
    mosaic.mosaic_fn_list(["tests/data/isgri_mosa_ima_1.fits.gz", "tests/data/isgri_mosa_ima_2.fits.gz"], out_fn, "first", mock=True)

    import os
    assert os.path.exists(out_fn)

def test_first_skyima_mock():
    import mosaic
    out_fn = "out.fits"
    mosaic.mosaic_fn_list(["tests/data/isgri_sky_ima_1.fits.gz", "tests/data/isgri_sky_ima_2.fits.gz"], out_fn, "first", mock=True)

    import os
    assert os.path.exists(out_fn)

def test_first_skyima():
    import mosaic
    out_fn = "out.fits"
    mosaic.mosaic_fn_list(["tests/data/isgri_sky_ima_1.fits.gz", "tests/data/isgri_sky_ima_2.fits.gz"], out_fn, "first", mock=False)

    import os
    assert os.path.exists(out_fn)

def test_first_mosaima():
    import mosaic
    out_fn = "out.fits"
    mosaic.mosaic_fn_list(["tests/data/isgri_mosa_ima_1.fits.gz", "tests/data/isgri_mosa_ima_2.fits.gz"], out_fn, "first", mock=False)

    import os
    assert os.path.exists(out_fn)

def test_first_mosaima_double():
    import mosaic
    out_p1_fn = "out_p1.fits"
    mosaic.mosaic_fn_list(["tests/data/isgri_mosa_ima_1.fits.gz", "tests/data/isgri_mosa_ima_2.fits.gz"], out_p1_fn, "first", mock=False)
    import os
    assert os.path.exists(out_p1_fn)

    out_fn = "out.fits"
    mosaic.mosaic_fn_list(["tests/data/isgri_mosa_ima_1.fits.gz", "tests/data/isgri_mosa_ima_2.fits.gz", out_p1_fn], out_fn, "first", mock=False)

def test_healpix_mosaima():
    import mosaic
    out_fn = "out_p1.fits"
    mosaic.mosaic_fn_list(["tests/data/isgri_mosa_ima_1.fits.gz", "tests/data/isgri_mosa_ima_2.fits.gz"], out_fn, "healpix", mock=False)

