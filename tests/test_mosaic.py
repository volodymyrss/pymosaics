
def test_healpy():
    pass

def test_first():
    import mosaic
    out_fn = "out.fits"
    mosaic.mosaic_fn_list(["tests/data/isgri_mosa_ima_1.fits.gz", "tests/data/isgri_mosa_ima_2.fits.gz"], out_fn, "first", mock=True)

    import os
    assert os.path.exists(out_fn)
