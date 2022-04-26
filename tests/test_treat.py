def test_treat_mosaima():
    import mosaic.treat
    out_fn = "out_p1.fits"
    
    oia = mosaic.treat.OSAMosaicImageAnalysis()
    oia.raw_mosaic_fn = "tests/data/isgri_mosa_ima_1.fits.gz"
    oia.source_analysis = True
    oia.exposure_fraction_cut = 100

    oia.main()
