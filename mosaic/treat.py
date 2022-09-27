# pylint: disable-all

#!/bin/env python

import astropy.io.fits as fits
import numpy as np
import pandas as pd
import subprocess
import os, re

from astropy import wcs
from scipy.optimize import curve_fit

import matplotlib.pylab as plt

# Our goal in this image analysis is to estimate sensitivity in the whole mosaic 
# We also find sources and estimate background component
# see examples on https://apc.u-paris.fr/Downloads/astrog/savchenk/imgb/ISGRI_deep.html

# TODO: find the config
sextractor_share = "/usr/local/share/sextractor"

def render(*args):
    print(*args)

## 2d gaussian

# todo not
from numpy import *
from scipy import optimize

# import crab


def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x, y: height * exp(
        -(((center_x - x) / width_x) ** 2 + ((center_y - y) / width_y) ** 2) / 2
    )


def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments"""
    total = data.sum()
    X, Y = indices(data.shape)
    x = (X * data).sum() / total
    y = (Y * data).sum() / total
    col = data[:, max(0, min(int(y), data.shape[1]-1))]
    width_x = sqrt(abs((arange(col.size) - y) ** 2 * col).sum() / col.sum())
    row = data[max(0, min(int(x), data.shape[0]-1)), :]
    width_y = sqrt(abs((arange(row.size) - x) ** 2 * row).sum() / row.sum())
    height = data.max()
    return height, x, y, width_x, width_y


def fitgaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)

    print("moments:", params)

    errorfunction = lambda p: ravel(gaussian(*p)(*indices(data.shape)) - data)
    # errorfunction = lambda p: ravel(gaussian(*p)(*indices(data.shape)) -
    # data)
    p, success = optimize.leastsq(errorfunction, params)
    print("optimized", p, success)
    return p


##


class SExtractor:
    back_size = 4
    threshold = 3

    def __init__(self, hostdir=None, back_size=4, threshold=3):
        self.hostdir = hostdir
        os.makedirs(hostdir, exist_ok=True)

    def set_mosaic(self, fn):
        self.mosaic_fn = fn

    def run(self):
        cwd = os.getcwd()
        os.chdir(self.hostdir)
        try:
            self._run()
            os.chdir(cwd)
        except:
            os.chdir(cwd)
            raise
            

    def _run(self):
        config = open(sextractor_share + "/default.sex").read()
        open("default.param", "w").write("\n".join([
                                        "X_IMAGE",
                                        "Y_IMAGE",
                                        "MAG_BEST",
                                        "ALPHA_SKY",
                                        "DELTA_SKY",
                                        "FLUX_ISO"]))

        #conv = open(sextractor_share + "/default.conv").read()
        open("default.conv", "w").write(
            "CONV NORM\n"
            "0 0 0\n"
            "0 1 0\n"
            "0 0 0\n")
        

        def set_key(config, key, value):
            return re.sub("(?<=" + key + "[ \t])(.*?)(?=\n)", value + " #", config)

        check_images = [
            ("BACKGROUND_RMS", "background_rms.fits"),
            ("BACKGROUND", "background.fits"),
            ("FILTERED", "filtered.fits"),
            ("-OBJECTS", "noobjects.fits"),
            ("-BACKGROUND", "nobackground.fits"),
            ("APERTURES", "apertures.fits"),
            ("SEGMENTATION", "segmentation.fits"),
        ]

        config = set_key(config, "CHECKIMAGE_TYPE", " ".join(list(zip(*check_images))[0]))
        config = set_key(config, "CHECKIMAGE_NAME", " ".join(list(zip(*check_images))[1]))
        config = set_key(config, "BACK_SIZE", str(self.back_size))
        config = set_key(config, "BACK_FILTERSIZE", "1")
        config = set_key(config, "DETECT_THRESH", str(self.threshold))
        config = set_key(config, "ANALYSIS_THRESH", str(self.threshold))
        open("default.sex", "w").write(config)

        # find binary
        subprocess.check_call(["sex",  "../" + self.mosaic_fn])

    def background_rms_ext(self):
        return fits.open(self.hostdir + "/background_rms.fits")[0]

    def background_ext(self):
        return fits.open(self.hostdir + "/background.fits")[0]

    def filtered_ext(self):
        return fits.open(self.hostdir + "/filtered.fits")[0]

    def noobjects_ext(self):
        return fits.open(self.hostdir + "/noobjects.fits")[0]




class DistributionAnalysis:
    def __init__(self, pixels):
        self.set_pixels(pixels)

    def set_pixels(self, pixels):
        self.pixels = pixels

    def show_all(self):
        pixels = self.pixels
        h = plot.p.hist(
            pixels[~isnan(pixels)].flatten(),
            linspace(-6, 20, 1000),
            lw=0,
            alpha=0.5,
            log=True,
            normed=True,
        )
        hist, bin_edges = h[0:2]
        bin_centres = (bin_edges[:-1] + bin_edges[1:]) / 2

        bin_centres, hist_fit = self.core_fit

        plot.p.plot(bin_centres, hist_fit)

        plot.p.xlim([-10, 10])
        plot.p.ylim([1e-5, 1])
        plot.p.grid()
        plot.plot("sig_histogram.png")

    def calibrate(self):
        pixels = self.pixels
        h = plot.p.hist(
            pixels[~isnan(pixels)].flatten(), linspace(-10, 10, 30), lw=0, alpha=0.5
        )
        hist, bin_edges = h[0:2]
        bin_centres = (bin_edges[:-1] + bin_edges[1:]) / 2

        # posp=hist[:]
        # posp[bin_centres<=0]=0

        symneg = hist[::-1].copy()
        symneg[bin_centres < 0] = 0
        symneg = symneg + symneg[::-1]
        symneg[symneg.shape[0] / 2] /= 2

        plot.p.plot(bin_centres, hist)
        plot.p.plot(bin_centres, symneg)
        plot.p.errorbar(bin_centres, hist - symneg, sqrt(symneg))

        rec_thresh = -(bin_centres[symneg > 0].min())

        print("recommended threshold for 1 noise source", rec_thresh)

        plot.p.semilogy()

        plot.p.grid()
        plot.plot("sig_histogram.png")

        self.calibration = None

        import scipy.stats

        print(render("{RED}kurtosis:{/}"), scipy.stats.kurtosis(pixels))

    def fit_core(self):
        pixels = self.pixels
        h = plot.p.hist(
            pixels[~isnan(pixels)].flatten(),
            linspace(-6, 6, 1000),
            lw=0,
            alpha=0.5,
            log=True,
            normed=True,
        )
        hist, bin_edges = h[0:2]
        bin_centres = (bin_edges[:-1] + bin_edges[1:]) / 2

        def gauss(x, *p):
            A, mu, sigma = p
            return A * exp(-((x - mu) ** 2) / (2.0 * sigma ** 2)) / sqrt(2 * pi * sigma)

        def lorentz(x, *p):
            A, mu, gamma = p
            return A * gamma ** 2 / (gamma ** 2 + (x - mu) ** 2)

        def profile(x, *p):
            return gauss(x, *p[:3]) + gauss(x, *p[3:])

        p0 = [1.0, 0.0, 1.0, 1.0, 0.0, 1.0]
        # p0 = [1., 0., 1.]
        coeff, var_matrix = curve_fit(profile, bin_centres, hist, p0=p0)
        hist_fit = profile(bin_centres, *coeff)

        self.core_fit_function = lambda x: profile(x, *coeff)
        self.core_fit = bin_centres, hist_fit

        print("fit result", coeff)

        plot.p.plot(bin_centres, hist_fit)

        plot.p.xlim([-10, 10])
        plot.p.ylim([1e-5, 1])
        plot.p.grid()
        plot.plot("sig_histogram_core.png")


class ImageAnalysis:
    back_size = 30
    # Size of background region estimate
    exposure_fraction_cut = 100
    #It cuts at 1/100 of max exposure
    threshold = 4
    # S/N threshold
    cached = True
    source_analysis = False
    # get source parameters
    hostdir = "./"
    #output dir ?
    version = "v1"
    i_band = 0
    #convenience variable
    source_radius = 5
    annulus_radius = 10
    save_plots = False
    additional_sources = None
    reference_catalog = None


    def __init__(self, in_fn='none', int_out_prefix='', back_size=30, exposure_fraction_cut=100,
                 source_analysis=False, hostdir="./", threshold=4, source_radius=5, annulus_radius=10,
                 save_plots=False):
        self.raw_mosaic_fn = in_fn
        self.out_prefix = int_out_prefix
        self.hostdir = hostdir
        self.back_size = back_size
        self.exposure_fraction_cut = exposure_fraction_cut
        self.source_analysis = source_analysis
        self.threshold = threshold
        self.source_radius = source_radius
        self.annulus_radius = annulus_radius
        self.save_plots = save_plots

    def fullpath(self, x):
        return x

    def raw_erange(self):
        h = self.raw_intensity_ext().header
        if 'E_MIN' in h.keys() and "E_MAX" in h.keys():
            return h["E_MIN"], h["E_MAX"]
        else:
            return -1, -1

    def set_mosaic_fn(self, fn):
        self.raw_mosaic_fn = fn

    def get_mosaic_fn(self):
        if hasattr(self, "input_image"):
            if hasattr(self.input_image, "skyima"):
                return self.input_image.skyima.get_path()

        return self.raw_mosaic_fn

    def e_range(self):
        return self.raw_erange()

    def get_n_ebands(self):
        return self.get_total_ebands()

    def get_total_ebands(self):
        gt = fits.open(self.get_mosaic_fn())
        # print("group:", gt)
        # print(len(gt))
        # print("found bands:", (len(gt) - 1) / 4)
        return int((len(gt) - 1) / 4)

    def get_imatype(self, myhead, default=None):
        imatype = myhead.get('IMATYPE', '')
        extname = myhead.get('EXTNAME', '')

        # Priority on IMATYPE !        
        if imatype != '':
            return imatype
        elif extname != '':
            return extname
        else:
            if default is None:
                raise KeyError('Both EXTNAME and IMATYPE are empty, impossible to get image type')
            else:
                return default

    def get_extension_by_type(self, exttype: str):
        with fits.open(self.get_mosaic_fn()) as ff:
            h_ret = None
            for hh in ff:
                if self.get_imatype(hh.header, '') == exttype:
                    if h_ret is not None:
                        raise NotImplementedError(f"found two extensions for requested type {exttype}")

                    h_ret = hh.copy()

        # Cannot close for way of dealing with input of this class
        #ff.close()
        if h_ret is not None:
            return h_ret
        else:
            return None

    def raw_exposure_ext(self):
        ret = self.get_extension_by_type('EXPOSURE')
        if ret is None:
            raise RuntimeError('Cannot find EXPOSURE in %s' % (self.get_mosaic_fn()))
        return ret

    def raw_intensity_ext(self):
        #IBIS
        ret = self.get_extension_by_type('INTENSITY')
        if ret is None:
            #JEM-X
            ret = self.get_extension_by_type('RECONSTRUCTED')
            if ret is None:
                raise RuntimeError('Cannot find INTENSITY/Reconstructed in %s' % (self.get_mosaic_fn()))
        return ret

    def raw_significance_ext(self):
        ret = self.get_extension_by_type('SIGNIFICANCE')
        if ret is None:
            raise RuntimeError('Cannot find SIGNIFICANCE in %s' % (self.get_mosaic_fn()))
        return ret

    def raw_variance_ext(self):
        ret = self.get_extension_by_type('VARIANCE')
        if ret is None:
            raise RuntimeError('Cannot find VARIANCE in %s' % (self.get_mosaic_fn()))
        return ret


    def inspect_raw(self):
        try:
            print("exposure", self.raw_exposure_ext().header["EXTNAME"])
            print("intensity", self.raw_intensity_ext().header["EXTNAME"])
            print("significance", self.raw_significance_ext().header["EXTNAME"])
            print("variance", self.raw_variance_ext().header["EXTNAME"])
        except:
            print("mosaic without extname?")
        try:
            print("exposure", self.raw_exposure_ext().header["IMATYPE"])
            print("intensity", self.raw_intensity_ext().header["IMATYPE"])
            print("significance", self.raw_significance_ext().header["IMATYPE"])
            print("variance", self.raw_variance_ext().header["IMATYPE"])
        except:
            print("mosaic without imatype?")

    def get_shortname(self):
        return "ImageAnalysis"

    def exposure_ext(self):
        return self.raw_exposure_ext()

    def intensity_ext(self):
        return self.raw_intensity_ext()

    def variance_ext(self):
        return self.raw_variance_ext()

    def main(self):
        print("About to analyze %s energy bands" % self.get_n_ebands())
        for self.i_band in range(self.get_n_ebands()):
            self.tag = "%.5lg_%.5lg" % self.raw_erange()
            #print(self.tag)
            self.mask_mosaic()
            self.extract_sources()
            if self.source_analysis:
                self.read_sources()
            self.estimate_sensitivity()
            # self.pixel_distribution()
        self.dump_statistics()

    def pixel_distribution(self):
        pixels = fits.open("significance.fits")[0].data
        D = DistributionAnalysis(pixels)
        D.calibrate()
        self.pixel_calibration = D.calibration
        # D.fit_core()
        # D.show_all()

    def estimate_sensitivity(self):
        print("opening mosaic...")
        exposure = self.exposure_ext().data
        variance = self.variance_ext().data
        image = self.sextractor.filtered_ext().data
        sigimage = self.sextractor2.filtered_ext().data

        print("opening rms...")

        rms = self.sextractor.background_rms_ext().data

        bkg = self.sextractor.background_ext().data

        print(rms)

        def avg_nonan(a):
            return average(a[where(~isnan(a))])

        # sensi=rms/0.015*1e-11

        e1, e2 = self.raw_erange()

        # sensi = crab.Crab().to_mcrab(rms, e1, e2)
        sensi = rms

        # crab_here = crab.Crab().counts_in(e1, e2)
        # cts2ecs = lambda x: crab.Crab().to_ecs(x, e1, e2)
        # cts2mcrab = lambda x: crab.Crab().to_mcrab(x, e1, e2)
        # mcrab2ecs = lambda x: cts2ecs(x * crab_here / 1000.0)
        # print("Crab is:", crab_here, "cts", cts2ecs(crab_here), "ecs in", e1, e2)
        cts2ecs = lambda x:x
        cts2mcrab = lambda x:x
        mcrab2ecs = lambda x:x

        print(render("{RED}enery band{/}"), e1, e2)
        print(
            render("{RED}avg_nonan{/}                             "),
            avg_nonan(sensi),
            render("{YEL}mcrab{/}"),
            mcrab2ecs(avg_nonan(sensi)),
            "ecs",
            avg_nonan(exposure) / 1e6,
            "Ms",
        )

        print(
            render("{RED}minimal sensitivity{/}                   "),
            sensi.min(),
            render("{YEL}mcrab{/}"),
            mcrab2ecs(sensi.min()),
            "ecs",
            exposure.flatten()[sensi.flatten().argmin()] / 1e6,
            "Ms",
        )

        print(
            render("{RED}exposure-corrected avg_nonan{/}          "),
            avg_nonan(sensi * sqrt(exposure / 1e6)),
            render("{YEL}mcrab-Ms{/}"),
            mcrab2ecs(avg_nonan(sensi * sqrt(exposure / 1e6))),
            "ecs",
            avg_nonan(exposure) / 1e6,
            "Ms",
        )

        print(
            render("{RED}exposure-corrected minimal sensitivity{/}"),
            sensi.min() * 1e-3 * sqrt(exposure.flatten()[sensi.flatten().argmin()]),
            render("{YEL}mcrab-Ms{/}"),
            mcrab2ecs(
                sensi.min() * 1e-3 * sqrt(exposure.flatten()[sensi.flatten().argmin()])
            ),
            "ecs",
            exposure.flatten()[sensi.flatten().argmin()] / 1e6,
            "Ms",
        )

        bmin, bmax = bkg.min(), bkg.max()
        print(
            render("{RED}smooth background range{/}                   "),
            cts2mcrab(bmin),
            cts2mcrab(bmax),
            render("{YEL}mcrab{/}"),
            cts2ecs(bmin),
            cts2ecs(bmax),
            "ecs",
        )

        statdict = dict(
            sensi_avg_nonan=avg_nonan(sensi),
            sensi_min=sensi.min(),
            rms_min=rms.min(),
            rms_avg_nonan=avg_nonan(rms),
            exposure_at_sensi_min=exposure.flatten()[rms.flatten().argmin()],
            exposure_avg_nonan=avg_nonan(exposure),
            sensi_min_exposure_corrected_Ms=sensi.min()
            * sqrt(exposure.flatten()[sensi.flatten().argmin()])
            * 1e-3,
            sensi_avg_nonan_exposure_corrected_Ms=avg_nonan(
                sensi * sqrt(exposure / 1e6)
            ),
            exposyre_max=exposure.max(),
            smooth_b_min=cts2mcrab(bmin),
            smooth_b_max=cts2mcrab(bmax),
            significance_max=sigimage.max(),
            rate_max=cts2mcrab(image.flatten()[sigimage.argmax()]),
        )

        self.save_statistics(statdict)

        # tosave=[avg_nonan(sensi),
        #        avg_nonan(exposure),
        #        sensi.min(),
        #        exposure.flatten()[rms.flatten().argmin()]]
        # open("statistic_raw.txt","w").write()
        # savetxt("statistic_raw.txt")

        self.save_as_fits(sensi, "sensitivity%s.fits" % self.tag)
        self.save_as_fits(
            sensi * sqrt(exposure) * 1e-3, "relative_sensitivity%s.fits" % self.tag
        )
        self.save_as_fits(rms ** 2 / variance, "excess_rms%s.fits" % self.tag)

        rms = self.sextractor2.background_rms_ext().data
        significance = self.sextractor2.filtered_ext().data / rms
        significance[significance < -1e10] = NaN
        self.save_as_fits(significance, "significance%s.fits" % self.tag)

        significance = self.sextractor2.noobjects_ext().data / rms
        significance[significance < -1e10] = NaN
        self.save_as_fits(significance, "significance_noobjects%s.fits" % self.tag)

    def save_statistics(self, statdict):
        e1, e2 = self.raw_erange()
        fnstat = "statistic_%.5lg_%.5lg.txt" % (e1, e2)
        statdict["E_MIN"] = e1
        statdict["E_MAX"] = e2
        open(fnstat, "w").write(repr(statdict))
        if not hasattr(self, "statistics"):
            self.statistics = []
        self.statistics.append(statdict)

    def dump_statistics(self):
        e1, e2 = self.raw_erange()
        keys = sorted(self.statistics[0].keys())
        fn = "statistics.txt"
        f = open(fn, "w")
        f.write(" ".join(keys) + "\n")
        for s in self.statistics:
            f.write(" ".join(["%.5lg" % s[k] for k in keys]) + "\n")
        f.close()
        # self.statfile = da.DataFile(fn)

    def save_as_fits(self, data, fn):
        h = fits.PrimaryHDU(data)
        wcsh = self.raw_wcs().to_header()
        h.header.extend(wcsh)
        h.writeto(self.hostdir+self.out_prefix+fn, overwrite=True)

    def extract_sources(self):
        self.sextractor = SExtractor("sextractor")
        self.sextractor.back_size = self.back_size

        print(self.mosaic_fn + "[1]")
        print(self.fullpath(self.mosaic_fn) + "[1]")

        self.sextractor.set_mosaic(self.fullpath(self.mosaic_fn + "[1]"))
        self.sextractor.threshold = self.threshold
        self.sextractor.run()

        rms = self.sextractor.background_rms_ext().data
        significance = self.sextractor.filtered_ext().data / rms
        significance[significance < -1e10] = NaN
        self.save_as_fits(significance, "significance_first_pass.fits")

        self.sextractor2 = SExtractor(self.hostdir + "/sextractor2")
        self.sextractor2.back_size = self.back_size
        self.sextractor2.threshold = self.threshold
        self.sextractor2.set_mosaic(self.fullpath(self.out_prefix+"significance_first_pass.fits"))
        self.sextractor2.run()

    def get_additional_sources(self):
        ra = self.additional_sources['ra']
        dec = self.additional_sources['dec']
        #print(ra, dec)
        x, y = self.raw_wcs().world_to_pixel_values(ra, dec)
        #print(x, y)
        cc = zeros_like(self.sources, shape=(len(x),))
        cc['X_IMAGE'] = x
        cc['Y_IMAGE'] = y
        cc['ALPHA_SKY'] = ra
        cc['DELTA_SKY'] = dec
        cc['FLUX_ISO'] = ones(len(x))
        # for a, b in zip(x, y):
        #     cc = array(dtype=self.sources.dtype)
        #     cc[''] = a
        #     cc[1] = b
        #     # print(array([[a], [b], [0], [0], [0]], dtype=self.sources.dtype))
        self.sources = append(self.sources, cc)
        #print(self.sources)

    @staticmethod
    def match_sources(src, isdc_sources, threshold=5. / 60.):
        from astropy import units as u
        from astropy.coordinates import SkyCoord
        c = SkyCoord(ra=src['ra'] * u.degree, dec=src['dec'] * u.degree)
        catalog = SkyCoord(ra=isdc_sources['RA_OBJ'] * u.degree, dec=isdc_sources['DEC_OBJ'] * u.degree)
        idx, d2d, d3d = c.match_to_catalog_sky(catalog)
        names = []
        i = 1
        for ii, dd in zip(idx, d2d.deg):
            if dd < threshold:
                names.append(isdc_sources['NAME'][ii])
            else:
                names.append('NEW_%02d' % i)
                i += 1
        return names

    def read_sources(self):
        cat = open(self.sextractor.hostdir + "/test.cat").read()
        cd = re.findall("#(.*?)\n", cat)

        try:
            keys, keypositions = list(
                zip(*[reversed(re.search(" *(\d+) (.*?) ", c).groups()) for c in cd])
            )
        except Exception as e:
            print("problem with keys: ", e)

        # print("Reading " + self.sextractor2.hostdir + "/test.cat")
        sources = genfromtxt(self.sextractor2.hostdir + "/test.cat", names=keys)
        # print(keys)
        # print(sources)
        self.sources = sources

        if self.additional_sources is not None:
            self.get_additional_sources()


        #######################################################################
        #                             get sources fluxes
        #######################################################################

        image = self.sextractor.filtered_ext().data
        rmsimage = self.sextractor.background_rms_ext().data
        sigimage = self.sextractor2.filtered_ext().data
        exposure = self.raw_exposure_ext().data
        ix, iy = meshgrid(arange(image.shape[1]), arange(image.shape[0]))

        print(ix.shape, iy.shape, image.shape)

        results = []
        results_extra = []

        image[image < -1e29] = NaN

        fns = ""

        flag = sources["X_IMAGE"] * 0

        #print(sources, sources.shape, sources.dtype)

        try:
            if len(sources) == 0:
                print("no sources!")
                return
        except Exception as e:
            print("problem with sources ", e)
            return

        for i, (x, y, ra, dec) in enumerate(
            zip(
                self.sources["X_IMAGE"],
                self.sources["Y_IMAGE"],
                self.sources["ALPHA_SKY"],
                self.sources["DELTA_SKY"],
            )
        ):
            r = (ix - x) ** 2 + (iy - y) ** 2
            source = r < self.source_radius ** 2
            ring = logical_and(r < self.annulus_radius ** 2, r > self.source_radius ** 2)

            print(ring.shape, source.shape, r.shape)

            print(
                "source", where(source)[0].shape, self.sources["FLUX_ISO"][i], x, y, ra, dec
            )
            print("ring", where(ring)[0].shape)

            if any(isnan(image[logical_or(source, ring)])):
                print("there is nan in the source region! skipping")
                flag[i] = 1
                continue

            expo = average(exposure[source])

            good = image > 1e-20
            source = logical_and(source, good)
            ring = logical_and(ring, good)

            if image[source].flatten().shape[0] == 0:
                print("empty image for this one!")
                continue

            peakcounts = image[source].max()
            peaksig = sigimage[source].max()
            sumcounts = image[source].sum() - average(image[ring]) * source.sum()

            rms = average(rmsimage[source])

            # print
            # print "position:",x,y,ra,dec
            # print "peak,sum,rms",peakcounts,sumcounts,rms
            # print "peak sig",peaksig
            # print "significance:",peakcounts/rms,sumcounts/rms
            # print "exposure:",expo

            results.append([x, y, ra, dec, peakcounts, rms, expo])

            try:
                #print('Pippo')
                image_shape = image.shape
                fr = fitgaussian(image[max(0,int(y) - 7): min(int(y) + 7, image_shape[0]-1),
                                        max(0,int(x) - 7): min(int(x) + 7,image_shape[1]-1)])
                #print('Pippo')
                #print(fr)
                wx, wy = fr[-2:]
                dx, dy = fr[1:3]
                # x+=dx-6.5
                # y+=dy-6.5
                ecc = (1 - (min(wx, wy) / max(wx, wy)) ** 2) ** 0.5

                results_extra.append(
                    [x, y, ra, dec, peakcounts, rms, expo] + list(fr) + [ecc]
                )
                if self.save_plots:
                    plt.clf()
                    plt.imshow(
                        sigimage[int(y) - 10 : int(y) + 10, int(x) - 10 : int(x) + 10], interpolation="none"
                    ).set_clim(1, 5)
                    plt.colorbar()
                    plt.title("rms: %.5lg" % rms)
                    plt.savefig("source_%.5lg.png" % peaksig)
                    fns += "source_%.5lg.png " % peaksig
            except Exception as e:
                print("problem with fitting: ", e)
                raise RuntimeError(e)


        #save(self.hostdir + "/sources.npy", sources)
        if self.reference_catalog is not None:
            #print("MATCHING SOURCES", type(results), len(results))
            self.source_names = ImageAnalysis.match_sources({'ra': [xx[2] for xx in results],
                                                             'dec': [xx[3] for xx in results]},
                                                            self.reference_catalog)

        if hasattr(self, 'source_names'):
            regionfile_extra = "\n".join(
                [
                    "circle(%.5lg,%.5lg,3) # text={%s}"
                    % (x[0], x[1], self.source_names[i])
                    for i, x in enumerate(results_extra)
                ]
            )
            regionfile = "\n".join(
                [
                    "circle(%.5lg,%.5lg,3) # text={%s}"
                    % (x[0], x[1], self.source_names[i])
                    for i, x in enumerate(results_extra)
                ]
            )
        # write region file
        else:
            regionfile_extra = "\n".join(
                [
                    "circle(%.5lg,%.5lg,3) # text={%.5lg,e=%.5lg}"
                    % (x[0], x[1], x[4] / x[5], x[-1])
                    for x in results_extra
                ]
            )
            regionfile = "\n".join(
                [
                    "circle(%.5lg,%.5lg,3) # text={%.5lg,e=%.5lg}"
                    % (x[0], x[1], x[4] / x[5], x[-1])
                    for x in results
                ]
            )
        # regionfile="\n".join(["circle(%.5lg,%.5lg,3) # text={%.5lg,e=%.5lg}"%(x,y,peaksig,ecc) for x,y in zip(sources[flag==0]["X_IMAGE"],sources[flag==0]["Y_IMAGE"])])
        regionfile = (
            """
# Region file format: DS9 version 4.1
# Filename: significance.fits
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
image
"""
            + regionfile
        )
        regionfile_extra = (
            """
# Region file format: DS9 version 4.1
# Filename: significance.fits
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
image
"""
            + regionfile_extra
        )
        open(self.hostdir+self.out_prefix+"source.reg", "w").write(regionfile)
        open(self.hostdir+self.out_prefix+"source_extra.reg", "w").write(regionfile_extra)

        self.source_results = pd.DataFrame(data=results,
            columns=['x', 'y', 'ra', 'dec', 'peakcounts', 'rms', 'expo'])

        self.source_results_extra = pd.DataFrame(data=results_extra,
            columns=['x', 'y', 'ra', 'dec', 'peakcounts', 'rms', 'expo', 'height', 'x_gauss', 'y_gauss',
                     'width_x', 'width_y', 'ecc'])

        if hasattr(self, 'source_names'):
            #print("ADDING Names to Pandas")
            self.source_results['name'] = self.source_names
            self.source_results_extra['name'] = self.source_names


        self.source_results.to_csv(self.hostdir+self.out_prefix+'source_results.csv')
        self.source_results_extra.to_csv(self.hostdir+self.out_prefix+'source_results_extra.csv')

        savetxt(self.hostdir+self.out_prefix+"source_results.txt", array(results))
        savetxt(self.hostdir+self.out_prefix+"source_results_extra.txt", array(results_extra))
        if self.save_plots:
            os.system("montage -geometry 300x300 " + fns + " sources.png")


    def mask_mosaic(self):
        self.inspect_raw()

        wcs = self.raw_wcs()

        #print("wcs:", wcs.to_header_string())

        exposure = self.raw_exposure_ext().data

        good_exposure = exposure[~isnan(exposure) & ~isinf(exposure) & (exposure != 0)]

        print("maximum exposure", good_exposure.max())
        print("average exposure", average(good_exposure))

        cut = exposure < good_exposure.max() / self.exposure_fraction_cut

        print("exposure cut below: ", exposure[~cut].min())

        hdu_list = []

        print("energy range:", self.raw_erange())

        sig_image = self.raw_significance_ext().data
        flux_image = self.raw_intensity_ext().data
        var_image = self.raw_variance_ext().data

        sig_image[cut] = NaN
        flux_image[cut] = NaN
        var_image[cut] = NaN

        print("band", self.i_band)

        # if self.i_band==0:
        hdu_list.append(fits.PrimaryHDU(sig_image))
        # else:
        #    hdu_list.append(fits.ImageHDU(sig_image))

        hdu_list.append(fits.ImageHDU(flux_image))
        hdu_list.append(fits.ImageHDU(var_image))

        hdu_list.append(fits.ImageHDU(exposure))

        wcsh = self.raw_wcs().to_header()
        e1, e2 = self.raw_erange()
        for h in hdu_list:
            h.header.extend(wcsh)
            h.header["E_MIN"] = e1
            h.header["E_MAX"] = e2
            print("updating header:", e1, e2)

        fn = self.hostdir+self.out_prefix+"masked_mosaic_%s.fits" % (self.tag)

        fits.HDUList(hdu_list).writeto(fn, overwrite=True)

        self.mosaic_fn = fn

    def raw_wcs(self):
        return wcs.WCS(self.raw_exposure_ext().header)


class VarmosaicImageAnalysis(ImageAnalysis):
    pass
    # def get_n_ebands(self):
    #     # for i in range((len(mosaic)-1)/4):
    #     return 1
    #
    # def set_raw_mosaic_fn(self, fn):
    #     self.raw_mosaic_fn = self.fullpath(fn)
    #
    # def raw_erange(self, i):
    #     h = self.raw_intensity_ext(i).header
    #     return h["E_MIN"], h["E_MAX"]
    #
    # def raw_exposure_ext(self):
    #     return fits.open(self.raw_mosaic_fn)[-1]
    #
    # def raw_significance_ext(self, i):
    #     return fits.open(self.raw_mosaic_fn)[i * 3]
    #
    # def raw_intensity_ext(self, i):
    #     return fits.open(self.raw_mosaic_fn)[1 + i * 3]
    #
    # def raw_variance_ext(self, i):
    #     return fits.open(self.raw_mosaic_fn)[2 + i * 3]


class SimpleImageAnalysis:
    back_size = 30
    exposure_fraction_cut = 2
    threshold = 4

    cached = True

    source_analysis = False

    hostdir = "./"

    version = "v1"

    out_prefix = ''

    def fullpath(self, x):
        return x

    def inspect_raw(self):
        try:
            print("exposure", self.raw_exposure_ext().header["EXTNAME"])
            print("intensity", self.raw_intensity_ext().header["EXTNAME"])
            print("significance", self.raw_significance_ext().header["EXTNAME"])
            print("variance", self.raw_variance_ext().header["EXTNAME"])
        except:
            print("mosaic without extname?")
        try:
            print("exposure", self.raw_exposure_ext().header["IMATYPE"])
            print("intensity", self.raw_intensity_ext().header["IMATYPE"])
            print("significance", self.raw_significance_ext().header["IMATYPE"])
            print("variance", self.raw_variance_ext().header["IMATYPE"])
        except:
            print("mosaic without imatype?")

    def get_shortname(self):
        return "ImageAnalysis"

    def exposure_ext(self):
        raise

    def intensity_ext(self):
        raise

    def main(self):
        for self.i_band in range(self.get_n_ebands()):
            self.tag = "_%.5lg_%.5lg" % self.raw_erange()
            self.estimate_sensitivity()
            # self.pixel_distribution()
        self.dump_statistics()

    def estimate_sensitivity(self):
        print("opening mosaic...")
        exposure = self.raw_exposure_ext().data

        variance = self.raw_variance_ext().data
        significance = self.raw_significance_ext().data
        intensity = self.raw_intensity_ext().data
        nce = self.raw_significance_ext().data

        mask = exposure > exposure.max() / self.exposure_fraction_cut

        from scipy.ndimage.filters import gaussian_filter as gf

        rms = (gf(intensity ** 2, 30) - gf(intensity, 30) ** 2) ** 0.5

        fits.PrimaryHDU(rms).writeto(self.hostdir+self.out_prefix+"rms_%i.fits" % self.i_band, overwrite=True)

        total_rms = std(significance[mask])

        def avg_nonan(a):
            return average(a[where(~isnan(a))])

        #e1, e2 = self.erange(self.i_band)

        statdict = dict(
            rms=total_rms,
            rms_min=rms.min(),
            variance_min=variance[variance > 0].min() ** 0.5,
            exposure_max=exposure.max(),
            significance_max=significance.max(),
        )

        self.save_statistics(statdict)

    def save_statistics(self, statdict):
        e1, e2 = self.raw_erange()
        fnstat = "statistic_%.5lg_%.5lg.txt" % (e1, e2)
        statdict["E_MIN"] = e1
        statdict["E_MAX"] = e2
        open(fnstat, "w").write(repr(statdict))
        if not hasattr(self, "statistics"):
            self.statistics = []
        self.statistics.append(statdict)

    def dump_statistics(self):
        e1, e2 = self.raw_erange()
        keys = sorted(self.statistics[0].keys())
        fn = "statistics.txt"
        f = open(fn, "w")
        f.write(" ".join(keys) + "\n")
        for s in self.statistics:
            f.write(" ".join(["%.5lg" % s[k] for k in keys]) + "\n")
        f.close()
        # self.statfile = da.DataFile(fn)

    def raw_wcs(self):
        return wcs.WCS(self.raw_exposure_ext().header)

    def erange(self):
        h = self.raw_intensity_ext().header
        return h["E_MIN"], h["E_MAX"]


class SimpleSingleScWImageAnalysis(SimpleImageAnalysis):
    def get_n_ebands(self):
        return self.get_total_ebands()

    def get_total_ebands(self):
        gt = fits.open(self.get_raw_mosaic_fn())[1].data
        print("group:", gt)
        print("found bands:", (len(gt)) / 5)
        return (len(gt)) / 5

    def set_raw_mosaic_fn(self, fn):
        self.raw_mosaic_fn = fn

    def get_raw_mosaic_fn(self):
        if hasattr(self, "input_image"):
            if hasattr(self.input_image, "skyima"):
                return self.input_image.skyima.get_path()

        return self.raw_mosaic_fn

    # def raw_erange(self):
    #     h = self.raw_intensity_ext().header
    #     return h["E_MIN"], h["E_MAX"]
    #
    # def raw_exposure_ext(self):
    #     # return fits.open(self.raw_mosaic_fn)[6+i*4]
    #     return fits.open(self.get_raw_mosaic_fn())[-1]
    #
    # def raw_intensity_ext(self):
    #     return fits.open(self.get_raw_mosaic_fn())[2 + self.i_band * 5]
    #
    # def raw_significance_ext(self):
    #     return fits.open(self.get_raw_mosaic_fn())[4 + self.i_band * 5]
    #
    # def raw_variance_ext(self):
    #     return fits.open(self.get_raw_mosaic_fn())[3 + self.i_band * 5]


class OSAMosaicImageAnalysis(ImageAnalysis):


    def set_raw_mosaic_fn(self, fn):
        self.raw_mosaic_fn = fn

    def get_mosaic_fn(self):
        if hasattr(self, "input_image"):
            if hasattr(self.input_image, "skyima"):
                return self.input_image.skyima.get_path()

        return self.raw_mosaic_fn

    # def raw_erange(self):
    #     h = self.raw_intensity_ext().header
    #     if 'E_MIN' in h.keys() and "E_MAX" in h.keys():
    #         return h["E_MIN"], h["E_MAX"]
    #     else:
    #         return -1, -1
    #
    # def raw_exposure_ext(self):
    #     # return fits.open(self.raw_mosaic_fn)[6+i*4]
    #     return fits.open(self.get_raw_mosaic_fn())[-1]
    #
    # def raw_intensity_ext(self):
    #     return fits.open(self.get_raw_mosaic_fn())[2 + self.i_band * 4]
    #
    # def raw_significance_ext(self):
    #     return fits.open(self.get_raw_mosaic_fn())[4 + self.i_band * 4]
    #
    # def raw_variance_ext(self):
    #     return fits.open(self.get_raw_mosaic_fn())[3 + self.i_band * 4]
