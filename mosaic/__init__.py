import click

from typing import List

from astropy import wcs as pywcs
from astropy.io import fits

import healpy

import numpy as np

def mosaic_fn_list(in_fn: List[str], out_fn: str, out_format: str, mock=False):
    print("input: %s output: %s"%( ", ".join(in_fn), out_fn))

    if out_format == "first":
        m = FITsMosaic(mock)
        for fn in in_fn:
            m.add(fn)

        m.writeto(out_fn)
    else:
        raise RuntimeError


class Mosaic:
    def writeto(self, fn: str):
        pass

    def add(self, fn: str):
        pass

    @staticmethod
    def m(a, s_a, b, s_b):
        c=dict()
        c['flux'] = s_a(a['flux']) + s_b(b['flux']) 
        c['var'] = ( 1 / s_a(a['var']) + 1 / s_b(b['var']) ) ** -1
        c['ex'] = s_a(a['ex']) + s_b(b['ex'])
        c['n'] = s_a(a['n']) + s_b(b['n'])
        return c

    @staticmethod
    def m_wm(a, s_a, b, s_b):
        c=dict()
        c['flux'] = s_a(a['flux']) / s_a(a['var']) + s_b(b['flux']) / s_b(b['var'])
        c['var'] = ( 1 / s_a(a['var']) + 1 / s_b(b['var']) ) ** -1
        c['flux'] = c['flux'] / c['var']
        c['ex'] = s_a(a['ex']) + s_b(b['ex'])
        c['n'] = s_a(a['n']) + s_b(b['n'])
        return c

class HealpixMosaic(Mosaic):
    def __init__(self, nsides=64, mock=False):
        self.nsides = nsides
        self.mosaic = dict(
                flux = np.zeros(healpy.nsides2npix(nsides)),
                var = np.zeros(healpy.nsides2npix(nsides)),
                ex = np.zeros(healpy.nsides2npix(nsides)),
                n = np.zeros(healpy.nsides2npix(nsides)),
            )

    def writeto(self, fn: str):
        healpy.write_map(fn, self.mosaic)

class FITsMosaic(Mosaic):
    def __init__(self, mock=False):
        self.mosaic = None
        self.mock = mock

    def add(self, fn):
        print("adding", fn)

        f = fits.open(fn)

        def imatype_selector(_f, n):
            for e in _f:
                if e.header.get('IMATYPE', None) == n:
                    return e


        img = dict(
            wcs = pywcs.WCS(imatype_selector(f, 'INTENSITY').header),
            flux = imatype_selector(f, 'INTENSITY').data,
            var = imatype_selector(f, 'VARIANCE').data,
            ex = imatype_selector(f, 'EXPOSURE').data,
        )
        img['n'] = np.ones_like(img['ex'])

        if self.mock:
            img['var'] = 10*np.ones_like(img['var'])
            img['flux'] = -10*np.ones_like(img['var'])
            img['flux'][
                    (int(img['flux'].shape[0]/2)-10):(int(img['flux'].shape[0]/2)+10),
                    (int(img['flux'].shape[1]/2)-10):(int(img['flux'].shape[1]/2)+10),
                    ] = 100
        

        if self.mosaic is None: 
            print("first mosaic")
            self.mosaic = img
        else:
            print("adding mosaic")
            m_j, m_i = np.meshgrid(
                    np.arange(self.mosaic['flux'].shape[0]),
                    np.arange(self.mosaic['flux'].shape[1]),
                )

            m_ra, m_dec = self.mosaic['wcs'].wcs_pix2world(m_i, m_j, 0)
            i, j = img['wcs'].wcs_world2pix(m_ra, m_dec, 0)
            i = i.astype(int)
            j = j.astype(int)

            print("MOSAIC:", self.mosaic['wcs'])
            print("IMG:", img['wcs'])

            tm = i - m_i
            #tm = img['flux'][i, j]
            fits.ImageHDU(tm, header=self.mosaic['wcs'].to_header()).writeto("t.fits", overwrite=True)

            # or interpolate

            def pickornan(x, i, j):
                f_x = x.flatten()
                f_i = i.flatten()
                f_j = k.flatten()

                m = f_i > 0
                m &= f_i < x.shape[0] 
                m &= f_j > 0
                m &= f_j < x.shape[1] 

                y = np.zeros_like(f_x)

                y[m] = x[i, j].flatten()
                y[~m.flatten()] = float('NaN')
                return y.reshape(x.shape)

            self.mosaic.update(self.m(
                        self.mosaic,
                        lambda x:pickornan(x, m_i, m_j),
                        img,
                        lambda x:pickornan(x, i, j),
                    ))


    def writeto(self, fn):
        fits.ImageHDU(self.mosaic['flux'], header=self.mosaic['wcs'].to_header()).writeto(fn, overwrite=True)

@click.command()
@click.argument("in_fn", nargs=-1)
@click.argument("out_fn", nargs=1)
@click.option("--out-format", default="first", help="test")
def mosaic(in_fn, out_fn, out_format):
    mosaic_fn_list(in_fn, out_fn, out_format)


if __name__ == "__main__":
    mosaic()
