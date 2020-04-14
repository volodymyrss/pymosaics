import click

from typing import List

from astropy import wcs as pywcs
from astropy.io import fits

import healpy

import numpy as np

debug = False

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
    def m_s(a, _s_a, b, _s_b):
        c=dict()

        _m = ~np.isnan(_s_b(a['flux']))

        s_a = lambda x:_s_a(x)[_m]
        s_b = lambda x:_s_b(x)[_m]

        c['flux'] = a['flux']
        c['var'] = a['var']
        c['ex'] = a['ex']
        c['n'] = a['n']

        c['flux'][_m] = s_a(a['flux']) + s_b(b['flux']) 
        c['var'][_m] = ( 1 / s_a(a['var']) + 1 / s_b(b['var']) ) ** -1
        c['ex'][_m] = s_a(a['ex']) + s_b(b['ex'])
        c['n'][_m] = s_a(a['n']) + s_b(b['n'])


        return c

    @staticmethod
    def m(a, _s_a, b, _s_b):
        _m = ~np.isnan(_s_b(b['flux']))

        s_a = lambda x:_s_a(x)[_m]
        s_b = lambda x:_s_b(x)[_m]

        c=dict()
        c['flux'] = a['flux']
        c['var'] = a['var']
        c['ex'] = a['ex']
        c['n'] = a['n']

        c['flux'][_m] = s_a(a['flux']) / s_a(a['var']) + s_b(b['flux']) / s_b(b['var'])
        c['var'][_m] = ( 1 / s_a(a['var']) + 1 / s_b(b['var']) ) ** -1
        c['flux'][_m] = c['flux'][_m] * c['var'][_m]
        c['ex'][_m] = s_a(a['ex']) + s_b(b['ex'])
        c['n'][_m] = s_a(a['n']) + s_b(b['n'])

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


def img_parser_osa(f):  #eminmax
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

    img['flux'][img['var']<=0] = float('NaN')

    return img

img_parsers = {
        'osa': img_parser_osa,
}


class FITsMosaic(Mosaic):
    def __init__(self, mock=False):
        self.mosaic = None
        self.mock = mock

    def add(self, fn):
        print("adding", fn)

        f = fits.open(fn)

        for n, img_parser in img_parsers.items():
            try:
                img = img_parser(f)
                print("img parser succeeded", n)
                break
            except Exception as e:
                print("img parser failed", n, e)

        if self.mock:
            img['var'] = 10*np.ones_like(img['var'])
            img['flux'] = -10*np.ones_like(img['var'])

            try:
                sj, si = map(int, img['wcs'].wcs_world2pix(83, 22, 0))
                sjd, sid = map(int, img['wcs'].wcs_world2pix(83, 24, 0))

                print("mock crab at", si, sj, "and", sid, sjd)

                img['flux'][
                        si-4:si+4,
                        sj-4:sj+4,
                        ] = 40
                img['flux'][
                        sid-4:sid+4,
                        sjd-4:sjd+4,
                        ] = 60
            except Exception as e:
                print("unable to add mock crab", e)


            img['flux'][
                    (int(img['flux'].shape[0]/2)-10):(int(img['flux'].shape[0]/2)+10),
                    (int(img['flux'].shape[1]/2)-10):(int(img['flux'].shape[1]/2)+10),
                    ] = 10
            img['flux'][
                    (int(img['flux'].shape[0]/2)-3):(int(img['flux'].shape[0]/2)+3),
                    (int(img['flux'].shape[1]/2)):(int(img['flux'].shape[1]/2)+30),
                    ] = 20
            fits.ImageHDU(img['flux'], header=img['wcs'].to_header()).writeto("mock-%s.fits"%fn.replace("/","_"), overwrite=True)
        

        if self.mosaic is None: 
            print("first mosaic")
            self.mosaic = img
        else:
            print("adding mosaic")
            m_j, m_i = np.meshgrid(
                    np.arange(self.mosaic['flux'].shape[1]),
                    np.arange(self.mosaic['flux'].shape[0]),
                )

            m_ra, m_dec = self.mosaic['wcs'].wcs_pix2world(m_j, m_i, 0)
            j, i = img['wcs'].wcs_world2pix(m_ra, m_dec, 0)
            i = i.astype(int)
            j = j.astype(int)

            print("MOSAIC:", self.mosaic['wcs'])
            print("IMG:", img['wcs'])

            def pickornan(x, i, j):
                ij_usable = i>0
                ij_usable &= j>0
                ij_usable &= i<x.shape[0]
                ij_usable &= j<x.shape[1]

                y = np.zeros_like(i, dtype=np.float32)
                y[ij_usable] = x[i[ij_usable], j[ij_usable]]
                y[~ij_usable] = float('nan')
                return y

            tm = pickornan(self.mosaic['flux'], m_i, m_j)
            fits.ImageHDU(tm - self.mosaic['flux'], header=self.mosaic['wcs'].to_header()).writeto("tmd.fits", overwrite=True)

            
            tm = pickornan(self.mosaic['flux'], m_i, m_j)
            fits.ImageHDU(tm, header=self.mosaic['wcs'].to_header()).writeto("tm.fits", overwrite=True)

            tm = pickornan(img['flux'], i, j)
            fits.ImageHDU(tm, header=self.mosaic['wcs'].to_header()).writeto("t.fits", overwrite=True)

            if debug:
                tm = m_i 
                fits.ImageHDU(tm, header=self.mosaic['wcs'].to_header()).writeto("mi.fits", overwrite=True)
                
                tm = m_j 
                fits.ImageHDU(tm, header=self.mosaic['wcs'].to_header()).writeto("mj.fits", overwrite=True)

                tm = self.mosaic['flux']
                fits.ImageHDU(tm, header=self.mosaic['wcs'].to_header()).writeto("tmf.fits", overwrite=True)
                
                tm = self.mosaic['flux'][m_i, m_j]
                fits.ImageHDU(tm, header=self.mosaic['wcs'].to_header()).writeto("tmr.fits", overwrite=True)

                tm = pickornan(img['flux'], i, j)
                fits.ImageHDU(tm, header=self.mosaic['wcs'].to_header()).writeto("t.fits", overwrite=True)

                tm=np.zeros_like(i)
                tm[ij_usable] = img['flux'][i[ij_usable], j[ij_usable]]
                tm[~ij_usable] = -20
                fits.ImageHDU(np.transpose(tm), header=self.mosaic['wcs'].to_header()).writeto("t.fits", overwrite=True)

                tm = i 
                tm[~ij_usable] = -100
                #tm = img['flux'][i, j]
                fits.ImageHDU(np.transpose(tm), header=self.mosaic['wcs'].to_header()).writeto("i.fits", overwrite=True)
                
                tm = j 
                tm[~ij_usable] = -100
                #tm = img['flux'][i, j]
                fits.ImageHDU(np.transpose(tm), header=self.mosaic['wcs'].to_header()).writeto("j.fits", overwrite=True)

            self.mosaic.update(self.m(
                        self.mosaic,
                        lambda x:pickornan(x, m_i, m_j),
                        img,
                        lambda x:pickornan(x, i, j),
                    ))


    def writeto(self, fn):
        self.mosaic['sig'] = self.mosaic['flux'] / self.mosaic['var']**0.5

        fits.HDUList([fits.PrimaryHDU()] + [
            fits.ImageHDU(self.mosaic[k], header=self.mosaic['wcs'].to_header())
                for k in self.mosaic if k != 'wcs'
        ]).writeto(fn, overwrite=True)

@click.command()
@click.argument("in_fn", nargs=-1)
@click.argument("out_fn", nargs=1)
@click.option("--out-format", default="first", help="test")
def mosaic(in_fn, out_fn, out_format):
    mosaic_fn_list(in_fn, out_fn, out_format)


if __name__ == "__main__":
    mosaic()
