import click

from typing import List, Union

from astropy import wcs as pywcs # type: ignore
from astropy.io import fits # type: ignore
from astropy.io.fits import HDUList # type: ignore

import healpy # type: ignore

import numpy as np # type: ignore

debug = False

def mosaic_list(in_fs: List[Union[str, fits.HDUList]], 
                out_fn: Union[str, None]=None, 
                pixels="healpix", 
                healpix_nsides=512, 
                mock=False):
    print("input: %s output: %s"%( ", ".join(map(str, in_fs)), out_fn))

    if pixels == "first":
        m = FITsMosaic(mock) # type: Mosaic
        for in_f in in_fs:
            m.add(in_f)

    elif pixels == "healpix":
        m = HealpixMosaic(nsides=healpix_nsides, mock=mock)
        for in_f in in_fs:
            m.add(in_f)

    else:
        raise RuntimeError

    if out_fn is not None:
        m.writeto(out_fn)

    return m

def mosaic_header_list(*a, **aa):
    raise RuntimeError("this is a compatibility placeholder, please use mosaic_list")

def mosaic_fn_list(*a, **aa):
    raise RuntimeError("this is a compatibility placeholder, please use mosaic_list")

def pickornan(x, i, j) -> np.ndarray:
    ij_usable = i>0
    ij_usable &= j>0
    ij_usable &= i<x.shape[0]
    ij_usable &= j<x.shape[1]

    # pylint: disable=unsupported-assignment-operation
    y = np.zeros_like(i, dtype=np.float32) 
    y[ij_usable] = x[i[ij_usable], j[ij_usable]]
    y[~ij_usable] = np.nan
    return y 


class Mosaic:
    mock = False

    def writeto(self, fn: str):
        raise RuntimeError

    def add(self, f: Union[str, HDUList] ):
        pass


    def parse_in(self, in_f):
        if isinstance(in_f, str):
            print("parsing", in_f)
            fn = in_f
            f = fits.open(in_f)

        elif isinstance(in_f, HDUList):
            f = in_f
            try:
                fn = in_f.filename()
                print("parsing hdulist", fn)
            except Exception as e:
                fn = "unknown"
                print("parsing hdulist with no detectable filename?")

        else:
            raise RuntimeError("unknown input " +repr(f))


        img = None
        failures=[]
        for n, img_parser in img_parsers.items():
            try:
                img = img_parser(f)
                print("img parser succeeded", n)
                break
            except Exception as e:
                print("img parser failed", n, e)
                failures.append("img parser failed %s %s"%(n, e))

        if img is None:
            raise RuntimeError("all image parsers failed for %s: %s"%(f, ", ".join(failures)))

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

        return img

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
        _m_a = ~np.isnan(_s_a(a['flux']))
        _m_b = ~np.isnan(_s_b(b['flux']))
        _m = _m_a & _m_b
        
        c=dict()

        print("m: _m shape", _m.shape)
        print("m: _m_a shape", _m_a.shape)
        print("m: _m_b shape", _m_b.shape)

        for k in 'flux', 'var', 'ex', 'n':
            c[k] = _s_a(a[k])
            c[k][~_m_a & _m_b] = _s_b(b[k])[~_m_a & _m_b]


        print("max exposure was", np.nanmax(a['ex']))
        print("max exposure to add", np.nanmax(b['ex']))

        s_a = lambda x:_s_a(x)[_m]
        s_b = lambda x:_s_b(x)[_m]


        c['flux'][_m] = s_a(a['flux']) / s_a(a['var']) + s_b(b['flux']) / s_b(b['var'])
        c['var'][_m] = ( 1 / s_a(a['var']) + 1 / s_b(b['var']) ) ** -1
        c['flux'][_m] = c['flux'][_m] * c['var'][_m]
        c['ex'][_m] = s_a(a['ex']) + s_b(b['ex'])
        c['n'][_m] = s_a(a['n']) + s_b(b['n'])

        print("max exposure become", np.nanmax(c['ex']))

        return c


def img_parser_own(f):  #eminmax
    return None

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

    img['flux'][img['var']<=0] = np.nan 

    return img

img_parsers = {
        'osa': img_parser_osa,
        'own': img_parser_own,
}


class FITsMosaic(Mosaic):
    def __init__(self, mock=False):
        self.mosaic = None
        self.mock = mock

    def add(self, in_f):
        super(FITsMosaic, self).add(in_f)

        img = self.parse_in(in_f)


        

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


            self.mosaic.update(self.m(
                        self.mosaic,
                        lambda x:pickornan(x, m_i, m_j),
                        img,
                        lambda x:pickornan(x, i, j),
                    ))

    def to_hdu_list(self):
        self.mosaic['sig'] = self.mosaic['flux'] / self.mosaic['var'] ** 0.5

        el = []

        for k in self.mosaic:
            if k == 'wcs': continue

            h = self.mosaic['wcs'].to_header()
            h['EXTNAME'] = dict(
                flux='INTENSITY',
                var='VARIANCE',
                sig='SIGNIFICANCE',
                ex='EXPOSURE',
                n='NIMAGE',
            )[k]
            h['IMATYPE'] = h['EXTNAME']
            e = fits.ImageHDU(self.mosaic[k], header=h)

            el.append(e)

        return fits.HDUList([fits.PrimaryHDU()] + el)


    def writeto(self, fn):
        self.to_hdu_list().writeto(fn, overwrite=True)


class HealpixMosaic(Mosaic):
    def __init__(self, nsides=512, mock=False):
        self.nsides = nsides
        self.mosaic = None
        self.mock = mock

    def radec_grid(self):
        px = np.arange(healpy.nside2npix(self.nsides))
        return healpy.pix2ang(self.nsides, px, lonlat=True)

    def add(self, in_f):
        img = self.parse_in(in_f)

        ra, dec = self.radec_grid()

        i, j = img['wcs'].wcs_world2pix(ra, dec, 0)
        i = i.astype(np.int32)
        j = j.astype(np.int32)
        

        def tohp(x):
            return pickornan(x, j, i).flatten() # pylint: disable=no-member

        if self.mosaic is None:
            self.mosaic = dict()

            for k, v in img.items():
                if k == 'wcs': continue
                self.mosaic[k] = tohp(img[k])
        else:
            self.mosaic.update(self.m(
                        self.mosaic,
                        lambda x:x,
                        img,
                        tohp,
                    ))


    def writeto(self, fn: str):
        # only import this as needed
        import matplotlib # type: ignore
        matplotlib.use('agg')
        from matplotlib import pylab as plt

        f = plt.figure()
        healpy.mollview(self.mosaic['flux'])
        plt.savefig("out.png")

  #      healpy.write_map(fn, self.mosaic['flux'], overwrite=True)

        self.writeto_reproj(fn)

    def to_hdu_list(self):
        #self.mosaic['flux'] overwrite=True)

        ra, dec = self.radec_grid()

        print("max exposure", np.nanmax(self.mosaic['ex']))

        m_best = self.mosaic['ex']>(np.nanmax(self.mosaic['ex'])/2.)
        c_ra = np.nanmean(ra[m_best])
        c_dec = np.nanmean(dec[m_best])
        c_rad_deg = 60

        print("center ra, dec", c_ra, c_dec)

        w = pywcs.WCS(naxis=2)
        ni, nj = 800, 800

        pxsize = c_rad_deg/ni

        w.wcs.crpix = int(ni/2), int(nj/2)
        w.wcs.crval = c_ra, c_dec
        w.wcs.cunit = 'deg', 'deg'
        w.wcs.cdelt = np.array([-pxsize, pxsize])
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        #w.wcs.set_pv([(2, 1, 45.0)])

        m_j, m_i = np.meshgrid(
                np.arange(nj),
                np.arange(ni),
            )

        m_ra, m_dec = w.wcs_pix2world(m_j, m_i, 0)

        px = healpy.ang2pix(self.nsides, m_ra, m_dec, lonlat=True)


        self.mosaic['sig'] = self.mosaic['flux'] / self.mosaic['var']**0.5
        el = []
        for k in self.mosaic:
            if k == 'wcs': continue

            h = w.to_header()
            h['EXTNAME']=dict(
                        flux='INTENSITY',
                        var='VARIANCE',
                        sig='SIGNIFICANCE',
                        ex='EXPOSURE',
                        n='NIMAGE',
                    )[k]
            h['IMATYPE'] = h['EXTNAME']

            mp = np.zeros((ni, nj))
            mp[:,:] = np.nan 
            mp[m_i, m_j] = self.mosaic[k][px]

            e = fits.ImageHDU(mp, header=h)

            el.append(e)
        return fits.HDUList([fits.PrimaryHDU()] + el)

    def writeto_reproj(self, fn):
        self.to_hdu_list().writeto(fn, overwrite=True)



@click.command()
@click.argument("in_fn", nargs=-1)
@click.argument("out_fn", nargs=1)
#@click.option("--out-format", default="first", help="from frist image, new map, healpix")
@click.option("--pixels", default="first", help="pixelization for building the map")
@click.option("--healpix-nsides", default=512, help="pixelization for building the map")
@click.option("--mock/--no-mock", default=False, help="mock or not")
def mosaic(in_fn, out_fn, pixels, healpix_nsides, mock):
    mosaic_fn_list(in_fn, out_fn, pixels=pixels, healpix_nsides=healpix_nsides, mock=mock)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    mosaic()
