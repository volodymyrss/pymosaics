import click

from typing import List

import healpy

import numpy as np

def mosaic_fn_list(input_fn_list: List[str]):
    pass


class Mosaic:
    def writeto(self, fn):
        pass

class HealpixMosaic(Mosaic):
    def __init__(self, nsides):
        self.nsides = nsides
        self.mosaic = np.zeros(healpy.nsides2npix(nsides))

    def writeto(self, fn: str):
        healpy.write_map(fn, self.mosaic)

class FITsMosaic(Mosaic):
    def writeto(self, fn):
        pass

@click.command()
@click.argument("in_fn", nargs=-1)
@click.argument("out_fn", nargs=1)
@click.option("--format", default="first", help="test")
def mosaic(in_fn, out_fn, format):
    click.echo("input: %s output: %s"%( ", ".join(in_fn), out_fn))


if __name__ == "__main__":
    mosaic()
