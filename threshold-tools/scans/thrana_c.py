#!/usr/bin/env python3

import os
import argparse
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid1 import make_axes_locatable


def gaus(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


class ThresholdScan:
    def __init__(self, thrmap):
        self.thrname = thrmap
        self.fthr = open(self.thrname, 'r')
        self.thrdata = None

    def main(self, vi, vf, xmin=-1, xmax=-1, rows=512, cols=1024, verbose=True, fname=None):

        if xmin == -1:
            xmin = vi
        if xmax == -1:
            xmax = vf

        self.thrdata = np.fromfile(self.fthr, dtype=np.float32)
        self.fthr.close()

        thrmap = np.zeros((3, 512,1024), dtype=np.float32)

        # fill maps
        for j in range(512):
            for k in range(1024*3):
                i = k // 1024
                thrmap[i][j][k - i*1024] = self.thrdata[j*3*1024 + k]

        xmax = xmax if xmax else vf
        print('filling histogram')

        for chip in range(0,3):
            n_hits = np.sum(thrmap[chip] > 0)
            if n_hits == 0:  # no data for this chip
                continue
            thrs = []
            npix = 0
            nan_count = 0
            for r in range(rows):
                for c in range(cols):
                    m = thrmap[chip,r,c]
                    if (m < 0) or (m > float(nsteps)):
                        continue
                    if (c % 32 == 0) and (r % 16 == 0):
                        npix += 1

                    if m is np.nan:
                        nan_count += 1
                        continue
                    thrs.append(float(m+vmin))

            thrs[:] = [x for x in thrs if (x >= xmin) and (x <= xmax)]
            nbins = 40
            npix = rows*cols
            plt.figure(f"threshold {chip}")
            plt.xlabel('Threshold (e$^-$)')
            plt.ylabel(f'# pixels / ({(xmax-xmin)/nbins:.1f} e$^-$)')
            plt.title(f'Chip {chip} Threshold distribution ({npix} pixels)')
            hist,bin_edges,_ = plt.hist(thrs, bins=nbins,
                                        label=f"Mean: {np.mean(thrs):5.1f} e$^-$\nRMS:  {np.std(thrs):5.1f} e$^-$")
            bin_mid = (bin_edges[:-1] + bin_edges[1:])/2
            try:
                popt,_  = curve_fit(gaus, bin_mid, hist, [10, np.mean(thrs), np.std(thrs)])
                plt.plot(np.arange(xmin, xmax, (xmax-xmin)/nbins), gaus(np.arange(xmin, xmax, (xmax-xmin)/nbins), *popt),
                         label=f'$\mu$:    {popt[1]:5.1f} e$^-$\n$\sigma$:    {popt[2]:5.1f} e$^-$')
            except Exception as e:
                if verbose:
                    print("Fitting error", e)

            plt.legend(loc="upper right", prop={"family":"monospace"})
            plt.savefig(fname+f"_{chip}_threshold.png")

            cmap = plt.colormaps["viridis"].copy()
            cmap.set_under(color='white')

            plt.figure(f"Chip {chip} Threshold map")
            thrmap[thrmap==0] = np.nan
            plt.subplots_adjust(left=0.085, right=0.85)
            ax = plt.gca()
            plt.imshow(thrmap[chip], cmap=cmap, vmin=2.0, vmax=25.0)
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            plt.colorbar(format='%.0e', cax=cax).set_label('Threshold (e$^-$)')
            ax.set_xlabel('Column')
            ax.set_ylabel('Row')
            ax.set_title(f'Chip {chip} Threshold map')
            plt.savefig(fname+f"_{chip}_thrmap.png")

        if not verbose:
            plt.close('all')


if __name__=="__main__":
    parser = argparse.ArgumentParser("Threshold analysis.",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-thr", help="Threshold file", required=True)
    parser.add_argument("-o", help="Output file name", required=True)
    parser.add_argument("-vi", help="Minimum threshold", type=int, default=120, required=False)
    parser.add_argument("-vf", help="Maximum threshold", type=int, default=169, required=False)
    parser.add_argument("-nsteps", help="Threshold step", type=int, default=50, required=False)
    parser.add_argument('-xmin', default=0, type=int, help="X axis low limit")
    parser.add_argument('-xmax', default=0, type=int, help="X axis high limit (0 = use vmax)")
    parser.add_argument('-row', default=512, type=int, help="Rows to scan")
    parser.add_argument('-col', default=1024, type=int, help="Columns to scan")
    parser.add_argument('-q', '--quiet', action='store_true', help="Do not display plots.")
    args = parser.parse_args()

    thrfile = args.thr
    output = args.o
    vmin = args.vi
    vmax = args.vf
    nsteps = args.nsteps
    xmin = args.xmin
    xmax = args.xmax
    row = args.row
    col = args.col
    verbose = not args.quiet

    print("Threshold file:", thrfile)
    print("Output file:", output)
    print("Minimum threshold:", vmin)
    print("Maximum threshold:", vmax)
    print("Threshold step:", nsteps)
    print("X axis low limit:", xmin)
    print("X axis high limit:", xmax)
    print("Rows to scan:", row)
    print("Columns to scan:", col)

    plotter = ThresholdScan(thrfile)

    plotter.main(vi=vmin, vf=vmax, xmin=xmin, xmax=xmax, rows=row, cols=col, verbose=verbose, fname=output)

    if not args.quiet:
        plt.show()
