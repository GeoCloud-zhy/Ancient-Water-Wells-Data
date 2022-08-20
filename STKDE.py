"""时空核密度模型（Spatiotemporal Kernel Density Estimation）
Author:zeng hongyun
email:hy_zeng@ynu.edu.cn
Describtions:spatiotemporal kernel density estimation toolkit
"""

import matplotlib.pyplot as plt
import geopandas as gpd
from scipy import spatial
import contextily as ctx

import numpy as np


# gaussian kernel
# parameters describe:
# (x,y,t):the coordinate of the kernel
# (xi,yi,ti):the coordinate of the caculated grid center
# hs：spatial bandwidth
# ht：temporal bandwidth
def gaussian(x, y, t, xi, yi, ti, hs, ht):
    u = (x - xi) / hs
    v = (y - yi) / hs
    w = (t - ti) / ht

    Ks = (1. / (2 * np.pi)) * np.exp(-(u ** 2) / 2.) * np.exp(-(v ** 2) / 2.)
    Kt = (1. / (np.sqrt(2 * np.pi))) * np.exp(-(w ** 2) / 2.)

    st = [Ks, Kt]

    return st


# An iterative plug-in bandwidth selector
# x: array-like, Array for which to get the bandwidth
def plugin_bandwidth(x):
    I_7 = 1. / 7.
    I_9 = 1. / 9.
    rt2pi = np.sqrt(2. * np.pi)
    rtpi2 = 2. * np.sqrt(np.pi)
    iter = 5

    nx = len(x)
    n2 = nx ** 2
    xiqr = x[(3 * nx) / 4 - 1] - x[nx / 4]  # = IQR(x[])
    # estimate inflation constant c
    h2 = (0.920 * xiqr) / np.power(nx, I_7)
    h3 = (0.912 * xiqr) / np.power(nx, I_9)

    s2 = s3 = 0.
    for i in range(nx - 1):
        for j in range(i + 1, nx):
            t = x[i] - x[j]
            a = t / h2
            d2 = a * a
            a = t / h3
            d3 = a * a

            if ((d2 > 50) and (d3 > 60)):
                break

            s2 += np.exp(-d2 / 2.) * (3. + d2 * (-6. + d2))
            s3 += np.exp(-d3 / 2.) * (-15. + d3 * (45. + d3 * (-15. + d3)))

    rhat2 = 2. * s2 / (rt2pi * n2 * np.power(h2, 5)) + 3. / (rt2pi * nx * np.power(h2, 5))
    rhat3 = -2. * s3 / (rt2pi * n2 * np.power(h3, 7)) + 15. / (rt2pi * nx * np.power(h3, 7))
    co1 = 1.357 * np.power(rhat2 / rhat3, I_7)
    co2 = 1. / rtpi2
    a = 1.132795764 / (np.power(rhat3, I_7) * np.sqrt(nx))

    # loop over iterations
    for it in range(1, iter + 1):
        s2 = 0.
        for i in range(nx - 1):
            for j in range(i + 1, nx):
                t = (x[i] - x[j]) / a
                d2 = t * t
                if d2 > 50: break
                s2 += np.exp(-d2 / 2.) * (3. + d2 * (-6. + d2))

        a = rt2pi * nx * np.power(a, 5)
        rhat2 = 2 * s2 / (nx * a) + 3. / a

        # estimate bandwidth by asymptotic formula
        a = co1 * np.power(co2 / (rhat2 * nx), I_7)
    hs = np.power(co2 / (rhat2 * nx), 0.2)

    return hs


# construct the grids
def construct_grids(xmin, ymin, xmax, ymax, grid_size):
    xgrid = np.arange(xmin, xmax, grid_size)
    ygrid = np.arange(ymin, ymax, grid_size)

    return (xgrid, ygrid)


class SpatiotemporalKernelDensity(object):
    """Spatiotemporal Kernel Density Estimation Toolkit

    """
    def __init__(self):
        self.__featureclass = ""
        self.__temporal_field = ""
        self.__output_raster_filename = ""
        self.__spatial_bandwidth = 0.
        self.__temporal_bandwidth = 0.

    @property.setter
    def FeatureClass(self, featureclass):
        if featureclass != "":
            self.__gdf = gpd.read_file(featureclass)
            if self.__gdf.geom_type == 'Point':
                self.__featureclass = featureclass
            else:
                raise IOError

    @property.setter
    def Temporal_Field(self, temporal_field):
        self.__temporal_field = temporal_field

    @property.setter
    def Output_Raster_FileName(self, output_raster_filename):
        self.__output_raster_filename = output_raster_filename

    def doSTKDE(self):
        if self.__gdf == None:
            return
        



