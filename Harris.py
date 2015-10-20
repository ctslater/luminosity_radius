#!/usr/bin/env python

from numpy import *
import matplotlib
matplotlib.use('Agg')
from pylab import *
import sys
import os
from astropy.io import ascii,fits
import astropy.units as u
from astropy.coordinates import SkyCoord


class Harris:
    def __init__(self, base_dir=""):
        self.harris_tableI = ascii.read(base_dir + "mwgc.dat",
                                        data_start=48, data_end=205,
                                        delimiter=' ', guess=False,
                                        format="fixed_width_no_header",
                                        names=("ID", "Name", "RA", "Dec", "l",
                                               "b", "R_sun", "R_gc", "X", "Y",
                                               "Z"),
                                        col_starts=(0,10,23,37,50,59,67,74,80,86,92),
                                        col_ends=(9,22,36,49,58,66,73,79,85,91,97))

        self.harris_tableII = ascii.read(base_dir + "mwgc.dat", data_start=221,
                                         data_end=378, delimiter=' ',
                                         format="fixed_width_no_header",
                                         names=("ID", "[Fe/H]", "[Fe/H]wt",
                                                "EBV", "V_HB", "(m-M)V", "V_t",
                                                "M_V,t", "U-B", "B-V", "V-R",
                                                "V-I", "spt", "ellip"),
                                         col_starts=(0,12,19,22,29,35,41,47,54,61,67,
                                                     73,79,83),
                                         col_ends=(11,18,21,28,34,40,46,53,60,66,
                                                   72,78,82,90))

        self.harris_tableIII = ascii.read(base_dir + "mwgc.dat", data_start=395,
                                          data_end=552, delimiter=' ',
                                          format="fixed_width_no_header",
                                          names=("ID", "v_r", "v_r_pm", "v_LSR",
                                                 "sig_v", "sig_v_pm", "c",
                                                 "r_c_flag", "r_c", "r_h",
                                                 "mu_V", "rho_0", "lg(tc)",
                                                 "lg(th)"),
                                          col_starts=(0,10,19, 25,33,41,47,54,
                                                      58,64,70,78,83,90),
                                          col_ends=(9,18,24,32,40,46,53,57,63,
                                                    69,77,82,89,96))


        self.names = self.harris_tableII['ID']
        self.names2 = self.harris_tableI['Name']
        self.M_V = self.harris_tableII['M_V,t']

        self.r_h_arcmin = self.harris_tableIII['r_h']
        # Physical units
        #self.dist = 10**(0.2*self.harris_tableI['R_sun'])/100.0
        self.dist = self.harris_tableI['R_sun']

        self.r_h = self.r_h_arcmin*pi/(60*180.0)*self.dist*1000.0
        self.glat = self.harris_tableI['b']
        self.glon = self.harris_tableI['l']

        coords = SkyCoord(self.glon, self.glat, frame="galactic",
                          unit=(u.deg, u.deg))
        equ_coords = coords.transform_to('icrs')

        self.ra = equ_coords.ra
        self.dec = equ_coords.dec





