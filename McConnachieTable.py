#!/usr/bin/env python

from numpy import *
import matplotlib
matplotlib.use('Agg')
from pylab import *
import sys
import os
from astropy.io import ascii,fits
import astropy.units as u
import astropy.constants as const
from astropy.coordinates import SkyCoord, Galactic, Distance


class McConnachieTable:
    def __init__(self, base_dir="", R_sun=8.0):
        mcconnachie_table  = ascii.read(base_dir + "NearbyGalaxies.dat",
                                        data_start=31, format='fixed_width',
                                        header_start=29,delimiter='|')
        self.names = mcconnachie_table['GalaxyName']
        self.distance_modulus = array(mcconnachie_table['(m-M)'])
        self.Vmag = array(mcconnachie_table['Vmag'])
        rh_ = array(mcconnachie_table['rh(arcmins)'])
        self.rh_arcmin = array([float(x.partition(' ')[0]) for x in rh_])

        self.HI = array(mcconnachie_table['MHI'])

        self.M_V = self.Vmag - self.distance_modulus

        ML = 1.0
        self.mass = ML*10**((self.Vmag - self.distance_modulus - 4.8)/-2.5) * u.M_sun
        self.distance = 10**(self.distance_modulus/5 + 1)/1e3 * u.kpc
        self.rh_phys = (self.rh_arcmin*(60/206265.0)*self.distance).to(u.pc)

        RA = [15*(int(x[0]) + int(x[1])/60.0 + float(x[2])/3600.0)
              for x in [y.split() for y in mcconnachie_table['RA']]]
        self.RA = RA

        dec = [(abs(int(x[0])) + int(x[1])/60.0 + float(x[2])/3600.0)*sign(int(x[0]))
              for x in [y.split() for y in mcconnachie_table['Dec']]]
        self.dec = dec

        coords = SkyCoord(ra=RA, dec=dec, unit=(u.degree, u.degree),
                          frame='icrs',
                          distance=Distance(self.distance))
        self.coords = coords
        coords_gal = coords.transform_to(Galactic)
        self.MW_dist = sqrt((coords_gal.cartesian.x - R_sun*u.kpc)**2 +
                            coords_gal.cartesian.y**2 + coords_gal.cartesian.z**2)
        M31 = SkyCoord(ra=10.6847929, dec=41.2690650, unit=(u.degree, u.degree),
                       frame="icrs",
                       distance=Distance(745 * u.kpc))
        M31 = M31.transform_to(Galactic)
        self.M31_dist = sqrt((coords_gal.cartesian.x - M31.cartesian.x)**2 +
                             (coords_gal.cartesian.y - M31.cartesian.y)**2 +
                             (coords_gal.cartesian.z - M31.cartesian.z)**2)
        self.min_dist = minimum(self.MW_dist, self.M31_dist)

        self.M_HI = mcconnachie_table['MHI']*1e6
        sel_bad_HI, = where(mcconnachie_table['MHI'] == 99.9)
        self.M_HI[sel_bad_HI] = -999

        #f_extend = ascii.read(base_dir + "NearbyGalaxies_extended.dat")
        #self.pre_SDSS = f_extend['PreSDSS'] == 1
        self.pre_SDSS = zeros(len(self.RA))

if __name__ == "__main__":
    t = McConnachieTable()

    sel_MW, = where(t.MW_dist < 300*u.kpc)
    sel_MW_preSDSS, = where((t.MW_dist < 300*u.kpc) & (t.pre_SDSS))
    print "Pre SDSS: ", len(sel_MW_preSDSS)
    print "Total: ", len(sel_MW)

    figure(1).clear()
    loglog(t.mass.value, t.M_HI, 'bx', ms=10, mew=3)

    Mvir = array([2.5e9, 7.6e9, 7.7e9, 7.7e9])
    Mstars = array([0.00002, 0.0004, 0.0003, 0.0003])*Mvir
    Mgas = array([0.0049, 0.018, 0.014, 0.011])*Mvir
    loglog(Mstars, Mgas, 'r+', mew=4, ms=15)
    xlabel("Stellar mass")
    ylabel("HI gas mass")
    xlim(1e4,1e10)
    ylim(1e4,1e10)
    savefig("dwarf_HI.png")

    figure(1).clear()
    loglog(t.mass.value, t.M_HI/t.mass.value, 'bx', ms=10, mew=3)

    loglog(Mstars, Mgas/Mstars, 'r+', mew=3, ms=20)
    sel, = where(t.M_HI/t.mass.value > 5)
    for n in sel:
        text(t.mass.value[n], t.M_HI[n]/t.mass.value[n], t.names[n])
    xlabel("Stellar mass")
    ylabel("M gas / M stars")
    xlim(1e4,1e10)
    ylim(1e-3,1e3)
    savefig("dwarf_HI_2.png")


    #sel, = where(t.M31_dist < 400*u.kpc)
    #sel_sorted = sel[argsort(t.M31_dist[sel])]
    #for n in sel:
    #    print t.names[n]


