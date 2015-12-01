#!/usr/bin/env python

from numpy import *
import matplotlib
matplotlib.use('Agg')
from pylab import *
import sys
import os
from astropy.io import ascii,fits
import astropy.units as u
from matplotlib.ticker import FormatStrFormatter

from McConnachieTable import McConnachieTable
from Harris import Harris

have_sb = False

def save_latest(n, sb=False):
    global have_sb
    if sb and not have_sb:
        draw_sb()
        have_sb = True
    legend(frameon=False, numpoints=1, loc="upper left", ncol=2,
           fontsize=10)
    xlim(1,5000)
    ylim(2,-19)
    xlabel("Half light radius (pc)")
    ylabel("$M_V$")
    gca().xaxis.set_major_formatter(FormatStrFormatter("%d"))
    savefig("lum_radius%02d.pdf" % n)
    print "Wroteout lum_radius%02d.pdf" % n

def draw_sb():
    x = logspace(0,4,30)
    mu = 24 # mag/arcsec^2
    plot(x, mu - 5*log10(x/10.0*206265) + 2.5*log(0.5), 'k--')
    text(10**((1.5 - mu - 2.5*log10(0.5))/-5.0)*10/206265.0, -1.0,
         "$\mu=%.0f$" % mu, fontsize=10, rotation=38)
    mu = 28 # mag/arcsec^2
    plot(x, mu - 5*log10(x/10.0*206265) + 2.5*log(0.5), 'k--')
    text(10**((1.5 - mu - 2.5*log10(0.5))/-5.0)*10/206265.0, -1.0,
         "$\mu=%.0f$" % mu, fontsize=10, rotation=38)
    mu = 30 # mag/arcsec^2
    plot(x, mu - 5*log10(x/10.0*206265) + 2.5*log(0.5), 'k--')
    text(10**((1.5 - mu - 2.5*log10(0.5))/-5.0)*10/206265.0, -1.0,
         "$\mu=%.0f$" % mu, fontsize=10, rotation=38)
    mu = 32 # mag/arcsec^2
    plot(x, mu - 5*log10(x/10.0*206265) + 2.5*log(0.5), 'k--')
    text(10**((1.5 - mu - 2.5*log10(0.5))/-5.0)*10/206265.0, -1.0,
         "$\mu=%.0f$" % mu, fontsize=10, rotation=38)

rc("lines", markersize=8, markeredgewidth=0)
if __name__ == '__main__':

    t = McConnachieTable()
    mackey = ascii.read("mackey.dat", fill_values=('-','-999'))
    mackey_M_V = mackey['col7']
    mackey_rh = mackey['col9']

    oleg_regions = ascii.read("oleg_regions.dat")

    portegieszwart = ascii.read("portegieszwart10.dat")

    # M_V, rh
    FFs = [(-4.7, 24), (-4.7, 19), (-4.9, 9),
           (-5.0, 6), (-5.05, 7), (-5.1, 20),
           (-5.25, 14), (-5.4, 15), (-5.6, 16),
           (-5.7, 14), (-5.75, 9), (-5.8, 17),
           (-6.7, 19), (-7.1, 10), (-9.5, 16)]


    # Evestigneeva 2008, Table 1
    ucd_M_V = array([-13.33, -12.58, -11.08, -10.86, -10.49, -11.08, -11.50,
                     -10.76, -11.28, -11.24, -10.87, -10.81, -11.09, -12.24,
                     -12.23, -12.59, -12.25, -12.32, -12.07, -13.42, -11.95,
                     -12.19, -12.27, -12.45, -11.99])

    ucd_r_h = array([86.5, 10.3, 6.4, 11.8, 7.0, 11.4, 6.9, 5.8, 5.1, 10.9, 8.7,
                     4.0, 9.3, 11.2, 13.1, 20.0, 23.2, 17.8, 17.4, 93.2, 23.5,
                     22.4, 23.1, 29.5, 25.0])

    harris = Harris()

    figure(1).clear()

    pre_sdss = (t.pre_SDSS == True)

    sel_M31, = where((t.M31_dist < t.MW_dist) & (t.min_dist < 300*u.kpc) &
                     pre_sdss)
    semilogx(t.rh_phys.value[sel_M31], t.M_V[sel_M31], 'bs', label="M31")

    sel_MW, = where((t.M31_dist >= t.MW_dist) & (t.min_dist < 300*u.kpc) &
                    pre_sdss)
    semilogx(t.rh_phys.value[sel_MW], t.M_V[sel_MW], 'ro', label="MW")

    sel_LG, = where((t.min_dist >= 300*u.kpc) & (t.min_dist < 1500*u.kpc) &
                    pre_sdss)
    semilogx(t.rh_phys.value[sel_LG], t.M_V[sel_LG], 'g^', label="LG")

    sel_NB, = where((t.min_dist >= 1500*u.kpc))
    semilogx(t.rh_phys.value[sel_NB], t.M_V[sel_NB], 'm*', ms=12, label="Nearby")

    #plot(mackey_rh, mackey_M_V, 'k^', label="M31 GC", markersize=4)
    plot(harris.r_h, harris.M_V, 'bo', label="MW GC",
         markersize=4)

    save_latest(1)

    figure(1).clear()

    pre_sdss = (t.pre_SDSS == True)

    sel_M31, = where((t.M31_dist < t.MW_dist) & (t.min_dist < 300*u.kpc) &
                     pre_sdss)
    semilogx(t.rh_phys.value[sel_M31], t.M_V[sel_M31], 'bs', label="M31")

    sel_MW, = where((t.M31_dist >= t.MW_dist) & (t.min_dist < 300*u.kpc))
    semilogx(t.rh_phys.value[sel_MW], t.M_V[sel_MW], 'ro', label="MW")

    sel_LG, = where((t.min_dist >= 300*u.kpc) & (t.min_dist < 1500*u.kpc))
    semilogx(t.rh_phys.value[sel_LG], t.M_V[sel_LG], 'g^', label="LG")

    sel_NB, = where((t.min_dist >= 1500*u.kpc))
    semilogx(t.rh_phys.value[sel_NB], t.M_V[sel_NB], 'm*', ms=12, label="Nearby")

    #plot(mackey_rh, mackey_M_V, 'k^', label="M31 GC", markersize=4)
    plot(harris.r_h, harris.M_V, 'bo', label="MW GC",
         markersize=4)

    save_latest(2)

    figure(1).clear()
    sel_M31, = where((t.M31_dist < t.MW_dist) & (t.min_dist < 300*u.kpc))
    semilogx(t.rh_phys.value[sel_M31], t.M_V[sel_M31], 'bs', label="M31")

    sel_MW, = where((t.M31_dist >= t.MW_dist) & (t.min_dist < 300*u.kpc))
    semilogx(t.rh_phys.value[sel_MW], t.M_V[sel_MW], 'ro', label="MW")

    sel_LG, = where((t.min_dist >= 300*u.kpc) & (t.min_dist < 1500*u.kpc))
    semilogx(t.rh_phys.value[sel_LG], t.M_V[sel_LG], 'g^', label="LG")

    sel_NB, = where((t.min_dist >= 1500*u.kpc))
    semilogx(t.rh_phys.value[sel_NB], t.M_V[sel_NB], 'm*', ms=12, label="Nearby")

    #plot(mackey_rh, mackey_M_V, 'k^', label="M31 GC", markersize=4)
    plot(harris.r_h, harris.M_V, 'bo', label="MW GC",
         markersize=4)

    save_latest(3)

    plot(mackey_rh, mackey_M_V, 'k^', label="M31 GC", markersize=4)
    save_latest(4)

    save_latest(5, sb=True)

    plot(ucd_r_h, ucd_M_V, 'mx',mew=1, ms=6, label="UCDs")

    save_latest(6, sb=True)


    plot(226, -9.85, 'y*', markersize=20,mew=1)

    plot(20, -5.0, 'y*', markersize=20,mew=1)

    plot(340, -10.3, 'y*', markersize=20,mew=1)

    save_latest(7, sb=True)

    fill_between([1.5e3,4.5e3], [-12.5,-12.5], [-16.0, -16.0], color='k', alpha=0.3)
    #text(2e3, -11.0, "Coma Galaxies\nfrom Dragonfly", fontsize=10)

    labeled_galaxies = ["Sagittarius dSph", "WLM", "NGC 3109", "IC 1613"]
    for n in range(len(t.names)):
        if t.names[n] in labeled_galaxies:
            annotate(t.names[n], (t.rh_phys[n].value, t.M_V[n]),
                         xytext=(3,3), textcoords='offset points', fontsize=8)


    save_latest(8, sb=True)

    plot(3e3, -12.85, 'y*', ms=20, mew=2)

    mag_vals = [-11.1, -11.2]
    size_vals =[1.75e3, 1.5e3]
    for mag,size in zip(mag_vals, size_vals):
        plot(size, mag, 'k_', mew=3)
        plot([size, size], [mag, mag + 1.0], 'k-')
        plot(size, mag + 1.0, 'kv')

    save_latest(9, sb=True)

    #annotate("Kim 1", (6.9, 0.3),
    #         xytext=(3,3), textcoords='offset points', fontsize=10)
    plot(6.9, 0.3, 'rx', ms=13, mew=3)

    #annotate("Kim 2", (12.8, -1.5),
    #         xytext=(3,3), textcoords='offset points', fontsize=10)
    plot(12.8, -1.5, 'rx', ms=13, mew=3)

    save_latest(10, sb=True)

    cluster_ML = 2.0
    solarMV = 4.7
    semilogx(portegieszwart['radius'],
         -2.5*log10(portegieszwart['mass']/cluster_ML) + solarMV, 'c.',
         label='OCs')

    save_latest(11, sb=True)

    GC_ML = 3.0
    plot(oleg_regions ['r_h'], -2.5*log10(oleg_regions['mass']/GC_ML) + 4.7, 'k-')

    save_latest(12, sb=True)

    annotate("Kim 1", (6.9, 0.3),
             xytext=(3,3), textcoords='offset points', fontsize=10)
    annotate("Kim 2", (12.8, -1.5),
             xytext=(3,3), textcoords='offset points', fontsize=10)

    annotate("G1", (226,-9.85), xytext=(5,5), textcoords='offset points',
             fontsize=14)
    annotate("G2", (20,-5.0), xytext=(5,5), textcoords='offset points',
             fontsize=14)
    annotate("253 Dw1", (340,-10.3), xytext=(5,5), textcoords='offset points',
             fontsize=14)

    for n in range(len(t.names)):
        annotate(t.names[n], (t.rh_phys[n].value, t.M_V[n]),
                     xytext=(3,3), textcoords='offset points', fontsize=8)

    save_latest(0, sb=True)

    sys.exit(0)

    # Jay Strader's "The Densest Galaxy"
    #plot(24, -14.2, 'y*', markersize=20,mew=1)
    #annotate("M60-UCD1", (24,-14.2), xytext=(5,5), textcoords='offset points',
    #         fontsize=14)


    plot([x[1] for x in FFs], [x[0] for x in FFs], "rs", label="FFs", ms=5)



    labeled_gc = ["NGC 5139", "NGC 2419"]
    for n in range(len(harris.names)):
        if harris.names[n] in labeled_gc:
            print harris.names[n]
            # This is ridiculous
            if harris.names[n] == "NGC 5139":
                annotate(r"$\omega$ Cen", (harris.r_h[n], harris.M_V[n]),
                         xytext=(3,3), textcoords='offset points', fontsize=10)
            else:
                annotate(harris.names[n], (harris.r_h[n], harris.M_V[n]),
                         xytext=(3,3), textcoords='offset points', fontsize=8)
