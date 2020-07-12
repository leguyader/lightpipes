# -*- coding: utf-8 -*-
"""
Example script to demonstrate Zernike filtering of a beam input field.

First, several abberations are addded to a circular beam of uniform
(flat-top) intensity by applying Zernike modes to the phase.

Next, the decomposition is done (see ZernikeFit example) and the specified
modes are subtracted from the phase, yielding a "corrected" input beam with
lower abberations and thus better beam quality / focusability.

The fit and filter only works reliably on a circular input beam whose
unwrapped phase is smooth (no 2pi jumps). If there are jumps in the unwrapped
phase, the filtering will give bad results.
Also, if the input beam is not circular, the Zernike polynomials do not form
an orthogonal set anymore, making the fit ambiguous if not all modes present
in the beam are fitted for. In that case, more modes can be fitted than are
actually subtracted from the beam. It is best to fit a different number of
modes and see if the decomposition "converges".

Created on Sun Jul 12 12:49:36 2020

@author: Leonard.Doyle
"""

import numpy as np
import matplotlib.pyplot as plt

import LightPipes as lp
from LightPipes import noll_to_zern
import LightPipes.plotutils as lpplot
from LightPipes.units import m, cm, mm, um, nm, PI


lam = 800*nm      #[m] wavelength
N = 128

beam_D = 10*cm
beam_r = beam_D/2
L = beam_D #can add padding around if we plan to propagate etc.

F = lp.Begin(L, lam, N)
F = lp.CircAperture(F, beam_r) #new syntax of LightPipes v2
F = lp.Zernike(F, 2, 2, beam_r, units='lam', norm=False) #A=1 lambda -> 2lam PtV
F = lp.Zernike(F, *noll_to_zern(3), beam_r, units='lam') # syntax see (*)
F = lp.Zernike(F, *noll_to_zern(9), beam_r, A=0.3, units='lam') #0.3lam
# (*) = unpack tuple returned by noll_to_zern to 2 args for n, m

fdpi = 100 #matplotlib default
figsiz = (1000/fdpi, 400/fdpi) #specify figsize in pixel to make plot nicer


plt.figure(figsize=figsiz)
j_fits, A_fits = lp.ZernikeFit(F, 10, beam_r, norm=True, units='opd')
plt.bar(j_fits, A_fits/um) #convert amplitudes from OPD[m] to [um]
plt.title('Zernike fit of input field')
plt.xlabel('Noll index of mode')
plt.ylabel('rms [um]')
plt.show()

ptv = lp.Wavefront_ptv(F)/um
rms = lp.Wavefront_rms(F)/um
strehl = lp.Strehl(F)

metrics = 'PtV:{:.3f}um ; rms:{:.3f}um ; Strehl:{:.3f}'.format(
    ptv, rms, strehl)
print(metrics)

lpplot.Plot(F, circ_r=beam_r, unwrap=True,
        title='Parallel beam: {}'.format(
             metrics), ph_units='lam',
         figsize=figsiz)

#***** Filtering of a certain mode ******
#supplying a number filters all modes up to this mode j_Noll:
# filter_which = 6  #modes 1-6
#*OR* supplying a tuple/list/array will filter only these:
filter_which = [9]  #only Noll mode 9 = Zernike (3, -3) = vertical trefoil
F = lp.ZernikeFilter(F, filter_which, beam_r)

"""CAUTION: if your beam is non-circular, filtering more modes than filtered
# changes the result since Zernikes do not form an orthogonal set. At
# sufficiently high number of fitted modes, it should converge, so play
# around and test."""
# F = lp.ZernikeFilter(F, filter_which, beam_r,
#                      j_terms_to_fit=40) #fit 40, filter only 9


#****** Repeat plot of beam and Zernike modes ****
plt.figure(figsize=figsiz)
j_fits, A_fits = lp.ZernikeFit(F, 10, beam_r, norm=True, units='opd')
plt.bar(j_fits, A_fits/um) #convert amplitudes from OPD[m] to [um]
plt.title('Zernike fit of filtered beam')
plt.xlabel('Noll index of mode')
plt.ylabel('rms [um]')
plt.show()

ptv = lp.Wavefront_ptv(F)/um
rms = lp.Wavefront_rms(F)/um
strehl = lp.Strehl(F)

metrics = 'PtV:{:.3f}um ; rms:{:.3f}um ; Strehl:{:.3f}'.format(
    ptv, rms, strehl)
print(metrics)

lpplot.Plot(F, circ_r=beam_r, unwrap=True,
        title='Parallel filtered beam: {}'.format(
             metrics), ph_units='lam',
         figsize=figsiz)