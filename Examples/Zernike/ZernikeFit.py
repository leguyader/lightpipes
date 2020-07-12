# -*- coding: utf-8 -*-
"""
Example script to demonstrate Zernike decomposition of a beam input field.

First, we define a circular beam of uniform (flat-top) intensity and add
abberations on purpose by applying Zernike modes to the phase.
Next, the decomposition is done which unwraps the phase of the field to
generate a wavefront and returns the magnitude of each Zernike mode 
contribution in the chosen units.

If norm=False, no normalization is applied and the specified amplitude is
the maximum amplitude, i.e. peak to valley or half peak-to-valley for
polynomials extending to plus and minus.
If norm=True, each contribution is normalized to the same rms deviation and
therefore the specified amplitudes correspond to rms deviations.

The units can be changed between
'lam' for multiples of lambda (e.g. 0.5, 1)
'opd' for optical path distance (e.g. 800*nm for 1 lam) or
'rad' for phase differences (e.g. 2*np.pi for 1 lam).
The default is OPD, so make sure to supply a small number for A (e.g. 800*nm
or 800e-9, not 1)

The fit only works reliably on a circular input beam whose unwrapped phase is
smooth (no 2pi jumps). If there are jumps in the unwrapped phase, the
results will not be meaningful.
Also, if the input beam is not circular, the Zernike polynomials do not form
an orthogonal set anymore, and the magnitude for 1 mode may vary with the
number of modes fitted! In that case, it is best to fit more modes and see
if the decomposition "converges".
So far, no intensity weighting is applied for e.g. Gaussian beams, but the
wavefront is fitted to Zernike modes regardless of intensity distribution.

Created on Sun Jul 12 12:39:03 2020

@author: Leonard.Doyle
"""

import numpy as np
import matplotlib.pyplot as plt

import LightPipes as lp
from LightPipes import noll_to_zern
import LightPipes.plotutils as lpplot
from LightPipes.units import m, cm, mm, um, nm


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
j_fits = j_fits[1:] #do not plot piston, meaningless
A_fits = A_fits[1:]/um #convert [m] to [um]
#with norm=True and units=lam, A are rms values in lambdas
#with norm=True and units=opd, A are rms values in [um]
#with norm=False and units=opd, A are PtV values in [um]
plt.bar(j_fits, A_fits)
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
