# -*- coding: utf-8 -*-

import numpy as np
from numpy import pi

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from LightPipes import *

f = Begin( 20*mm, 1*um, 1000 );  # grid size is 20*mm, grid dimension is N=1000, wavelength is 1um
c = CircAperture( f, 3*mm );
tx = 26*mrad
c = Tilt( c, tx, 0*mrad );    # the planewave tilts 'tx' in x direction
I0 = Intensity(c)
Ph0 = Phase(c)
Ph0[I0<I0.max()*0.001] = np.nan #blank phase where intensity so low it's meaningless

# c = Fresnel( c, 100*mm );
c = Forvard( c, 100*mm );
I1 = Intensity( c )
Ph1 = Phase(c)
Ph1[I1<I1.max()*0.001] = np.nan

fig, axs = plt.subplots(2,2,
                        figsize=[7,8], #fiddle size to make nice PNG
                        sharex=True, sharey=True #zoom all simultaneously
                        )
[ax_tl, ax_tr], [ax_bl, ax_br] = axs #unpack tuple to access each axis

ax_tl.imshow(I0);
fig.suptitle( "size={}mm, lambda={}um, N={}, tx={}mrad".format(c.siz/mm,
                                                               c.lam/um,
                                                               c.N,
                                                               tx/mrad));

ax_tr.imshow(Ph0, vmin=-pi, vmax=pi)

ax_bl.imshow(I1)
ax_br.imshow(Ph1)

plt.show()

