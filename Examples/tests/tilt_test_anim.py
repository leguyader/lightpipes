# -*- coding: utf-8 -*-
"""
Animation example utilising LightPipes and matplotlib.

The field propagation is calculated inside each update of the plot,
significantly limiting the update speed live on screen. If the refresh
rate should be faster, the Intensities etc. can be pre-calculated and saved
in arrays/lists.
In addition the 'interval' of the animation can be fine tuned.

A video of the animation can also be saved. The playback speed of the
finished video is completely independent of the refresh speed on screen
and can be set with 'fps' below.

Created on Sun Jun 21 14:56:44 2020

@author: Leonard.Doyle
"""


import numpy as np
from numpy import pi

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from LightPipes import *

def getfields(tx):
    """A complete simulation run for the given tilt, returning the
    intensity and phase before and after propagation, and a reference
    to the final field for further use."""
    f = Begin( 20*mm, 1*um, 1000 );  # grid size is 20*mm, grid dimension is N=1000, wavelength is 1um
    c = CircAperture( f, 3*mm );
    c = Tilt( c, tx, 0*mrad );    # the planewave tilts 'tx' in x direction
    I0 = Intensity(c)
    Ph0 = Phase(c)
    Ph0[I0<I0.max()*0.001] = np.nan #blank phase where intensity so low it's meaningless
    
    c = Fresnel( c, 100*mm );
    # c = Forvard( c, 100*mm );
    I1 = Intensity( c )
    Ph1 = Phase(c)
    Ph1[I1<I1.max()*0.001] = np.nan
    return I0, Ph0, I1, Ph1, c

# txs = np.arange(0.25,50,0.25)*mrad #simple range
# more customized range by mixing manual and automatic ranges:
txs = np.array([0.15,
                0.25,
                0.5,
                *range(1, 25),
                24.5,
                24.75,
                25,
                25.25,
                25.5,
                *range(26, 50),
                49.5,
                49.75,
                50,
                50.25,
                50.5,
                *range(51,70),
                ]) * mrad #scale each value by mrad once converted to numpy array

I0, Ph0, I1, Ph1, c = getfields(txs[0]) #init plots with first round,
    # then update in animation with txs[1,2,...]

fig, axs = plt.subplots(2,2,
                        figsize=[7,8], #fiddle size to make nice PNG
                        sharex=True, sharey=True #zoom all simultaneously
                        )
[ax_tl, ax_tr], [ax_bl, ax_br] = axs #unpack tuple to access each axis

im_tl = ax_tl.imshow(I0);
fig.suptitle( "size={}mm, lambda={}um, N={}, tx={}mrad".format(c.siz/mm,
                                                               c.lam/um,
                                                               c.N,
                                                               txs[0]/mrad));

im_tr = ax_tr.imshow(Ph0, vmin=-pi, vmax=pi)
im_bl = ax_bl.imshow(I1)
im_br = ax_br.imshow(Ph1)

plt.show()


def updatefig(ii):
    tx = txs[ii]
    I0, Ph0, I1, Ph1, c = getfields(tx)
    
    fig.suptitle("size={}mm, lambda={}um, N={}, tx={}mrad".format(c.siz/mm,
                                                               c.lam/um,
                                                               c.N,
                                                               tx/mrad));
    im_tl.set_data(I0)
    im_tr.set_data(Ph0)
    im_bl.set_data(I1)
    im_br.set_data(Ph1)
    return [im_tl, im_tr, im_bl, im_br]
        

ani = animation.FuncAnimation(fig, updatefig, frames = len(txs),
                              interval=500, blit=False, repeat=False)

plt.show()

savevideo = True
"""Note: it might be necessary to set up matplotlib and ffmpeg together.
There are tutorials and troubleshoots online, or maybe a different writer
can be selected which is available on your machine. Check matplotlib
documentation.
If it does not work, set savevideo=False to skip video writer
"""
if savevideo:
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=1.5, metadata=dict(artist='Me'), bitrate=4000)
    ani.save('tilt_anim_fresnel.gif', writer=writer)
