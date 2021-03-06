# LightPipes

**Simulations of optical phenomena where diffraction is essential**

LightPipes is a set of functions written in Python (Before version 2.0.0 these functions are in C++). It is designed to model coherent optical devices when the diffraction is essential. We put the C++ based version of LightPipes in another repository: opticspy/clightpipes.
The pure Python version is as fast as the C++ version due to the use of the numpy, scipy and pyFFTW packages.

The toolbox consists of a number of functions. Each function represents an optical element or a step in the light propagation. There are apertures, intensity filters, beam-splitters, lenses and models of free space diffraction. There are also more advanced tools for manipulating the phase and amplitude of the light. The program operates on a large data structure, containing square two-dimensional arrays of complex amplitudes of the optical field of the propagating light beam.

The LightPipes routines are modifications of the LightPipes C routines written by Gleb Vdovin for Unix, Linux, DOS and OS2 workstations.

Visit the website of **Flexible Optical**: [http://www.okotech.com](http://www.okotech.com), where you can find the source code of *LightPipes* and a *manual*.

## Credits.

**Gleb Vdovin**

LightPipes was written by Gleb Vdovin in 1993 for MS DOS.
the output of a command was the input for the next command. Gleb used the "pipe" feature of MS DOS, that's why it is called LightPipes. The source code of LightPipes for UNIX can be obtained for free from http://www.okotech.com.

**Fred van Goor**

Fred rewrote the commands for using them in Mathcad and Matlab in 1996. Later he made a version for Python (2017).
The [Mathcad](http://www.ptc.com/engineering-math-software/mathcad) and [Matlab](https://www.mathworks.com/) versions (not free) can be obtained from [Flexible Optical](http://www.okotech.com). 
The Python version can be obtained for free from [https://github.com/opticspy/lightpipes](https://github.com/opticspy/lightpipes).
Fred also developed the [LightPipes for Python website](https://opticspy.github.io/lightpipes/) using [Sphinx](http://www.sphinx-doc.org).
Fred used LightPipes for his course "Introduction to Optics" at the [University of Twente](https://www.utwente.nl/en/education/bachelor/programmes/applied-physics/)
as an introduction to the optics lab.
See how Fred is [blowing out a HeNe laser](https://youtu.be/sm7MrA8Usuw?t=58s) during this lesson.
(The laser stops because the outcoupling mirror was covered with vapor from Fred's breath)

**Guyskk**

Guyskk' s contribution to the development of the LightPipes C++ package was very important.
His knowledge of Python helped a lot to get the package operative for the windows,
macintosh and several linux platforms.

**Leonard Doyle**

Leonard translated the C++ code of LightPipes into pure Python, using the numpy,
scipy and pyFFTW packages.
Thanks to the matrix routines of numpy and the fast pyFFTW Fourier transform
the speed of the pure python version is as fast as the C++ version.
Also parallel processing is possible now, which enhances the speed even more.
From LightPipes version 2.0.0 LightpIpes is pure python. As usual, it can be installed with pip:

    pip install LightPipes

The older C++ versions are still on [PyPi](https://pypi.python.org/pypi/LightPipes/) 
and can be installed by typing (for the latest C++ version, 1.2.0):

    pip install LightPipes==1.2.0

## Install

The pure Python version of LightPipes (from version 2.0.0) can be installed on any platform with Python version 3.+ installed. We tested it on a number of computers: Windows, Macintosh and Linux.

     


The packages are on [PyPi](https://pypi.python.org/pypi/LightPipes/), so simply open a terminal window and type at the prompt:

	pip install LightPipes

The C++ versions are still on [PyPi](https://pypi.python.org/pypi/LightPipes/) and can be installed using:

	pip install LightPipes==1.2.0

for the latest C++ version.

Use pip3 to install for Python 3.7.
You can also download packages from [Releases](https://github.com/opticspy/lightpipes/releases).

**Raspberry Pi:**

We encountered problems when installing it on a Raspberry Pi 4.0. It seems that pyFFTW package is not compatible with the Raspberry (ARM processor). Maybe they will solve that in the future. 
In the mean time you can follow the instructions in the [install section](https://opticspy.github.io/lightpipes/install.html#known-installation-problems) of our documentation to install pyFFTW on a Raspberry PI, or install the latest C++ version, 1.2.0, of LightPipes when Python 3.7 is installed on your Raspberry Pi.
Type at a terminal prompt:
    
	sudo pip3 install LightPipes==1.2.0

## Document

[https://opticspy.github.io/lightpipes/](https://opticspy.github.io/lightpipes/)

## Example: Young interferometer

A plane wave is diffracted by two small holes, separated a distance, d. So two more or less spherical waves will propagate from these holes.

![](img/twoholesSetUp.png)

The resulting interference pattern on a screen at distance z looks like:

![](img/twoholesPattern.png)

The Python program [Young.py](Examples/Interference/Young.py) described in detail.

The first step in *Python* is to import the *LightPipes* library:

```python
from LightPipes import *
```
Besides the *LightPipes* library, we import units and a few more in this way.

If the *LightPipes* library is successful installed on your computer Python can proceed with the next step.
You probably want to plot the results, so import *matplotlib*:

```python
import matplotlib.pyplot as plt
```

Next we define some variables: a wavelength of 20 micrometer , a 30 x 30 mm2 square grid with 500 x 500 pixels.

```python
wavelength = 20*um
size = 30.0*mm
N = 500
```

Now we are ready to start the simulation. The *Begin* command generates a field with amplitude 1.0 and phase zero, a plane wave. So, all the 500 x 500 elements of array, F, contain the complex number: 1.0 + j0.0.
The next commands generate two waves, F1 and F2, which are apertured by the two circular apertures and combined (simply added) by the *BeamMix* command. The combined wavefront is propagated a distance z=30 cm by the *Fresnel* command. After that the intensity is caculated and normalized to 255 (2 -> 255, 1 -> 1.0, 0 -> not normalized) by the *Intensity* command.

```python
F = Begin(size,wavelength,N)
F1 = CircAperture(0.15*mm, -0.6*mm,0, F)
F2 = CircAperture(0.15*mm, 0.6*mm,0, F)    
F = BeamMix(F1,F2)
F = Fresnel(10*cm,F)
I = Intensity(2,F)
```

The result is plotted using the fantastic *matplot* routines. We are not interested in axis around the pattern and we like to write a title above the plot.

```python
plt.imshow(I, cmap='rainbow');
plt.axis('off');
plt.title('intensity pattern')
plt.show()
```
