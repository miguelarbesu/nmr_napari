#!/usr/bin/env python

import sys
import napari
import numpy as np
import nmrglue as ng
import skimage


expn = sys.argv[1]
dic, data = ng.fileio.pipe.read(expn)
print('Loaded spectrum with dimensions {}'.format(data.shape))

# Rescale intensity to 0-100 range
data = skimage.exposure.rescale_intensity(data, out_range=(0,100))
# Determine noise from most common intensity in histogram
hist, bins = skimage.exposure.histogram(data)

noise = bins[hist.argmax()]
maxint = data.max()


minlev = noise
maxlev = maxint
gamma = 1

print('noise={:.2e}\tmax intensity={:.2e}'.format(noise, maxint))
print('min contour={:.2e}\tmax contour={:.2e}\tgamma={:.2f}'.format(minlev,
                                                                    maxlev,
                                                                    gamma))

v = napari.Viewer()
v.add_image(data,
            colormap='turbo',
            name='4D',
            contrast_limits=(minlev, maxlev),
            gamma=gamma)
            
