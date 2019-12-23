#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import napari
import numpy as np
import nmrglue as ng

expn = sys.argv[1]
dic, data = ng.fileio.pipe.read(expn)
print('Loaded spectrum with dimensions {}'.format(data.shape))

noise = np.median(data)
maxint = data.max()

minlev = (1e3)*noise
maxlev = 0.1*maxint
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
