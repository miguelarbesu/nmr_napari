#!/usr/bin/env python
# -*- coding: utf-8 -*-#

import re
import os
import sys
import napari
import numpy as np
import nmrglue as ng
import matplotlib.pyplot as plt
from skimage import feature, exposure

"""This module contains functions to load and display NMR spectra using nmrglue
and napari, along with other utilities.

When run as a script:
Arguments:
    exp_path {str} -- path to experiment folder

Returns:
    viewer {qt instance} -- napari viewer object
"""


def natural_sort(l):
    """Sort list naturally

    Arguments:
        l {list} -- List of items.

    Returns:
        {list} -- List sorted naturally.
    """
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


def load_pipe(exp_path, pseudo=False):
    """Load nmrpipe formatted experiment

    Arguments:
        exp_path {str} -- path to experiment folder

    Returns:
        parameters {dict} -- experimental parameters
        data {np.array} -- nD array containing intensities
    """
    if pseudo is True:
        ftdirs = natural_sort(os.listdir(exp_path))
        totalfts = len(ftdirs)
        datalist = []
        i = 1
        for ft in ftdirs:
            ftpath = os.path.join(exp_path, ft, 'test%03d.dat')
            parameters, data = ng.fileio.pipe.read(ftpath)
            print('Loaded spectrum {} of {}'.format(i, totalfts))
            datalist.append(data)
            i += 1
        data = np.stack(datalist)
    else:
        parameters, data = ng.fileio.pipe.read(exp_path)
    print('Loaded {}D spectrum sized {}'.format(data.ndim, data.shape))

    return parameters, data


def rescale(data, scale=(0, 100)):
    """Rescale spectrum intensity to fit a defined scale

    Arguments:
        data {np.ndarray} -- input spectrum

    Keyword Arguments:
        scale {tuple} -- desired output scale (default: {(0,100)})

    Returns:
        data {np.array} -- rescaled spectrum
    """
    init_vals = (data.min(), data.max())
    print('Min intensity={:.2e}\tMax intensity={:.2e}'.format(*init_vals))
    print('Re-scaled intensities to {} - {}'.format(*scale))
    data = exposure.rescale_intensity(data, out_range=(0, 100))

    return data


def calc_hist(data):
    """Calculate the histogram of a spectrum

    Arguments:
        data {np.array} -- nD array containing intensities

    Returns:
        counts {np.ndarray} -- Intensity value occurrences
        bins {np.ndarray} -- Center positions of computed bins
    """
    counts, bins = exposure.histogram(data)

    return counts, bins


def plot_hist(data):
    """Plot spectrum histogram

    Arguments:
        data {np.array} -- nD array containing intensities
    """
    counts, bins = calc_hist(data)
    plt.plot(bins, counts)
    plt.show()


def calc_noise(data):
    """Estimate noise from most repeated value in histogram

    Arguments:
        counts {np.ndarray} -- Intensity value occurrences
        bins {np.ndarray} -- Center positions of computed bins

    Returns:
        noise {float} -- Estimated noise level
    """
    counts, bins = calc_hist(data)
    noise = bins[counts.argmax()]
    print('Noise level estimated at {:.2f}'.format(noise))

    return noise


def view_spectrum(exp_path):
    """Load and visualize a spectrum using napari.
    Requires running from an IPython session with qt backend:
    `ipython --gui qt`
    or, inside IPython, use magic
    `%gui qt`

    Arguments:
        exp_path {str} -- Path to experiment file

    Returns:
        viewer {qt instance} -- napari viewer object
        parameters {dict} -- experimental parameters
        data {np.ndarray} -- nD array containing spectrum intensities
    """
    # Load spectrum, rescale, and estimate noise
    name = os.path.basename(exp_path)
    parameters, data = load_pipe(exp_path)
    data = rescale(data)
    noise = calc_noise(data)
    # Define plotting parameters
    minlev = noise
    maxlev = data.max()  # To be refined
    gamma = 1
    print('Min contour={:.2f}\tMax contour={:.2f}\tgamma={:.2f}'.format(minlev,
                                                                        maxlev,
                                                                        gamma))
    # Create viewer and add spectrum layer
    viewer = napari.Viewer()
    viewer.add_image(data,
                     colormap='turbo',
                     name=name,
                     contrast_limits=(minlev, maxlev),
                     gamma=gamma)

    return viewer, parameters, data


def pick_peaks(data):
    """Pick peaks from a spectrum

    Arguments:
        spectrum {np.array} -- nD array containing intensities

    Returns:
        peaks {np.ndarray} -- peak coordinates
    """
    noise = calc_noise(data)
    # NOTE: The peak picking parameters are harcoded for now assuming:
    # 1. That the noise distribution is normal, centered at the estimated value
    # 2. The average lenght of a protein construct is 120 aa
    threshold = 2*noise
    max_peaks = 120
    peaks = feature.peak_local_max(data,
                                   threshold_abs=threshold,
                                   num_peaks=max_peaks,
                                   min_distance=1)
    print('{} peaks found with minimum threshold at {:.2f}'.format(len(peaks),
                                                                   threshold))
    return peaks


def draw_peaks(peaks, viewer):
    """Draws peak list

    Arguments:
        peaks {np.ndarray} -- peak coordinates
        viewer {qt instance} -- napari viewer object
    """
    viewer.add_points(peaks,
                      name='peaks',
                      size=4,
                      symbol='x',
                      n_dimensional=True)
    print('{} peaks drawn'.format(len(peaks)))


class Spectrum():
    def __init__(self, parameters, data, name='spectrum'):
        self.parameters = parameters
        self.data = data
        self.name = name
    pass


if __name__ == "__main__":
    exp_path = sys.argv[1]
    viewer, parameters, data = view_spectrum(exp_path)
