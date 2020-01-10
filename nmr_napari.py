#!/usr/bin/env python
# -*- coding: utf-8 -*-#

import os
import sys
import napari
import nmrglue as ng
import skimage

"""This module contains functions to load and display NMR spectra using nmrglue
and napari, along with other utilities.

When run as a script:
Arguments:
    exp_path {str} -- path to experiment folder

Returns:
    viewer {qt instance} -- napari viewer object
"""


def load_pipe(exp_path):
    """Load nmrpipe formatted experiment

    Arguments:
        exp_path {str} -- path to experiment folder

    Returns:
        specpars {dict} -- experimental parameters
        spectrum {np.array} -- nD array containing spectrum intensities
    """
    specpars, spectrum = ng.fileio.pipe.read(exp_path)
    shape = spectrum.shape
    print('Loaded {}D spectrum sized {}'.format(len(shape), shape))

    return specpars, spectrum


def rescale(spectrum, scale=(0, 100)):
    """Rescale spectrum intensity to fit a defined scale

    Arguments:
        spectrum {np.ndarray} -- input spectrum

    Keyword Arguments:
        scale {tuple} -- desired output scale (default: {(0,100)})

    Returns:
        spectrum {np.array} -- rescaled spectrum
    """
    init_vals = (spectrum.min(), spectrum.max())
    print('Min intensity={:.2e}\tMax intensity={:.2e}'.format(*init_vals))
    print('Re-scaled intensities to {} - {}'.format(*scale))
    spectrum = skimage.exposure.rescale_intensity(spectrum, out_range=(0, 100))

    return spectrum


def calc_hist(spectrum):
    """Calculate the histogram of a spectrum

    Arguments:
        spectrum {np.ndarray} -- [description]

    Returns:
        counts {np.ndarray} -- Intensity value occurrences
        bins {np.ndarray} -- Center positions of computed bins
    """
    counts, bins = skimage.exposure.histogram(spectrum)

    return counts, bins


def calc_noise(counts, bins):
    """Estimate noise from most repeated value in histogram

    Arguments:
        counts {np.ndarray} -- Intensity value occurrences
        bins {np.ndarray} -- Center positions of computed bins

    Returns:
        noise {float} -- Estimated noise level
    """
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
        specpars {dict} -- experimental parameters
        spectrum {np.array} -- nD array containing spectrum intensities
    """
    # Load spectrum, rescale, and estimate noise
    name = os.path.basename(exp_path)
    specpars, spectrum = load_pipe(exp_path)
    spectrum = rescale(spectrum)
    counts, bins = calc_hist(spectrum)
    noise = calc_noise(counts, bins)
    # Define plotting parameters
    minlev = noise
    maxlev = spectrum.max()  # To be refined
    gamma = 1
    print('Min contour={:.2f}\tMax contour={:.2f}\tgamma={:.2f}'.format(minlev,
                                                                        maxlev,
                                                                        gamma))
    # Create viewer and add spectrum layer
    viewer = napari.Viewer()
    viewer.add_image(spectrum,
                     colormap='turbo',
                     name=name,
                     contrast_limits=(minlev, maxlev),
                     gamma=gamma)

    return viewer, specpars, spectrum


if __name__ == "__main__":
    exp_path = sys.argv[1]
    viewer, specpars, spectrum = view_spectrum(exp_path)
