# nmr_napari

A collection of utilities for visualization and analysis of n-dimensional Nulear Magnetic ResonanceNMR data sets using [**napari**](https://www.napari.org/) and [**nmrglue**](https://www.nmrglue.com/)<sup>1</sup>.

---

This project is currently conceived as a single script where to collect useful functions that simplify loading, manipulation and analysis of NMR data using napari.

Since napari is under development, the code here is expected to break and be subject to change.

## Motivation and purpose

High dimensional (>=3) NMR spectra are nowadays common in biomolecular NMR.

However, these datasets are still difficult to explore with the most common software tools (e.g. Bruker Topspin, MestReNova, CCPNMR, Sparky, nmrDraw, etc.), which are based on 2D contour projections of the nD array. 3D rendering is, when available, clunky and see-only.

The recent appearance of the Python-based multidimensional image viewer **napari** permits interface directly with modern image analysis and machine learning tools. The **nmrglue** library can smoothly parse the most common formats of NMR data into **numpy nD arrays**, which can be then passed to **napari** and analyzed as multidimensional image stacks.

This approach enables fast, interactive analysis of high dimensional data sets, such as 4D or time-resolved 3D experiments. Image processing-bassed methods are tested here in order to perform of the typical tasks in NMR, such as noise estimation, peak picking, assignment, etc.

## Usage

1. (Optional) Create a Python virtual environment for nmr-napari.

If you use [`conda`](https://docs.conda.io/en/latest/), a dependency file (`env.yml`) is provided.

2. Install [napari](https://napari.org/tutorials/installation.html) and [nmrglue](https://nmrglue.readthedocs.io/en/latest/install.html).

napari is installable using `pip`. If you used the environment above, it has
been already installed.

3. Start an IPython session with Qt as GUI backend

Two options are possible:

- Start `ipython --gui qt`
- Start `ipython`, then use magic `%gui qt`

4. Import or run nmr_napari

All functions can be imported from `nmr_napari.py`. The documentation should make use evident for now.

For quick visualization of a nmrPipe file (see next section for more details), `nmr_napari.py` can be run as a script:

```python
run nmr-napari.py PATH_TO_EXPFILE
```

This will load and display the selected experiment file.

Alternatively, `import nmr_napari` and use the `nmr_napari.load_pipe()` or
`nmr_napari.load_bruker()` functions to quickly load and inspect data sets.

Both modes return the following objects, that can be used to progammatically
interact with the data set:

- `viewer`: A napari.Viewer instance.
- `parameters`: A dictionary containing the NMR-speficic parameters of the spectrum.
- `data`: A n-dimensional array of intensity values.

Example:

```python
viewer, parameters, data = nmr_napari.load_bruker(path_to_expno, ndim=4)
```


## Current limitations and next steps

Currently, nmrPipe and Bruker binary file formats are implemented.

Dimensionality is parsed automatically for nmrPipe data sets, but has to be
provided for Bruker experiments via the ndim keyword argument. 3 and 4
dimensions are supported, but extension is trivial.

Time-resolved experiments can be read as collections of nmrPipe `tes%03d.dat`
planes stored in consecutive `ft*/` folders.

At this stage, spectra are treated as arrays of pixels containing intensities. Referencing the dimensions to the corresponding chemical shift scales is a desirable feature under development.

## References

<sup>1</sup> J.J. Helmus, C.P. Jaroniec, Nmrglue: An open source Python package for the analysis of multidimensional NMR data, J. Biomol. NMR 2013, 55, 355-367.