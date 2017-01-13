PhotoPipe
=========

Photometry Pipeline for RATIR and RIMAS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

|astropy|

PhotoPipe is a pipeline for automated reduction, photometry and
astrometry, written in Python (2.7), and designed to process imaging
data from the following instruments:

-  `RATIR`_: the Reionization and Transients Infrared/Optical Project.
-  `RIMAS`_: the Rapid infrared IMAger Spectrometer.

| This project is based on the previous versions of the pipeline by
  `cenko`_ and `vickitoy`_.
| Built at the NASA Goddard Space Flight Center, in collaboration with
  the University of Maryland.
| 
| |NASA GSFC - RATIR - UMD|

Full Documentation
------------------

See the `Wiki`_ for full documentation, examples, operational details
and other information.

Prerequisites
-------------

The pipeline can be installed with very little additional software.
However, depending on the installation, different libraries and software
may be required. Please see the `Installation`_ section below, and the
dedicated `wiki page`_ for a full list of prerequisites.

Installation
------------

The easiest way to get up and running with the pipeline is to download
the ready-to-use **virtual machine** box and run it with Vagrant.

| Alternatively, it can be installed on Linux or macOS either with
  ``pip``, or with the provided installation scripts.
| **WARNING**: *compiling the dependencies can take more than 6
  hours.*

1) Download a pre-configured Virtual Machine
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Please refer to the virtual machine `repository`_.

2) Install on your machine from PyPI or git
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  Run ``pip install --user photopipe`` to install the latest stable
   version from `PyPI`_.

-  Or clone from ``git``:

``$ git clone git@github.com:maxperry/photometrypipeline.git  $ cd photometrypipeline/  $ python setup.py install``

3) Manual Installation
^^^^^^^^^^^^^^^^^^^^^^

If you to run the pipeline from the Python enviroment rather than using
the ``photopipe`` command as described in the `Usage`_ section, please
follow the step by step instructions in the wiki pages below to install
all the dependencies manually.

-  `Instructions for macOS`_
-  `Instructions for Linux-Debian`_


Bugs and Feedback
-----------------

.. _RATIR: http://butler.lab.asu.edu/RATIR/
.. _RIMAS: https://lowell.edu/research/research-facilities/4-3-meter-dct/rimas/
.. _cenko: https://github.com/cenko/RATIR-GSFC
.. _vickitoy: https://github.com/vickitoy/photometry_pipeline
.. _Wiki: https://github.com/maxperry/photometrypipeline/wiki
.. _Installation: #installation
.. _wiki page: https://github.com/maxperry/photometrypipeline/wiki/Prerequisites
.. _repository: https://github.com/maxperry/photometrypipeline-vm
.. _PyPI: https://pypi.python.org/pypi/photopipe
.. _Usage: #usage
.. _Instructions for macOS: https://github.com/maxperry/photometrypipeline/wiki/Manual-Installation-(macOS)
.. _Instructions for Linux-Debian: https://github.com/maxperry/photometrypipeline/wiki/Manual-Installation-(Linux-Debian)
.. _sample data: 

.. |astropy| image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
   :target: http://www.astropy.org/
.. |NASA GSFC - RATIR - UMD| image:: https://github.com/maxperry/photometrypipeline/raw/master/docs/readme-logos.jpg