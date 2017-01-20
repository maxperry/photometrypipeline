
# PhotoPipe

### Photometry Pipeline for RATIR and RIMAS

[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)
[![PyPI version](https://badge.fury.io/py/photopipe.svg)](https://badge.fury.io/py/photopipe)

PhotoPipe is a pipeline for automated reduction, photometry and astrometry, written in Python (2.7), and designed to process imaging data from the following instruments:

* [RATIR](http://butler.lab.asu.edu/RATIR/): the Reionization and Transients Infrared/Optical Project.
* [RIMAS](https://lowell.edu/research/research-facilities/4-3-meter-dct/rimas/): the Rapid infrared IMAger Spectrometer.

This project is based on the previous versions of the pipeline by [cenko](https://github.com/cenko/RATIR-GSFC) and [vickitoy](https://github.com/vickitoy/photometry_pipeline).  
Built at the NASA Goddard Space Flight Center, in collaboration with the University of Maryland.
<br><br><br>
![NASA GSFC - RATIR - UMD](https://github.com/maxperry/photometrypipeline/raw/master/docs/readme-logos.jpg)


## Full Documentation

See the [Wiki](https://github.com/maxperry/photometrypipeline/wiki) for full documentation, examples, operational details and other information.


## Prerequisites

The pipeline can be installed with very little additional software. However, depending on the installation, different libraries and software may be required. Please see the [Installation](#installation) section below, and the dedicated [wiki page](https://github.com/maxperry/photometrypipeline/wiki/Prerequisites) for a full list of prerequisites.


## Installation

The easiest way to get up and running with the pipeline is to download the ready-to-use **virtual machine** box and run it with Vagrant.

Alternatively, it can be installed on Linux or macOS either with `pip`, or with the provided installation scripts.
_(**WARNING:** compiling the dependencies can take more than 6 hours.)_

#### 1) Download a pre-configured Virtual Machine (PhotoPipe-VM)

Please refer to the virtual machine [repository](https://github.com/maxperry/photometrypipeline-vm).

#### 2) Install on your machine from PyPI or git

* Run `sudo -H pip install photopipe` to install the latest stable version from [PyPI](https://pypi.python.org/pypi/photopipe). 

* Or, clone from `git`:

 ```
 $ git clone git@github.com:maxperry/photometrypipeline.git
 $ cd photometrypipeline
 $ sudo python setup.py install
 ```
 
**NOTE (macOS)**: If the installation fails with `sudo: port: command not found` make sure that [MacPorts](https://guide.macports.org/#installing) is installed and `/opt/local/bin` is in the $PATH (e.g. `export PATH=/opt/local/bin:/opt/local/sbin:$PATH`).

#### 3) Manual Installation
If you to run the pipeline from the Python enviroment rather than using the `photopipe` command as described in the [Usage](#usage) section, please follow the step by step instructions in the wiki pages below to install all the dependencies manually.

* [Instructions for macOS](https://github.com/maxperry/photometrypipeline/wiki/Manual-Installation-(macOS))
* [Instructions for Linux-Debian](https://github.com/maxperry/photometrypipeline/wiki/Manual-Installation-(Linux-Debian))


## Usage

The following steps can be reproduced using the test data downloadable [here](https://drive.google.com/file/d/0BzMOBEOpFL9LaHpkWnFXc0IzRmM/view?usp=sharing), either from your host machine or from the virtual machine.

####1. Create a new directory with the following structure:
 
 **NOTE for PhotoPipe-VM**: Create the `imdata` directory in the VM's shared data dir (should be `./photopipe-vm/data`) 
 
 ```
imdata  
│
└───bias
│   │   20160628T032914C0b.fits
│   │   20160628T032914C1b.fits
│   │   ...
│   
└───dark
│   │   20160628T040211C0d.fits
│   │   20160628T040211C1d.fits
│   │   ...   
│
└───flat
│   │   20160628T024207C0f.fits
│   │   20160628T024207C1f.fits
│   │   ...    
│
└───science
│   │   20160628T043940C0o.fits
│   │   20160628T043940C1o.fits
│   │   ...    
│
└───science_selected
│
└───reduced
```
 - **If running on your Host Machine**: Open the terminal and do `cd ./imdata/reduced`
 - **If running on PhotoPipe-VM**: Launch the VM first ([see instructions](https://github.com/maxperry/photometrypipeline-vm#usage)), and `cd /vagrant_data/imdata/reduced`
 - **NOTE**: It's important to perform the steps below from the `reduced` folder, otherwise make sure to move the master frames there after running `mkmaster` (it saves to the current dir!). 

####2. Run preprocessing functions
 1. Enter Python environment: `$ python`
 2. Run the following script:
  ```python
  from photopipe.reduction import preproc
  
  # Bias frames calibration 
  
  bias_calib = preproc.choose_calib('ratir', 
                                    'bias', 
                                    workdir='/vagrant_data/imdata/bias/', 
                                    cams=[0,1], 
                                    auto=True, 
                                    amin=0.0, amax=1.0, 
                                    reject_sat=False, 
                                    save_select=True, 
                                    noplot=False)
 
 
  # Dark frames calibration
  
  dark_calib = preproc.choose_calib('ratir', 
                                    'dark', 
                                    workdir='/vagrant_data/imdata/dark/', 
                                    cams=[0,1], 
                                    auto=True, 
                                    amin=0.0, amax=1.0, 
                                    reject_sat=False, 
                                    save_select=True, 
                                    noplot=False)

  # Flat frames calibration 
  
  flat_calib = preproc.choose_calib('ratir', 
                                    'flat', 
                                    workdir='/vagrant_data/imdata/flat/', 
                                    cams=[0,1,2,3], 
                                    auto=True, 
                                    amin=0.2, amax=0.8, 
                                    reject_sat=False, 
                                    save_select=True, 
                                    noplot=False)
                                    
  # WARNING: RATIR flats are often bad even when 
  # the median value is in the acceptable range. 
  # Auto mode is only recommended for bias frame 
  # selection.     
  
  
  # Select science frames
  # (selected frames will be copied to target_dir)
  
  science_dict = preproc.choose_science('ratir', 
                                        workdir='/vagrant_data/imdata/science, 
                                        targetdir='/vagrant_data/imdata/science_selected', 
                                        cams=[0,1,2,3], 
                                        auto=True, 
                                        save_select=True, 
                                        calibrate=False, 
                                        noplot=False) 

  # WARNING: When auto is True, all science frames
  # are selected. Since the telescope occasionally
  # has tracking issues, it is recommended to check
  # all frames.  
  
  
  # Make master frames
  # (saves to the current dir)
  
  preproc.mkmaster('ratir', bias_calib, 'bias')

  preproc.mkmaster('ratir', dark_calib, 'dark')
 
  preproc.mkmaster('ratir', flat_calib, 'flat')  
 ```
 
 **WARNING**:
  1. Always use absolute paths terminated by a `/` when setting `dir` parameters. 
  2. Make sure to set the correct number of **cameras** where required (e.g. `cams[0,1,2,3`).
  3. `amin/amax` set the min and max **saturation** values. Make sure to select the appropriate values for each type of frame, only frames with median values in this range will be selected.
 
See [wiki page](https://github.com/maxperry/photometrypipeline/wiki/preproc.py) for `preproc` functions reference.
  
####3. Run auto reduction
 1. Start a new python environment
 2. Execute the script below:
  ```python
   from photopipe.reduction.auto.autoproc import autoproc

   autoproc(datadir='/vagrant_data/imdata/science_selected/', 
            imdir='/vagrant_data/imdata/reduced/', 
            redo=1)
  ```
  
See [wiki page](https://github.com/maxperry/photometrypipeline/wiki/autoproc.py) for `autoproc` function reference.
  
## Bugs and Feedback

For bugs, questions and discussions please use the [Github Issues](https://github.com/maxperry/photometrypipeline/issues).


## License
This project is licensed under the GNU GPLv3 License - see the [LICENSE](https://github.com/scrapy/scrapy/blob/master/LICENSE) file for details.
