"""
    Purpose: repository for general purpose functions used by the preprocessing pipeline.  
    Usage:      
        *)  most of these functions are usually called by other scripts, not by the user
    Notes:
        - is numpy mean() buggy? found cases where mean(array) < min(array)
        - added show_list() to aid users in reviewing results
"""

# standard modules
import os
import itertools

# installed modules
import astropy.io.fits as pf
import matplotlib.pylab as pl
import numpy as np
from scipy.ndimage.interpolation import zoom

# custom modules/functions
from zscale import zscale

"""
ANSI escape sequences to print to terminal with color
http://stackoverflow.com/questions/287871/print-in-terminal-with-colors-using-python
"""
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def print_blue(msg):
    print bcolors.OKBLUE + msg + bcolors.ENDC

def print_bold(msg):
    print bcolors.BOLD + msg + bcolors.ENDC

def print_under(msg):
    print bcolors.UNDERLINE + msg + bcolors.ENDC

def print_head(msg):
    print bcolors.HEADER + msg + bcolors.ENDC

def print_warn(msg):
    print bcolors.WARNING + msg + bcolors.ENDC

def print_err(msg):
    print bcolors.FAIL + msg + bcolors.ENDC

def imcombine(indata, type='median', ret_std=False):

	"""
	Written by John Capone (jicapone@astro.umd.edu).
	
	Purpose:        combine stack of frames
	
	Input:
		indata:     stack of frames to be combined
		type:       function used to combine the stack.  currently only mean or median
		ret_std:    if set to true, return the standard deviation of each pixel.  default is false
		
	Output:
		combined:   combined stack
		sigma:      standard deviation of each pixel (optional)
		
	Notes:
		- added normalization before combining
		
	Future Improvements:
		- add outlier rejection
	
	"""
	
	if indata.ndim != 3:
		print_warn("Warning: data should be 3D stack of frames.")
	
	if type is 'mean':
		combined = np.mean(indata, axis=0)	
	else:
		combined = np.median(indata, axis=0)
		
	if ret_std:
		sigma = np.std(indata, axis=0)
		return combined, sigma
	else:
		return combined

def robust_sigma(y, zero=False):

	"""
	Converted from IDL ROBUST_SIGMA function by John Capone (jicapone@astro.umd.edu).
	
	Purpose:  Calculate a resistant estimate of the dispersion of a distribution. 
		For an uncontaminated distribution, this is identical to the standard deviation.

    Input:
        y:    Vector of quantity for which the dispersion is to be calculated
        zero: if set, the dispersion is calculated w.r.t. 0.0 rather than the central value of the vector. 
        	If Y is a vector of residuals, this should be set.
    
    Output:   robust_sigma returns the dispersion. In case of failure, returns value of -1.0

    Notes:
        - 
	"""
	
	eps = 1.0e-20
	if zero:
		y0 = 0.
	else:
		y0 = np.median(y)
		
	# first, the median absolute deviation about the median:
	mad = np.median(np.abs(y-y0))/.6745
	# if the mad=0, try the mean absolute deviation:
	if mad < eps:
		mad = np.average(np.abs(y-y0))/.80
	if mad < eps:
		sigma = 0.
		return sigma
	
	# now the biweighted value:
	u   = (y-y0)/(6.*mad)
	uu  = u*u
	q = uu <= 1.0
	count = np.sum(q)
	if count < 3:
		print_err('Error: robust_sigma - this distribution is too weird! returning -1.')
		sigma = -1.
		return sigma
	numerator = np.sum((y[q] - y0)**2. * (1 - uu[q])**4.)
	n = y.size
	den1 = np.sum((1. - uu[q]) * (1. - 5.*uu[q]))
	sigma = n * numerator / (den1 * (den1 - 1.))
	if sigma > 0.:
		sigma = np.sqrt(sigma)
	else:
		sigma = 0.
	return sigma

def show_list(fits_fns, nx=5, ny=3, size_mult=3.2, zoom_lvl=None, fontsize=8):

	"""
    Written by John Capone (jicapone@astro.umd.edu).

    Purpose:        displays images in specified list file.

    Input:
        fits_fns:   list of file names
        nx:         number of images to display simultaneously in x
        ny:         number of images to display simultaneously in y
        size_mult:  multiple to determine image sizes
        zoom_lvl:   amount to decrease image resolution by (to save memory)
        fontsize:   fontsize for plots

    Usage:
        1)  enter python or ipython environment
        2)  load function -> 'from rat_preproc import show_list'
        3)  run function -> 'show_list(fits_fns)'
            - decreasing zoom_lvl (i.e. from 0.5 to 0.1) decreases the size of the displayed image, thus decreasing the amount of memory required
        4)  function will display arrays of images in list file for inspection by user

    Notes:
        - added escape character
        - user can now change font size
	"""
	
	nx = int(nx); ny = int(ny) # force parameter types to int
	
	if type(fits_fns) not in [list, dict]:
		print_err("Error: fits_fns must be a list or dictionary of fits file names. Exiting...")
		return
	if type(fits_fns) is dict: # convert dictionary to list
		fits_fns = list(itertools.chain.from_iterable(fits_fns.values()))
		
	nfits = len(fits_fns) # number of fits files listed
	
	pl.ion() # pylab in interactive mode
	
	# create figures of subplots for review
	nfigs = int(np.ceil(nfits / float(nx * ny)))
	fig = pl.figure(figsize=(nx*size_mult,ny*size_mult), tight_layout=True)
	for i in range(nfigs):
	
		start_fits = i*nx*ny
		if (i + 1)*nx*ny <= nfits:
			stop_fits = (i + 1)*nx*ny - 1
			nsubplts = nx*ny
		else:
			stop_fits = nfits
			nsubplts = nfits - start_fits
			
		print "Displaying frames {} - {} ({} total).".format(start_fits, stop_fits, nfits)
		
		# display image in each subplot
		for j in range(nsubplts):
			
			ax = fig.add_subplot(ny, nx, j+1) # new subplot
			
			fits_fn = fits_fns[start_fits + j]
			dpath, fits_fn = os.path.split(fits_fn)
			if len(dpath) != 0:
				dpath += '/'
			temp = fits_fn.split('.')
			if len(temp) == 1:
				fits_fn += '.fits'
			elif len(temp) == 2:
				if temp[1].lower() != 'fits':
					print_err("Error: invalid \"{}\" file type detected.  Should be \"fits\" file type. Exiting...".format(temp[1]))
					pl.close('all') # close image to free memory
					return
			else:
				print_err("Error: file names should not include \".\" except for file type extention.  Exiting...")
				pl.close('all') # close image to free memory
				return
			fits_id = temp[0] # fits file name with extention removed
			
			# open data
			hdulist = pf.open(dpath + fits_fn)
			im = hdulist[0].data
			h = hdulist[0].header
			hdulist.close() # close FITs file
			
			# display
			
			if zoom_lvl is not None:
				imdisp = zoom(im, zoom_lvl)	
			else:
				imdisp = im
			z1, z2 = zscale(im)
			ax.imshow(imdisp, vmin=z1, vmax=z2, origin='lower', cmap=pl.cm.gray, interpolation='none')
			ax.set_xticks([])
			ax.set_yticks([])
			
			try:
				ax.set_title("{} - {} filter".format(fits_id, h['FILTER']), fontsize=fontsize) # title with identifier
			except:
				ax.set_title("{}".format(fits_id), fontsize=fontsize) # title with identifier
			
			
		fig.canvas.draw()
		usr_select = raw_input("Press any key to continue or \"q\" to quit: ") # prompt user to continue
		fig.clear() # clear image
		
		if usr_select.lower() == 'q':
			print_bold("Exiting...")
			pl.close('all') # close image to free memory
			return
			
	pl.close('all') # close image to free memory

def zsview(im, cmap=pl.cm.gray, figsize=(8,5), contours=False, ccolor='r'):
    z1, z2 = zscale(im)
    pl.figure(figsize=figsize)
    pl.imshow(im, vmin=z1, vmax=z2, origin='lower', cmap=cmap, interpolation='none')
    if contours:
        pl.contour(im, levels=[z2], origin='lower', colors=ccolor)
    pl.tight_layout()
