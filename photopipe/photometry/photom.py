"""
NAME:
	photom

PURPOSE:
	Creates a stacked image with all filters (saved to multicolor.fits and multicolor.weight.fits). 
	Crops all of the filter images and combined image to the same size and coordinates (file).ref.multi.fits.
	Then finds sources using the master combined image, and calculates the magnitude based on
	just the filtered images (with weight file using sextractor.  Saves new sextractor values to 
	"fluxes_(FILTER).txt'. Calculates absolute magnitude errors based on fluxes_(FILTER).txt and keyword 
	in file for absolute zeropoint RMS. Saves final magnitudes to 'finalphot(FILTER).am'

OUTPUTS:
	multicolor.[weights.]fits - files with all filter images stacked
	(file).ref.[weights.]fits - files resampled from swarp using all files
	(file).ref.multi.[weights.]fits - files cropped so that they include all filters
	coords(FILTER) - RA and DEC coordinates from sextractor from cropped images
	fluxes_(FILTER).txt - sextractor output from cropped images
	finalphot(FILTER).am - file containing pixel and coordinate location, mag and 
						   corrected mag error, flux and flux error of sextractor sources found (fluxes_(FILTER).txt)
	
Translated from icoords.pro by John Capone (jicapone@astro.umd.edu).
Modified by Vicki Toy (vtoy@astro.umd.edu) 5/21/2014
"""

import numpy as np
import os
import astropy.io.fits as pf
import photprocesslibrary as pplib
from string import index
from numpy import shape
from astropy import wcs


def photom(prefchar='coadd'):

	#Identify files (must have same number of images files as weight files)
	zffiles     = pplib.choosefiles(prefchar + '*_?.fits')
	weightfiles = pplib.choosefiles(prefchar + '*_?.weight.fits')
	
	if len(zffiles) > len(weightfiles):
		print 'Must have matching weight file to each image file to run automatic crop.'
		print 'To use manual crop user manualcrop keyword and change crop values by hand'
		return -1
		
	numfiles = len(zffiles)
		
	#Resample all images using SWarp to a reference image called multicolor using weight files
	swarpstr = ''
	for i in range(numfiles):		
		swarpstr = swarpstr + zffiles[i] + ' '

	stackcmd = 'swarp ' + swarpstr + '-DELETE_TMPFILES N -WRITE_XML N -SUBTRACT_BACK N -WEIGHT_TYPE MAP_WEIGHT -IMAGEOUT_NAME multicolor.fits -WEIGHTOUT_NAME multicolor.weight.fits'
	stackcmd = stackcmd + ' -COPY_KEYWORDS OBJECT,TARGNAME,TELESCOP,FILTER,INSTRUME,OBSERVAT,PIXSCALE,ORIGIN,CCD_TYPE,JD,DATE-OBS,DATE1,DATEN,AIRMASS,TOTALEXP,FLATFLD,FLATTYPE,SEEPIX,ABSZPT,ABSZPTSC,ABSZPRMS'
	print stackcmd
	os.system( stackcmd )

	#Rename all the resampled files to crop files
	for i in range(numfiles):
		tmp = zffiles[i].split('.')[0]
		ifile = tmp + '.resamp.fits'
		ofile = tmp + '.ref.fits'
		mvcmd = 'mv -f ' + ifile + ' ' + ofile
		os.system(mvcmd)
		
		ifile = tmp + '.resamp.weight.fits'
		ofile = tmp + '.ref.weight.fits'
		mvcmd = 'mv -f ' + ifile + ' ' + ofile
		os.system(mvcmd)

	coaddfiles = pplib.choosefiles(prefchar+'*_?.ref.fits')

	ra1arr  = []
	dec1arr = []
	ra2arr  = []
	dec2arr = []

	#Finds the RA and DEC of the first and the last pixel of each cropped coadded file
	for files in coaddfiles:

		fitsfile = pf.open(files)
		fitsheader = fitsfile[0].header
		data = fitsfile[0].data
	
		imSize = shape(data)
	
		#Converts pixel value to RA and DEC using header information (AstroPy function)
		w = wcs.WCS(fitsheader)
		pixcrd = [[0.,0.], [imSize[1]-1.0, imSize[0]-1.0]]
		[[ra1,dec1],[ra2,dec2]] = w.wcs_pix2world(pixcrd, 0)
	
		#Stores information into arrays
		ra1arr.append(ra1)
		dec1arr.append(dec1)
		ra2arr.append(ra2)
		dec2arr.append(dec2)

	#Finds the coordinates that fit all of the data
	raleft  = min(ra1arr)
	raright = max(ra2arr)
	decbot  = max(dec1arr)
	dectop  = min(dec2arr)

	#Crops data so the size of every filter image matches and saves to file 'coadd*.multi.fits'
	#Same for weight file
	for files in coaddfiles:

		newfile = files[:-4]+'multi.fits'
		fitsfile = pf.open(files)
		fitsheader = fitsfile[0].header
		data = fitsfile[0].data
	
		w = wcs.WCS(fitsheader)
		[[x1,y1],[x2,y2]] = w.wcs_world2pix([[raleft, decbot], [raright, dectop]], 0)
		
		pplib.hextractlite(newfile, data, fitsheader, x1+1, x2, y1+1, y2)
		
		wnewfile = files[:-4]+'multi.weight.fits'
		wfitsfile = pf.open(files[:-4]+'weight.fits')
		wfitsheader = wfitsfile[0].header
		wdata = wfitsfile[0].data
	
		w = wcs.WCS(wfitsheader)
		[[x1,y1],[x2,y2]] = w.wcs_world2pix([[raleft, decbot], [raright, dectop]], 0)
		
		pplib.hextractlite(wnewfile, wdata, wfitsheader, x1+1, x2, y1+1, y2)

	#Crops the multicolor (data file and weight) fits files to match the same coordinates
	mixfile = 'multicolor.fits'
	mixfitsfile = pf.open(mixfile)
	mixfitsheader = mixfitsfile[0].header
	mixdata = mixfitsfile[0].data

	mixw = wcs.WCS(mixfitsheader)
	[[mx1,my1],[mx2,my2]] = mixw.wcs_world2pix([[raleft, decbot], [raright, dectop]], 0)
	pplib.hextractlite(mixfile, mixdata, mixfitsheader, mx1+1, mx2, my1+1, my2)

	wmixfile = mixfile[:-4]+'weight.fits'
	wmixfitsfile = pf.open(wmixfile)
	wmixfitsheader = wmixfitsfile[0].header
	wmixdata = wmixfitsfile[0].data

	wmixw = wcs.WCS(wmixfitsheader)
	[[wmx1,wmy1],[wmx2,wmy2]] = mixw.wcs_world2pix([[raleft, decbot], [raright, dectop]], 0)
	pplib.hextractlite(wmixfile, wmixdata, wmixfitsheader, wmx1+1, wmx2, wmy1+1, wmy2)	

	#Find directory where this python code is located
	propath = os.path.dirname(os.path.realpath(__file__))
	
	#Make sure configuration file is in current working directory, if not copy it from
	#location where this code is stored
	if not os.path.exists('ratir_weighted.sex'): 
		os.system('cp '+propath+'/defaults/ratir_weighted.sex .')
	
	if not os.path.exists('weightedtemp.param'): 
		os.system('cp '+propath+'/defaults/weightedtemp.param .')
		
	if not os.path.exists('ratir.conv'): 
		os.system('cp '+propath+'/defaults/ratir.conv .')
				
	if not os.path.exists('ratir_nir.nnw'): 
		os.system('cp '+propath+'/defaults/ratir_nir.nnw .')

	coaddfiles = pplib.choosefiles(prefchar+'*_?.fits')
	#Uses sextractor to find the magnitude and location of sources for each file
	#Saves this information into 'fluxes_*.txt' files
	for files in coaddfiles:

		hdr = pf.getheader(files)

		try:
			filter = hdr['FILTER']
			abszpt = hdr['ABSZPT']
			abszprms = hdr['ABSZPRMS']
			pixscale = hdr['PIXSCALE']
		except KeyError as error:
			print error
			continue

		#Finds filter name and makes sure it is capitalized correctly		
		filter = files.split('_')[-1].split('.')[0]
	
		if filter.lower() in ('j','h','k'):
			filter = filter.upper()
		else:
			filter = filter.lower()
		
		compfile = files

		#Use individual image, not multi image now for individual detections
		#Call to sextractor in double image mode (image1 used for detection of sources, image2 only for measurements - must be same size) 
		os.system('sex -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE '+compfile[:-4]+'weight.fits' + \
			' -c ratir_weighted.sex -SEEING_FWHM 1.5 -PIXEL_SCALE '+str(pixscale)+' -DETECT_THRESH 1.5 -ANALYSIS_THRESH 1.5 -PHOT_APERTURES ' + \
			str(hdr['SPIX'])+ ' -MAG_ZEROPOINT ' + str(hdr['ABSZPT'])+' ' + compfile)
		os.system('mv -f temp.cat fluxes_'+filter+'.txt')
		
		#Columns unpacked for fluxes*.txt are: (x,y,ra,dec,mag,magerr,flux,fluxerr,e,fwhm,flags)
		sexout = np.loadtxt('fluxes_'+filter+'.txt', unpack=True)
		
		magcol = 4
		magerrcol = 5
		
		tout = np.transpose(sexout[0:8,:]) #Only include through fluxerr
		for i in np.arange(len(tout[:,magerrcol])):
			tout[i,magerrcol] = max(tout[i,magerrcol], 0.01)
		
		tout[:,magerrcol] = np.sqrt(tout[:,magerrcol]**2 + abszprms**2)

		tsorted =  tout[np.argsort(tout[:,magcol])]
			
		#Creates Absolute Magnitude file with coordinates
		amfile = 'finalphot'+filter+'.am'		
		np.savetxt(amfile, tsorted, fmt='%15.6f', 
			header='X\t Y\t RA\t DEC\t CAL_MAG\t CAL_MAG_ERR\t CAL_FLUX\t CAL_FLUX_ERR\t')
