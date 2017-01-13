"""
NAME:
	plotratir
PURPOSE:
	Read in prefchar*(FILTER).multi.ref.fits and finalphot(FILTER).am.  Find number
	of unique stars and saves each filter's magnitude and error to the same star. 
	Saves this to finalmags.txt.  Creates color plot with all of the images overlaid
	(red = J/H, green = z/y, blue = r/i) and creates plot of each filter field with
	green circles around each object.  Saves as coadd*(FILTER).png
	Calls printhtml which create HTML to view all information
	Can now take in different filters
OUTPUTS:
	finalmags.txt      - text file with all the magnitudes for each filter of each source
	color.png          - overlay of all filters (red = J/H, green = z/y, blue = r/i)
	coadd*(FILTER).png - images of each filter field with green circles over source
	ratir.html         - html showing all information
	
Slowest part is saving image using pl.savefig, can't find a workaround
Need to modify colors based on different types of filters being provided

Translated from plotratir.pro by John Capone (jicapone@astro.umd.edu).
Modified on 7/31/2014 by Vicki Toy (vtoy@astro.umd.edu)
"""

import numpy as np
import astropy.io.fits as pf
from scipy.ndimage.interpolation import zoom
from scipy.misc import bytescale
import pylab as pl
import scipy as sp
from astropy import wcs

import photprocesslibrary as pplib
import printhtml
import os
import sys

def plotphotom(prefchar='coadd'):

	#Initialize arrays
    filters = ['r','i','z','y','J','H']
    arr_size = 10000
    
    #Creates a list of names: ['ra', 'dec', '(FILTER)mag', '(FILTER)magerr'] including
    #all filters in filters list
    #Also initializes filter index dictionary for use later
    names = ['ra', 'dec']
    ifiltdict = {}
    for filter in filters:
    	names.extend([filter+'mag', filter+'magerr'])
    	ifiltdict[filter] = -1

    #Create dictionary with keys from names list and all set to np.zeros(arr_size)
    #easily changed for other filters
    plotdict = {}
    for name in names:
    	plotdict[name] = np.zeros(arr_size)

    # retrieve detection files
    zffiles = pplib.choosefiles( prefchar+'*_?.ref.multi.fits' )

    #Save filter and data from each file into arrays and find overlapping stars
    #between filter files using final photometry file (comparing distances from RA and DEC)
    cfilterarr = []
    imgarr = []
    harr   = []

    for i in range(len(zffiles)):
    
    	#Find filter of each file and make sure that has the right capitalization
        cfilter = zffiles[i].split('_')[-1].split('.')[0]
        if cfilter == 'Z' or cfilter == 'Y':
            cfilter = cfilter.lower()
        cfilterarr.append(cfilter)
        
        #Save data scaled by scale factor to imgarr
        ifile = zffiles[i]
        hdulist = pf.open(ifile)
        h = hdulist[0].header
        img = hdulist[0].data
        im_size = np.shape(img)
        scalefactor = 1.
        img = zoom( img, 1./scalefactor, order=0 )
        imgarr.append(img)
        harr.append(h)
        
        #Read in finalphot[FILTER].am which has the instrument corrected photometry
        pfile = 'finalphot' + cfilter + '.am'
        try:
            x, y, ra, dec, mag, magerr = np.loadtxt(pfile, unpack=True)
        except IOError as error:
            print error
            continue

        clen = len(mag)
        
        #For first file initialize variables for following files, use temporary variables
        #for comparison
        if i == 0:
            plotdict['ra'][0:clen]  = ra
            plotdict['dec'][0:clen] = dec
            
            plotdict[cfilter+'mag'][0:clen]    = mag
            plotdict[cfilter+'magerr'][0:clen] = magerr
            nstars = clen
        else:
            compra     = ra
            compdec    = dec
            compmag    = mag
            compmagerr = magerr
            
            #For each source in file find any sources that are within 1 arcsecond
            #if these exist then store information in same index but different filter's magnitude
            #array.  If these don't exist, put on the end of filter's (and position) arrays
            #to signify a new source

            for j in range(len(compra)):
                smatch = pplib.nearest( compra[j]*np.cos(compdec[j]*np.pi/180.), compdec[j], 
                						plotdict['ra']*np.cos(plotdict['dec']*np.pi/180.), plotdict['dec'], maxdist=1./3600. )
                
                if any(smatch):
                	plotdict[cfilter+'mag'][smatch]    = compmag[j]
                	plotdict[cfilter+'magerr'][smatch] = compmagerr[j]
                else:
                    plotdict['ra'][nstars]  = compra[j]
                    plotdict['dec'][nstars] = compdec[j]
                    plotdict[cfilter+'mag'][nstars]    = compmag[j]
                    plotdict[cfilter+'magerr'][nstars] = compmagerr[j]
                    nstars += 1
    
    imgarr = np.array(imgarr)
	
	#Save stars to finalmags.txt with correct format and removes zeros
    store = np.zeros(nstars)
    for name in names:
    	store = np.vstack( (store,plotdict[name][:nstars]) )
    
    store = store[1:, :] #Removes 0's from initialization

    #Finds sources that are cut off on at least one filter using weight maps
    #and checking if 25% of the number of pixels in a 10 pixel radius circle
    #are 0 in the weight map of each filter (i.e. cut off)
    sra = store[0]
    sdec = store[1]   
    removesource = []
    
    for s in np.arange(len(sra)):
    	for file in zffiles:
    		nzero = 0
    		wmultifile = file[:-4]+'weight.fits'
    		hlist = pf.open(wmultifile)
    		wdata = hlist[0].data
    		wh    = hlist[0].header
    		w     = wcs.WCS(wh)
    		spix  = w.wcs_world2pix(sra[s],sdec[s], 1)
    		
    		wcir  = pplib.circle(spix[0],spix[1],10)
    		ncir  = len(wcir)
    		for coord in wcir:
    			if (wh['NAXIS1'] > coord[0] > 0) & (wh['NAXIS2'] > coord[1] > 0):
    				if wdata[int(coord[1])][int(coord[0])] == 0: nzero += 1
    			else:
    				ncir -= 1
    		if nzero >= 0.25*ncir: 
    			removesource.append(s)
    			break    

    store = np.delete(store, removesource, axis=1)
    np.savetxt('finalmags.txt', store.T, fmt='%12.6f')
    
    #Find the index of the file that corresponds to each filter and save 
    #to ifiltdict (initialized to -1)
    for item in ifiltdict:
    	try:
    		ifiltdict[item] = cfilterarr.index(item)
    	except ValueError:
    		pass

	#Determines colors based on which filters are present.  
	#Red = J/H, green = z/y, blue = r/i
	#If neither filter present, set to 0, if one present, use imgarr of data from that filter
	#if both present use half from imgarr of data from each filter	
	def fcolor(filt1, filt2, ifiltdict, imgarr):
	
		if filt1 and filt2 in ifiltdict:
			if ifiltdict[filt1] >= 0 and ifiltdict[filt2] >= 0:
				x = imgarr[ifiltdict[filt1],:,:] * 0.5 + imgarr[ifiltdict[filt2],:,:] * 0.5
			if ifiltdict[filt1] >= 0 and ifiltdict[filt2] < 0:
				x = imgarr[ifiltdict[filt1],:,:]
			if ifiltdict[filt2] >= 0 and ifiltdict[filt1] < 0:
				x = imgarr[ifiltdict[filt2],:,:]
			if ifiltdict[filt2] < 0 and ifiltdict[filt1] < 0:
				x = 0
		else:
			print 'Valid filters were not supplied, set color to 0'
			x = 0
			
		return x 
        
    red   = fcolor('J', 'H', ifiltdict, imgarr)    
    green = fcolor('z', 'y', ifiltdict, imgarr) 
    blue  = fcolor('r', 'i', ifiltdict, imgarr)

	#Determine image size base on if color filter exists (priority: red, green, blue in that order)
    if np.size(red) > 1:
        im_size = np.shape(red)
    elif np.size(green) > 1:
        im_size = np.shape(green)
    elif np.size(blue) > 1:
        im_size = np.shape(blue)

    def bytearr( x, y, z ):
        return np.zeros((x,y,z)).astype(np.uint8)
    
    #Create color of image and save to color.png    
    color = bytearr( im_size[0], im_size[1], 3 )     
        
    #Changes color into bytescale range and saves to color array
    if np.size(blue) > 1:
        blue  = bytescale(blue,  0, 8, 250)
        color[:,:,2] = blue * 0.5
    if np.size(green) > 1:
        green = bytescale(green, 0, 8, 250)
        color[:,:,1] = green * 0.5
    if np.size(red) > 1:
        red   = bytescale(red,   0, 8, 250)
        color[:,:,0] = red * 0.5
        
    color = color[::-1,:]
    fig = pl.figure('color image')
    pl.axis('off')
    pl.imshow( color, interpolation='None', origin='lower' )
    sp.misc.imsave( 'color.png', color )
    
    #Find aperture size to make circles around sources that match sextractor aperture
    pline=''
    sfile = open('ratir_weighted.sex', 'r')
    for line in sfile:
    	if 'PHOT_APERTURES' in line: pline=line
    sfile.close()
    
    bpline = pline.split()	
    
    if len(bpline) != 0:
    	aper = int(bpline[1])/2.0 #radius
    else:
    	aper = 10
    
    objra  = store[0]
    objdec = store[1]
    
    #Plot each image with circles on star identification
    for i in range(len(zffiles)):
    	ifile   = zffiles[i]
    	img     = imgarr[i]
    	h       = harr[i]
    	cfilter = cfilterarr[i]
    	
    	scale   = bytescale(img, 0, 10, 255)
    	dpi     = 72. # px per inch
    	figsize = (np.array(img.shape)/dpi)[::-1]
    	fig     = pl.figure(i)
    	
    	pl.imshow( scale, interpolation='None', cmap=pl.cm.gray, origin='lower' )
    	xlims   = pl.xlim()
    	ylims   = pl.ylim()
    	
    	# Parse the WCS keywords in the primary HDU    	
    	w       = wcs.WCS(h)
    	world   = np.transpose([objra, objdec])
    	pixcrd  = w.wcs_world2pix(world, 1)
    	
    	fs = 20
    	fw = 'normal'
    	lw = 2
    	
    	#For each star create a circle and plot in green
    	#If pixel coordinates of star (from WCS conversion of RA and DEC) and within the 
    	#x and y limits, then put text on right side, otherwise put on left
    	for j in range(len(objra)):
    		ctemp = pplib.circle( pixcrd[j][0]/scalefactor, pixcrd[j][1]/scalefactor, aper ).T
    		pl.plot( ctemp[0], ctemp[1], c='#00ff00', lw=lw )
    		if pixcrd[j][0]/scalefactor+40 < xlims[1] and pixcrd[j][1]/scalefactor+20 < ylims[1]:
    			pl.text( pixcrd[j][0]/scalefactor+15, pixcrd[j][1]/scalefactor, `j`, color='#00ff00', fontsize=fs, fontweight=fw )
    		else:
    			pl.text( pixcrd[j][0]/scalefactor-45, pixcrd[j][1]/scalefactor-20, `j`, color='#00ff00', fontsize=fs, fontweight=fw )
    			
    	#Label plot and remove axes, save to filename+.png
    	pl.text( 0.2*xlims[1], 0.9*ylims[1], cfilter+'-Band', color='r', fontsize=fs, fontweight=fw )
        a = pl.gca()
        a.set_frame_on(False)
        a.set_xticks([]); a.set_yticks([])
        pl.axis('off')
        pl.xlim(xlims)
        pl.ylim(ylims)
        fig.set_size_inches(figsize[0], figsize[1])
        ofile = ifile.split('.')[0] + '.png'
        pl.savefig( ofile, bbox_inches='tight', pad_inches=0, transparent=True, dpi=dpi )
   		
    #Create HTML to do quick look at data
    printhtml.printhtml(filters, names)