import os
import glob
import pyfits as pf
import numpy as np
import datetime
import astrometrystats as astst
import cosmics
import matplotlib.pyplot as plt
import timeit
import scipy

def pipeprepare(filename, outname=None, biasfile=None, darkfile=None, verbose=1):

    """
    NAME:
        pipeprepare
    PURPOSE:
        Adds additional header keywords needed later on in the pipeline and removes 
        unnecessary header keywords by looking through a list of mandatory keywords.  
        Also runs bias and dark subtraction for filters with an existing master bias/dark 
        (CCDs).  The prepared images are written to disk with outname
    INPUTS:
        filename - name of FITS file, or array of filenames, or file w/list of filenames 
    OPTIONAL KEYWORDS:
        outname  - specify output file to write to disk
        biasfile - full name (including path) of master bias
        darkfile - full name (including path) of master dark
        verbose  - print out comments        
    EXAMPLE:
        pipeprepare(filename, outname=outname, biasfile=biasfile, darkfile=darkfile, verbose=1)
    DEPENDENCIES:
        autoproc_depend.pipeprepare()
    """
    
    # ------ Process input filenames(s) ------
    
    # Check for empty filename
    if len(filename) == 0:
        print 'No filename specified'
        return

    # If string, check if a file of items or wildcards
    # otherwise store all files
    if isinstance(filename,str):
        fileext = os.path.splitext(filename)[1][1:]
        
        files = [filename]
        
        if fileext in ['cat', 'lis', 'list', 'txt']:
            f = open(filename,'r')
            files = f.read().splitlines() 
            f.close()            
                      
        if '?' in filename or '*' in filename:
            files = glob.glob(filename)
            if len(files) == 0: 
                print 'Cannot find any files matching ', filename
                return
    else:
        files = filename


    # ------ Read data and process header information ------        
    for file in files: 
        f = pf.open(file)
        head = f[0].header
        data = f[0].data
        f.close()
        
        # If these keys exist keep, otherwise delete all extraneous keywords
        mandatorykey = ['SIMPLE','BITPIX','NAXIS','NAXIS1','NAXIS2',
 					    'HISTORY','DATE-OBS','EXPTIME','INSTRUME',
 					    'LATITUDE','LONGITUD','BINNING','BINY','BINX',
 					    'CAMERA','TARGNAME','UTC','OBJECT', 'OBJNAME','AIRMASS', 					    
 					    'GAIN','SATURATE','PIXSCALE','FILTER','WAVELENG',
 					    'CD1_1','CD1_2','CD2_1','CD2_2',
 					    'CRPIX1','CRPIX2','CRVAL1','CRVAL2','CTYPE1','CTYPE2', 
 					    'PV1_1','PV2_1','PV1_17','PV2_17','PV1_19','PV2_19','PV1_21','PV2_21',
 					    'PV1_31','PV2_31','PV1_33','PV2_33','PV1_35','PV2_35','PV1_37','PV2_37']
 					    
 	    # Finds list of unnecessary keywords, then deletes extraneous
        newhead = head
        for oldkey in head.keys():
            if oldkey not in mandatorykey:
                try:
                    del newhead[oldkey]
                except KeyError:
                    pass
        
        # If biasfile keyword set subtract master bias from current file with given master bias file
        # If they are not the same size, quit program without saving with preparation prefix (will not move
        # on in following processing steps)
        if biasfile != None:
            bias = pf.getdata(biasfile)
            
            if np.shape(data) != np.shape(bias):
                
                print file + ' could not be bias subtracted because it is not the same' +\
                             ' size as the master bias, remove file to avoid confusion'
                return
            
            if verbose > 0: print '    bias subtracting'
            
            newdata = data - bias
            
            # If darkfile keyword set subtract master dark from current file with given master dark file
            # If they are not the same size, quit program without saving with preparation prefix (will not move
            # on in following processing steps)
            if darkfile != None:
                dark = pf.getdata(darkfile) * newhead['EXPTIME']
                
                if np.shape(data) != np.shape(dark):
                    print ' '
                    print file + ' could not be dark subtracted because it is not the same' +\
                                 ' size as the master dark, remove file to avoid confusion'
                    return  
                          
                if verbose > 0: print '    dark subtracting'
                
                newdata = newdata - dark
            else:
                print file, 'could not be dark subtracted because the master dark file was not provided'
        else:
            newdata = data
        
        # Write changes to disk
        pf.writeto(outname, newdata, newhead, clobber=True)
        
        if verbose > 0: print file, '-> ', outname
        
def flatpipeproc(filename, flatname, flatminval=0, flatmaxval=0):

    """
    NAME:
        flatpipeproc
    PURPOSE:
        Checks if flat is same size as data, then divides for correct filter
    INPUTS:
        filename - name of FITS file, or array of filenames, or file w/list of filenames 
        flatname - name of FITS master flat file
    OPTIONAL KEYWORDS:
        flatminval - if not set to 0 below this value will set to NaNs
        flatmaxval - if not set to 0 above this value will set to NaNs
    EXAMPLE:
        flatpipeproc(filename, flatname, flatminval=0.3)
    """

    # ------ Process input filenames(s) ------
    
    # Check for empty filename
    if len(filename) == 0:
        print 'No filename specified'
        return

    # If string, check if a file of items or wildcards
    # otherwise store all files
    if isinstance(filename,str):
        fileext = os.path.splitext(filename)[1][1:]
        
        files = [filename]
        
        if fileext in ['cat', 'lis', 'list', 'txt']:
            f = open(filename,'r')
            files = f.read().splitlines() 
            f.close()            
                      
        if '?' in filename or '*' in filename:
            files = glob.glob(filename)
            if len(files) == 0: 
                print 'Cannot find any files matching ', filename
                return
    else:
        files = filename
        
    flat = pf.getdata(flatname)
    
    med = np.median(flat)
    if (med < 0.5) or (med > 2.0): print 'Warning: flat is not normalized to one'
    
    for file in files:
        f = pf.open(file)
        data = f[0].data
        head = f[0].header
        f.close()
        
        if np.shape(data) != np.shape(flat):
            print file + ' could not be dark subtracted because it is not the same' +\
                         ' size as the master dark, remove file to avoid confusion'
            return  
        
        # Set values too low/high to NaNs
        if flatminval > 0:
            flat[flat < flatminval] = float('NaN')
        goodsignal = np.where(flat-1.0 < 0.1)
   
        if flatmaxval > 0:
            flat[flat > flatminval] = float('NaN')
        
        # Divides out flattened field and adds keywords to header to show change
        fdata = data / flat         
        
        head['FLATFLD'] = flatname
        skycts = np.median(fdata[goodsignal])        
        head['SKYCTS']  = (skycts, 'Sky counts')
        
        try:
            head['CTRATE'] = (skycts/head['EXPTIME'], 'Sky counts per second')
        except:
            print 'No EXPTIME keyword'
            
        date = datetime.datetime.now().isoformat()
        head.add_history('Processed by flatproc ' + date)
        
        fileroot = os.path.basename(file)
        filedir  = os.path.dirname(file)
        outnameim = filedir + '/f' + fileroot
        
        pf.writeto(outnameim, fdata, head, clobber=True)
        
def skypipecombine(filelist, outfile, filt, pipevar, removeobjects=None, 
    objthresh=6, algorithm='median', trimlo=None, trimhi=None, mincounts=1, 
    maxcounts=55000, satlevel=30000, type=None):
    
    """
    NAME:
        skypipecombine
    PURPOSE:
        Create sigma clipped median sky flat.  Scales each file based on the overall 
        sigma clipped median, then removes objects selected with sextractor (uses flux 
        fraction radius) in each file. Removes saturated pixels.  Calculates sigma clipped 
        median of each pixel and saves anything with non-finite values (saturated or 
        source) to the median of the entire frame.  Save with outfile name.
    INPUTS:
        filelist - files to be processed
        outfile  - name for output fits file
        filt	 - filter of files
        pipevar  - pipeline parameters in dictionary   
    OPTIONAL KEYWORDS:
        removeobjects 	- specifies if you want objects removed
        objthresh		- sets sigma in removeobjects (default is 6)
        algorithm       - algorithm to solve (mean or median, default is median)
        trimlo			- trim off bottom of data in mean algorithm mode (default is 25%)
        trimhi			- rim off top of data in mean algorithm mode (default is 25%)
        mincounts		- sets minimum counts allowed (default is 1)
        maxcounts		- sets maximum counts allowed (default is 55000)
        satlevel		- sets saturation level (default is 30000)
        type			- sets 'SKYTYPE' keyword in header of outfile to this string
    EXAMPLE:
        skypipecombine(filelist, 'sky-filt.fits', removeobjects=True, type='sky')
    DEPENDENCIES:
        medclip, Sextractor
    FUTURE IMPROVEMENTS:
        medclip slow find faster solution
        Need to take saturation level from header?
        Saved header is from middle file, maybe use blank?
    """
    
    # Sets defaults for trimming (25% of list)
    if trimlo != None: (len(filelist)+1)/4
    if trimhi != None: trimlo
    
    # If given list, then grab all filenames, saved to files
    if len(filelist) == 1:
        f = open(filelist,'r')
        files = f.read().splitlines() 
        f.close()
    else:
        files = filelist
    
    nfiles = len(files)
    nmid = len(files)/2
    
    # Read in middle file and initialize arrays
    f = pf.open(files[nmid])
    data_m = f[0].data
    head_m = f[0].header
    f.close()
    
    nx = head_m['NAXIS1']
    ny = head_m['NAXIS2']
    
    data     = np.zeros((nfiles, ny, nx)) + float('NaN')
    skymeds  = []
    usefiles = []
    
    z = 0
    # For each file and make sure size matches middle file, calculate sigma clipped 
    # median (3sig, 6 iter), then if within counts limit save data into 3d data cube 
    # and save clipped median into skymed, and mark file as usable
    # Increment z by one when this is true
    for file in files:
        f = pf.open(file)
        data_i = f[0].data
        head_i = f[0].header
        f.close()        
        
        inx = head_i['NAXIS1']
        iny = head_i['NAXIS2']        

        if (inx != nx) or (iny != ny):
            print 'File ' + file + ' has wrong dimensions ('+str(inx)+ \
                  ' x '+ str(iny)+'; should have '+str(nx)+' x '+str(ny)+')'
        
        # Perform 3 sigma clipped median and save to inmeds
        inmed, instd = medclip(data_i, clipsig=3, maxiter=3)
        
        
        # If median is within limits save data, otherwise exclude files
        if inmed >= mincounts and inmed <= maxcounts:
            if pipevar['verbose'] > 0:
                print file + ' ('+str(inmed) + ' counts/pix)'
            
            skymeds  += [inmed]
            usefiles += [file]
            data[z, :,:] = data_i
            z += 1
        else:
            if inmed < mincounts:
                print file + ' (' + str(inmed) + ' counts/pix) - too few counts; excluding'
            if inmed > maxcounts:
                print file + ' (' + str(inmed) + ' counts/pix) - too many counts; excluding'            
        
    if z < 2:
        print 'ERROR - Not enough counts to make a flat with these data!'
        return
    
    # Median of sigma clipped medians
    medsky = np.median(skymeds)
    
    # Scale each file by median of sigma clipped medians divided by median of data
    # Corrects for each flat's changing sky background
    for f in np.arange(z-1):
        factor = medsky / skymeds[f]
        data[f,:,:] = data[f,:,:]*factor
    
    # Removes extraneous indexes in data for skipped files
    if z != nfiles: data = data[0:z,:,:]
    
    # Removes objects from field by calculating iterative median sigma clipping 
    # (5 sigma, 5 iter) and using the calculated stddev to remove 6sigma (or non-default 
    # object threshold) data from the median along with values above the saturation limit.
    # Find 3sigma clipped median of each pixel from remaining values or mean of middle 50%
    
    if removeobjects != None:
        if pipevar['verbose'] > 0: print '  Identifying objects...'
                
        for f in np.arange(z-1):
            
            indata = data[f,:,:]
            
            # Set sources above objthresh  limit to NaN
            datamed, datastd = medclip(indata, clipsig=5, maxiter=5)
            
            sourcepixels = np.where( abs(indata-datamed) >= objthresh*datastd)
                        
            satpixels = np.where( indata >= satlevel )
            
            if len(sourcepixels[0]) > 0:
                indata[sourcepixels] = float('NaN')

            if len(satpixels[0]) > 0:
                indata[satpixels] = float('NaN')
            
            data[f,:,:] = indata
        
    reflat = np.zeros((ny, nx)) + float('NaN')

    # If algorithm set to median, find 3 sigma clipped median of each pixel 
    # excluding NaN values (which are eventually set to median)
    if algorithm == 'median':
        if pipevar['verbose'] > 0: print '  Median-combining...'
            
        for y in np.arange(ny):
            
            vector = data[:,y,:]
            temp = np.isfinite(vector)
            me, st = medclip2d(vector, clipsig=3, maxiter=5, overaxis=0)
            reflat[y,:] = me

        # Replace bad pixels with median of entire sky
        good = np.isfinite(reflat)
        allmed = np.median(reflat[good])
        bad = ~good # Opposite of boolean array good
        reflat[bad] = allmed
        
    # If algorithm set to mean, takes mean of trimmed sorted values. Default is to 
    # trim 25% off top and bottom, if not enough good data, set trimming to 0
    if algorithm == 'mean':
        
        if pipevar['verbose'] > 0: print '  Combining via trimmed mean...'
            
        for y in np.arange(ny):
            for x in np.arange(nx):
                slice = data[:,y,x]
                good = np.isfinite(slice)
                    
                cslice = slice[good]
                ctgood = len(cslice)
                    
                if ctgood == 0:
                    reflat[y,x] = 1
                    
                itrimlo = trimlo
                itrimhi = trimhi
                    
                while ctgood-itrimlo-itrimhi < 1:
                    itrimlo = max(itrimlo - 1, 0)
                    itrimhi = max(itrimhi - 1, 0)
                        
                cslice = np.sort(cslice)
                cslice = cslice[itrimlo:ctgood-itrimhi]
                reflat[y,x] = np.mean(cslice)
        
    # Adds header information to signify what files we used 
    for f in np.arange(z-1):
        head_m['SKY'+str(f)] = usefiles[f]
            
    if type != None:
        head_m['SKYTYPE'] = type
        
    date = datetime.datetime.now().isoformat()
    head_m.add_history('Processed by skypipecombine ' + date) 
        
    if pipevar['verbose'] > 0: print '  Written to ' + outfile
        
    pf.writeto(outfile, reflat, head_m, clobber=True)       
                    
def skypipeproc(filename, flatname, outfile, flatminval=None, flatmaxval=None):

    """
    NAME:
        skypipeproc
    PURPOSE:
        Subtracts sky flat from data and then subtracts median of that from remaining data. 
    INPUTS:
    	filename - file or list of files to be sky subtracted
    	flatname - sky flat fits file 
    	outfile  - name of output file
    OPTIONAL KEYWORDS:
        flatminval - minimum required value in flat (default for skycts calculation is 0.1)
        flatmaxval - maximum required value in flat
    EXAMPLE:
        skypipeproc(filename, flatname, outfile)        
    """
    
    # ------ Process input filenames(s) ------
    
    # Check for empty filename
    if len(filename) == 0:
        print 'No filename specified'
        return

    # If string, check if a file of items or wildcards
    # otherwise store all files
    if isinstance(filename,str):
        fileext = os.path.splitext(filename)[1][1:]
        
        files = [filename]
        
        if fileext in ['cat', 'lis', 'list', 'txt']:
            f = open(filename,'r')
            files = f.read().splitlines() 
            f.close()            
                      
        if '?' in filename or '*' in filename:
            files = glob.glob(filename)
            if len(files) == 0: 
                print 'Cannot find any files matching ', filename
                return
    else:
        files = filename  
    
    # Open flat    
    flat = pf.getdata(flatname) 
    
    med = np.median(flat)
    
    # For each input file check if same size as flats (required). If there is a minimum 
    # or maximum flat value set, forces values outside of that range to NaN. Use finite 
    # values above 0.1 to determine skycounts, and subtract flat along with median of 
    # flattened data. Saves to new fits file
    for file in files:
        f = pf.open(file)
        data = f[0].data
        head = f[0].header
        f.close()
        
        if np.shape(data) != np.shape(flat):
            print file + ' could not be flat subtracted because it is not the same' +\
                         ' size as the master flat, remove file to avoid confusion'
            return

        if flatmaxval != None:
            w = np.where(flat > flatminval) 
            if len(w[0]) != 0:
                flat[w] = float('NaN')  
                
        if flatminval != None:
            w = np.where(flat < flatminval)  
            if len(w[0]) != 0:
                flat[w] = float('NaN')
            goodsignal = np.where((flat >= flatminval) & (np.isfinite(flat)))
        else:
            goodsignal = np.where((flat >= 0.1) & (np.isfinite(flat)))          
                
        # Scale skyflat, subtract scaled skyflat, and subtract median of subsequent flat 
        # subtracted data. Calculate skycounts from data (above minimum, or 
        # by default above 0.1)
        flattmp = np.median(flat[np.isfinite(flat)])
        imgtmp  = np.median(data[np.isfinite(data)])
        
        scalefr = imgtmp/flattmp
        fdata   = data - scalefr * flat
        
        tmp     = np.median(fdata[np.isfinite(fdata)])
        fdata   = fdata - tmp

        skycts  = np.median(fdata[goodsignal])
        
        # Adds header keywords to denote new median counts and file we used to flatfield
        head['SFLATFLD'] = flatname
        head['SKYCTS']   = (skycts, 'Sky counts')
        
        try:
            head['CTRATE'] = (skycts/head['EXPTIME'], 'Sky counts per second')
        except:
            print 'No EXPTIME keyword'        

        date = datetime.datetime.now().isoformat()
        head.add_history('Processed by skypipeproc ' + date)
        
        pf.writeto(outfile, fdata, head, clobber=True) 

def cosmiczap(filename, outname, sigclip=6.0, maxiter=3, verbose=True):

    """
    NAME:
        cosmiczap
    PURPOSE:
        Removes cosmic rays using Laplacian cosmic ray identification written in cosmics.py 
    INPUTS:
    	filename - file or list of files to be cosmic ray zapped
    	outfile  - name of output file
    OPTIONAL KEYWORDS:
        sigclip  - sigma to clip
        maxiter  - maximum number of times to iterate loop
        verbose  - quiet?
    EXAMPLE:
        cosmiczap(filename, outname)
    DEPENDENCIES:
        cosmic.py (described in http://arxiv.org/pdf/1506.07791v3.pdf)  
    FUTURE IMPROVEMENTS:
        Read readnoise from header?    
    """
    
    data, head = cosmics.fromfits(filename, verbose=False)
    
    gain = head['GAIN']
    c = cosmics.cosmicsimage(data, gain=gain, readnoise=18, sigclip=sigclip,
        sigfrac = 0.5, objlim = 5.0, verbose=False)
    
    tot = c.run(maxiter=maxiter, verbose=False)
    
    head['NPZAP'] = (tot, "Num. of pixels zapped by cosmiczap")
    date = datetime.datetime.now().isoformat()
    head.add_history('Processed by cosmiczap ' + date)    
    
    if verbose: print '  Zapped %d total affected pixels (%.3f%% of total)' \
                      %(tot,tot*100.0/np.size(data))
    
    cosmics.tofits(outname, c.cleanarray, head, verbose=False)

def astrometry(atfimages, scamprun=1, pipevar=None):

    """
    NAME:
        astrometry
    PURPOSE:
        Run sextractor and scamp to refine astrometric solution
    INPUTS:
        atfimages - list of files to run through scamp
        scamprun  - the first run does a LOOSE run with distortion degree 1, any
                    other run will look for high distortion parameters, if it
                    finds it will use distortion degree 7, otherwise 3 (will also cut out
                    FLXSCALE on runs after 1)
    EXAMPLE:
        astrometry(atfimages, scamprun=2, pipevar=pipevar)
    FUTURE IMPROVEMENTS:
        Better difference between scamp runs.
    """
        
    acatlist = ''
    scat = {'sdss': 'SDSS-R7', 'tmpsc': '2MASS', 'tmc': '2MASS', 'ub2': 'USNO-B1'} 

    for cfile in atfimages:
        head = pf.getheader(cfile)
        pixscale  = head['PIXSCALE']
        sourcecat = head['ASTR_CAT']
                
        trunfile = os.path.splitext(cfile)[0]
                
        if pipevar['verbose'] > 0:
            sexcom = pipevar['sexcommand'] + ' -CATALOG_NAME ' + trunfile + \
                    '.cat -CATALOG_TYPE FITS_LDAC -FILTER_NAME astrom.conv ' + \
                    '-PARAMETERS_NAME astrom.param -DETECT_THRESH 2.0 ' + \
                    '-ANALYSIS_THRESH 2.0 -PIXEL_SCALE ' + str(pixscale) + \
                    ' ' + cfile
            print sexcom
        else:
            sexcom = pipevar['sexcommand'] + ' -CATALOG_NAME ' + trunfile + \
                    '.cat -CATALOG_TYPE FITS_LDAC -FILTER_NAME astrom.conv ' + \
                    '-PARAMETERS_NAME astrom.param -DETECT_THRESH 2.0 ' + \
                    '-ANALYSIS_THRESH 2.0 -VERBOSE_TYPE QUIET -PIXEL_SCALE ' + \
                    str(pixscale) + ' ' + cfile
                
        os.system(sexcom)
                
        if head['ASTR_NUM'] > 0: acatlist += ' ' + trunfile + '.cat'
                
    if sourcecat in scat:
        cat_u = scat[sourcecat]
    else:
        cat_u = 'NONE'
        print 'No valid catalogs available for SCAMP, check that ' +\
                'vlt_autoastrometry.py ran correctly'
        return
    
    if scamprun == 1:
        loose = ' -MOSAIC_TYPE LOOSE'
        distdeg = 1
    else:
        loose = ' '
        try:
            distort = head['PV1_37']
            distdeg = 7
        except:
            distdeg = 3 
                
    if pipevar['verbose'] > 0:
        scampcmd = "scamp -POSITION_MAXERR 0.2 -DISTORT_DEGREES " + str(distdeg)+\
                    loose + " -ASTREF_CATALOG " + cat_u + \
                    " -SOLVE_PHOTOM N -SN_THRESHOLDS 3.0,10.0 " + \
                    "-CHECKPLOT_DEV NULL -WRITE_XML N -VERBOSE_TYPE FULL " +\
                    acatlist
        print scampcmd
    else:
        scampcmd = "scamp -POSITION_MAXERR 0.2 -DISTORT_DEGREES " + str(distdeg)+\
                    loose + " -ASTREF_CATALOG " + cat_u + \
                    " -SOLVE_PHOTOM N -SN_THRESHOLDS 3.0,10.0 " + \
                    "-CHECKPLOT_DEV NULL -WRITE_XML N -VERBOSE_TYPE QUIET " +\
                    acatlist                                    
                
    os.system(scampcmd)
    os.system('rm ' + acatlist)
            
    # Adds header information to file and delete extra files
    for cfile in atfimages:
        trunfile = os.path.splitext(cfile)[0]
                
        if pipevar['verbose'] > 0:
            os.system('missfits -WRITE_XML N ' + cfile)
        else:
            os.system('missfits -WRITE_XML N -VERBOSE_TYPE QUIET' + cfile)
                    
        os.system('rm ' + trunfile + '.head ' + cfile + '.back')

        if scamprun != 1:
            him  = pf.getheader(cfile)
            data = pf.getdata(cfile)
            del him['FLXSCALE']
            pf.update(cfile, data, him)

def findsexobj(file, sigma, pipevar, masksfx=None, zeropt=25.0, maptype='MAP_WEIGHT',
               wtimage=None, fwhm=1.5, pix=0.3787, aperture=5.0, elong_cut=1.5, 
               quiet=0):
    """
    NAME:
        findsexobj
    PURPOSE:
        Finds sextractor objects with optional inputs. Estimates seeing from stars found. 
    INPUTS:
    	file    - fits file to run sextractor on
    	sigma   - detection threshold and analysis threshold for sextractor
    	pipevar - pipeline parameters (typically set in autopipedefaults or autoproc)
    	
    OPTIONAL KEYWORDS:
    	masksfx   - text string identifier for sextractor CHECKIMAGE_NAME
    	zeropt    - input value for sextractor MAG_ZEROPOINT
    	wtimage   - file for input for sextractor WEIGHT_IMAGE
    	fwhm      - input value for sextractor SEEING_FWHM
    	pix       - input value for sextractor PIXEL_SCALE
    	aperture  - input value for sextractor PHOT_APERTURES
    	elong_cut - cutoff limit for FWHM calculation of elongation to eliminate non-stars
    	quiet     - no output from sextractor if set
    EXAMPLE:
        findsexobj(file, 3.0, pipevar, aperture=20.0)
    DEPENDENCIES:
        sextractor
    FUTURE IMPROVEMENTS:
        More keywords to sextractor?
    """
    
    # Move necessary sextractor configuration files if they are not in current directory
    if not os.path.isfile('coadd.param'): 
        os.system('cp ' + pipevar['defaultspath'] + '/coadd.param .')
    if not os.path.isfile('coadd.conv'): 
        os.system('cp ' + pipevar['defaultspath'] + '/coadd.conv .') 
    if not os.path.isfile('coadd.config'): 
        os.system('cp ' + pipevar['defaultspath'] + '/coadd.config .') 
    if not os.path.isfile('default.nnw'): 
        os.system('cp ' + pipevar['defaultspath'] + '/default.nnw .')  

    if quiet > 0: 
        verbosetype = 'QUIET'
    else:
        verbosetype = 'NORMAL'
        
    # Run sextractor with given input parameters. Saves temp.cat as 
    # starfile, saves starmask, and calculates seeing from starlike objects. Saves 
    # necessary parameters to header
    if file == '': return
    
    if not os.path.isfile(file): return
    starfile = file + '.stars'
        
    trunfile = os.path.splitext(file)[0]
        
    sexcmd = pipevar['sexcommand'] + ' -c coadd.config -DETECT_THRESH ' +\
             str(sigma) + ' -ANALYSIS_THRESH ' + str(sigma) + ' -PHOT_APERTURES ' +\
             str(aperture) + ' -MAG_ZEROPOINT ' + str(zeropt) + ' -PIXEL_SCALE ' +\
             str(pix) + ' -SEEING_FWHM ' + str(fwhm) + ' -VERBOSE_TYPE ' +verbosetype
    
    if masksfx != None:
        mskimg = trunfile + '_' + masksfx + '.fits'
        sexcmd += ' -CHECKIMAGE_TYPE OBJECTS' + ' -CHECKIMAGE_NAME ' + mskimg
        
    if wtimage != None:
        sexcmd += ' -WEIGHT_TYPE '+maptype+' -WEIGHT_IMAGE ' + wtimage + ' '
        
    sexcmd += ' ' + file
    if quiet == 0: print sexcmd
    os.system(sexcmd)
        
    if quiet == 0: print 'mv -f test.cat ' + starfile
    os.system('mv -f test.cat ' + starfile)
    
    num = 0    
    # Calculates seeing with starlike objects
    if os.path.isfile(starfile):
        vars   = np.loadtxt(starfile, unpack=True)
        num    = vars[0,:]
        flag   = vars[5,:]
        elon   = vars[8,:]
        fwhmim = vars[9,:]
        keep = (flag ==0) & (elon < elong_cut) & (fwhmim > fwhm) & (fwhmim < 20.0)

        if sum(keep) <= 1: 
            seepix = None
        else:
            seepix = np.median(fwhmim[keep])        
    else:
        print 'Failed to find Sextractor output file!'
        seepix = None
	 
    head = pf.getheader(file)
    
    if masksfx != None:
        head['MASKNAME'] = (mskimg, "Object mask image from Sextractor")
    
    head['STARFILE'] = (starfile, "Objects file from Sextractor" )
    head['ZEROPT']   = (zeropt, "Photometric zero-point used for Sextractor")
    if seepix != None:
        head['SEEPIX']   = (seepix, "Estimated seeing from Sextractor objects (pix)")
    head['NSTARS']   = (len(num), "Estimated number of objects from Sextractor")
    
    data = pf.getdata(file)
    pf.update(file, data, head)
    
    # Removes config files after done
    os.system('rm -f coadd.param')
    os.system('rm -f coadd.conv')
    os.system('rm -f coadd.config')
    os.system('rm -f default.nnw')
    
def calc_zpt(catmag, obsmag, wts, sigma=3.0, plotter=None):

    """
    NAME:
        calc_zpt
    PURPOSE:
        Find zeropoint using robust scatter
    INPUTS:
        catmag  - 2d array with catalog magnitudes catmag[nobs,nstar]
        obsmag  - 2d array with observed magnitudes obsmag[nobs,nstar]
        wts     - 2d array with weights wts[nobs,nstar]
    OPTIONAL KEYWORDS:
        sigma   - sigma value for how far values can be from robust scatter value
        plotter - filename to save zeropoint plot
    OUTPUTS:
        z2     - zeropoint correction
        scats  - robust scatter of each observation
        rmss   - standard deviation (without bad weight points) of each observation
    EXAMPLE:
        zpt,scats,rmss = calc_zpt(catmag,obsmag,wts, sigma=3.0)     
    """
  
    # Find difference between catalog and observed magnitudes
    diff = catmag - obsmag

    print diff
    # Find number of observations and stars	
    sz = np.shape(obsmag)
    nobs   = sz[0]
    nstars = sz[1]
    
    # For each observation (i.e. frame) find the weighted difference and store zeropoint
    # and new magnitude with zeropoint correction
    z = []
    modmag = np.copy(obsmag)
    for i in np.arange(nobs):
        indz = sum(diff[i,:]*wts[i,:])/sum(wts[i,:])
        z += [indz]
        modmag[i,:] = obsmag[i, :] + indz
    
    # Find difference of catalog and zeropoint corrected values. Remove any values with 
    # weights set to 0 or lower.  Calculate robust scatter on these values.  If difference 
    # with these weights is not within sigma*robust scatter then set weight to 0
    adiff1 = catmag - modmag
    scats, rmss = robust_scat(adiff1, wts, nobs, nstars, sigma)

    z2 = []
    # Recalculate zeropoint using corrected weights (difference still same)        
    modmag2 = np.copy(obsmag)
    for i in np.arange(nobs):
        indz = sum(diff[i,:]*wts[i,:])/sum(wts[i,:])
        z2 += [indz]
        modmag2[i,:] = obsmag[i, :] + indz
    
    adiff2 = catmag - modmag2
    # Recalculate robust scatter and rms scatter value on twice zeropoint corrected mags
    scats, rmss = robust_scat(adiff2, wts, nobs, nstars, sigma)   
    
    if plotter != None:

        keep = np.where(wts != 0)
        plt.plot(catmag[keep], adiff2[keep], '*')
        plt.errorbar(catmag[keep], adiff2[keep], yerr = 1.0/np.sqrt(wts[keep]), fmt='.')
        
        plt.ylabel('Difference between Catalog and Observed')
        plt.xlabel('Catalog magnitude')
        plt.savefig(plotter)
        plt.clf()

    return z2, scats, rmss

def robust_scat(diff, wts, nobs, nstars, sigma):

    """
    NAME:
        robust_scat
    PURPOSE:
        Calculate robust scatter and set the weight of those above this limit to 0
    INPUTS:
        diff   - values to calculate robust scatter over
        wts    - weighting (0 is bad)
        nobs   - number of observations to iterate over
        nstars - number of stars to iterate over
        sigma  - sigma*robust scatter that is acceptable
    OUTPUTS:
        scats  - robust scatter of each observation
        rmss   - standard deviation (without bad weight points) of each observation
    EXAMPLE:
        robust_scat(diff, wts, 1, 12, 3)
    """
    
    scats = np.zeros(nobs)
    rmss  = np.zeros(nobs)
    for i in np.arange(nobs):
        goodwts = np.where( wts[i,:] > 0 )
        if len(goodwts[0]) == 0: continue
        gooddiff = diff[i,goodwts]
        
        # Median absolute deviation
        scat = 1.48 * np.median(abs(gooddiff-np.median(gooddiff)))
        for j in np.arange(nstars):
            if abs(diff[i,j] - np.median(gooddiff)) > (sigma*scat):
                wts[i,j] = 0
        scats[i] = scat
        rmss[i]  = np.std(gooddiff)
    return scats, rmss

def medclip(indata, clipsig=3.0, maxiter=5, verbose=0):

    """
    NAME:
        medclip
    PURPOSE:
        Median iterative sigma-clipping
    INPUTS:
        indata - array to be clipped
    OPTIONAL KEYWORDS:
        clipsig - sigma to clip around
        maxiter - maximum number of times to clip
        verbose - allow to print messages
    EXAMPLE:
        med, sigma = medclip(indata, sigma=5.0)
    """
    
    # Flatten array
    skpix = indata.reshape( indata.size, )
 
    keep = np.isfinite(skpix)
    skpix = skpix[keep]
    ct = indata.size
    iter = 0
    numrej = len(skpix)
    ndata  = len(skpix)
    
    while (iter < maxiter) and (numrej > min(ndata*0.01, 50)):
        lastct = ct
        medval = np.median(skpix)
        sig = np.std(skpix)
        wsm = np.where( abs(skpix-medval) < clipsig*sig )
        ct = len(wsm[0])
        if ct > 0:
            skpix = skpix[wsm]
 
        numrej = abs(ct - lastct)
        if ct <=2: return 'Too few remaining'
        iter += 1
 
    med   = np.median( skpix )
    sigma = np.std( skpix )
 
    if verbose:
        print '%.1f-sigma clipped median' % (clipsig)
        print 'Mean computed in %i iterations' % (iter)
        print 'Mean = %.6f, sigma = %.6f' % (med, sigma)
 
    return med, sigma

def medclip2d(indata, clipsig=3.0, maxiter=5, verbose=0, overaxis=0):

    """
    NAME:
        medclip2d
    PURPOSE:
        Median iterative sigma-clipping over 2d array
    INPUTS:
        indata - array to be clipped
    OPTIONAL KEYWORDS:
        clipsig - sigma to clip around
        maxiter - maximum number of times to clip
        verbose - allow to print messages
        overaxis - axis that we want to take median over
    EXAMPLE:
        med, sigma = medclip2d(indata, sigma=5.0, overaxis=0)
    """
        
    # Flatten array
    skpix = np.ma.masked_array(indata)
    skpix = np.ma.masked_invalid(skpix)
    iter = 0
    
    while (iter < maxiter):
        medval = np.ma.median(skpix, axis=overaxis)
        sig = np.ma.std(skpix, axis=overaxis)
        
        if overaxis == 0:
            mask = ( abs(skpix-medval) < clipsig*sig )
        else:
            mask = ( abs(skpix.T-medval) < clipsig*sig )
        if (mask == skpix.mask).all:
            break
        skpix.mask = mask

        #if ct <=2: return 'Too few remaining'
        iter += 1
 
    med   = np.ma.median( skpix, axis=overaxis )
    sigma = np.ma.std( skpix, axis=overaxis )
 
    if verbose:
        print '%.1f-sigma clipped median' % (clipsig)
        print 'Mean computed in %i iterations' % (iter)
        print 'Mean = %.6f, sigma = %.6f' % (med, sigma)
 
    return med, sigma  
    
def identify_matches( queried_stars, found_stars, match_radius=3. ):
    '''
    Use a kd-tree (3d) to match two lists of stars, using full spherical coordinate distances.
    
    queried_stars, found_stars: numpy arrays of [ [ra,dec],[ra,dec], ... ] (all in decimal degrees)
    match_radius: max distance (in arcseconds) allowed to identify a match between two stars.
    
    Returns two arrays corresponding to queried stars:
    indices - an array of the indices (in found_stars) of the best match. Invalid (negative) index if no matches found.
    distances - an array of the distances to the closest match. NaN if no match found.
    '''
    # make sure inputs are arrays
    queried_stars = np.array(queried_stars)
    found_stars = np.array(found_stars)
    
    ra1, dec1 = queried_stars[:,0], queried_stars[:,1]
    ra2, dec2 = found_stars[:,0], found_stars[:,1]
    dist = 2.778e-4*match_radius # convert arcseconds into degrees 
    
    cosd = lambda x : np.cos(np.deg2rad(x))
    sind = lambda x : np.sin(np.deg2rad(x))
    mindist = 2 * sind(dist/2) 
    getxyz = lambda r, d: [cosd(r)*cosd(d), sind(r)*cosd(d), sind(d)]
    xyz1 = np.array(getxyz(ra1, dec1))
    xyz2 = np.array(getxyz(ra2, dec2))
    
    tree2 = scipy.spatial.KDTree(xyz2.transpose())
    ret = tree2.query(xyz1.transpose(), 1, 0, 2, mindist)
    dist, ind = ret
    dist = np.rad2deg(2*np.arcsin(dist/2))
    ind[ np.isnan(dist) ] = -9999
    
    return ind, dist