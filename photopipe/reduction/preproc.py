"""
Purpose:    this is a collection of preprocessing functions for general use.
"""
import os
from matplotlib.patches import Rectangle
import shutil
from glob import glob
import datetime
import pickle
import sys

# installed modules
import pyfits as pf
import matplotlib.pylab as pl
import numpy as np
from scipy.ndimage.interpolation import zoom

# custom modules/functions
from photopipe.reduction.dependencies.zscale import zscale
from photopipe.reduction.dependencies import astro_functs as af
from photopipe.instruments.specific_instruments import instrument_dict

# Preprocessing constants
FITS_IN_KEY = lambda n: 'IMCMB{:03}'.format(int(n)) # function to make FITS keywords to store file names of combined frames

def choose_calib(instrument, ftype, workdir='.', cams=[0,1,2,3], auto=False, reject_sat=True, amin=0.2, amax=0.8, save_select=True, figsize=(8,5), noplot=False):
    """
    NAME:
        choose_calib
    PURPOSE:
        Either auto-select or display calibration images for user verification
    INPUT:
        instrument  - instrument name defined in instrument_dict (ex. 'ratir')
        ftype       - type of calibration frames (ex. 'flat', 'bias', 'dark')
        workdir     - directory where function executed
        cams        - camera numbers (default is all)
        auto        - automated frame selection. If 'bias', will select all, if 'flat' 
                      will select non-saturated frames with sufficient counts
        reject_sat  - reject frames with saturated pixels
        amin        - minimum fraction of saturation value for median (automated)
        amax        - maximum fraction of saturation value for median (automated)
        save_select - save dictionary of selected frames to python pickle file
    EXAMPLE:
        file_dict = choose_calib('ratir', ftype = bias, dark or flat name, 
            workdir = 'path/to/data/', cams = [#,#,...])
        *** call mkmaster using dictionary or pickle file ***
    """

    instrum = instrument_dict[instrument]

    if auto and (ftype is instrum.flatname):
        temp = raw_input(af.bcolors.WARNING+"Warning: automated selection of flats is not recommended! Continue? (y/n): "+af.bcolors.ENDC)
        if (temp.lower() != 'y') and (temp.lower() != 'yes'):
            af.print_bold("Exiting...")
            return

    # check for non-list camera argument
    if type(cams) is not list:
        cams = [cams] # convert to list

    # check for non-integer camera designators
    if type(cams[0]) is not int:
        af.print_err("Error: cameras must be specified by an integer. Exiting...")
        return
    
    if not noplot: pl.ion() # pylab in interactive mode

    # move to working directory
    start_dir = os.getcwd()
    os.chdir(workdir)
    wordir = os.getcwd()
    
    d = os.getcwd().split('/')[-1] # name of current directory
    if not auto:
        af.print_head("\nDisplaying {} frames in {} for selection:".format(ftype, d))
    else:
        af.print_head("\nAutomatically selecting {} frames in {}:".format(ftype, d))
    
    # dictionary to store selected fits files by camera or filter
    fits_list_dict = {}

    fits_check = glob('????????T??????C??.fits')

    if len(fits_check) == 0:
        files = instrum.change_file_names(glob(instrum.original_file_format()))
        
        fits_check_2 = glob('????????T??????C??.fits')
    
        if len(fits_check_2) == 0:
            af.print_err("Error: no files with correct format after file name changes")
            return

    # open figure for images if not auto
    if not auto and not noplot:
        fig = pl.figure(figsize=figsize)

    # work on FITs files for specified cameras
    for cam_i in cams:

        # print current camera number
        af.print_under("\n{:^50}".format('CAMERA {}'.format(cam_i)))

        # check for valid calibration request
        if instrum.has_cam_bias(cam_i) == False and ftype is instrum.biasname:
            af.print_warn("Warning: Camera C{} does not have {} frames.  "+
                "Skipping...".format(cam_i, instrum.biasname))
            continue
        if instrum.has_cam_dark(cam_i) == False and ftype is instrum.darkname:
            af.print_warn("Warning: Camera C{} does not have {} frames.  "+
                "Skipping...".format(cam_i, instrum.darkname))
            continue
          
        # find raw files of selected type for this camera
        fits_list = glob('????????T??????C{}{}.fits'.format(cam_i,
            instrum.ftype_post[ftype]))
        
        if len(fits_list) == 0:
            af.print_warn("Warning: no fits files found.  Skipping camera {}.".format(cam_i))
            continue

        # look at FITs data sequentially
        for fits_fn in fits_list:
            
            fits_id = fits_fn.split('.')[0] # fits file name with extention removed
            print '{}'.format(fits_fn)

            # open data
            hdulist = pf.open(fits_fn, mode='update')
            im = hdulist[0].data
            h  = hdulist[0].header

            sat_pt = instrum.get_cam_sat(h, cam_i)
                
            if reject_sat:
                if np.any(im == sat_pt):
                    af.print_warn("Warning: saturated pixels in frame.  Skipping frame {}.".format(fits_fn))
                    continue

            # check if instrument camera is split, if it is make sure correct specs being used
            if instrum.is_cam_split(cam_i) == True:
                print '\t* Split filter used'
            else:
                if (instrum.get_filter(h, 'C{}'.format(cam_i)) not in instrum.possible_filters()) and (ftype is instrum.flatname):
                    af.print_warn("Warning: invalid filter detected.  Skipping {} band.".format(instrum.get_filter(h, 'C{}'.format(cam_i))))
                    continue
            
                if ftype is instrum.flatname:
                    print '\t* Filter used: {}'.format(instrum.get_filter(h, 'C{}'.format(cam_i)))
            
            h = instrum.change_header_keywords(h, 'C{}'.format(cam_i)) 
            
            # return quick image summary    
            [im1, m, s, sfrac] = image_summary(im, sat_pt, cam_i, instrum, split=instrum.is_cam_split(cam_i))
            
            if instrum.is_cam_split(cam_i) == True:
                [m1,m2] = m; [s1,s2] = s; [im1,im2] = im1; [sfrac1, sfrac2] = sfrac  
            
            # if auto select then find values with correct ranges
            if auto:

                # all bias and dark frames are selected
                if ftype in [instrum.biasname, instrum.darkname]:
                    addtodict(dict=fits_list_dict, key='C{}'.format(cam_i), value='{}{}'.format(workdir, fits_fn))

                # flats are selected based on median value
                elif ftype is instrum.flatname:

                    vmin = amin * sat_pt; vmax = amax * sat_pt
                    
                    if instrum.is_cam_split(cam_i) == True:
                                                
                        # check whether median values are in specified range
                        # bottom side
                        if m1 > vmin and m1 < vmax:
                            print '\t* Filter used: {}'.format(instrum.get_filter(h,'C{}a'.format(cam_i)))
                            af.print_blue("\t* Bottom side selected.")
                            imfits_1 = savefile(fits_id, im1, instrum.get_filter(h,'C{}a'.format(cam_i)), h)
                            addtodict(dict=fits_list_dict, 
                                key=instrum.get_filter(h,'C{}a'.format(cam_i)), value='{}{}'.format(workdir, imfits_1))

                        else:
                            if m1 < vmin:
                                af.print_warn("\t* Bottom side rejected:\tUNDEREXPOSED.")
                            else:
                                af.print_warn("\t* Bottom side rejected:\tSATURATED.")

                        # top side
                        if m2 > vmin and m2 < vmax:
                            print '\t* Filter used: {}'.format(instrum.get_filter(h,'C{}b'.format(cam_i)))
                            af.print_blue("\t* Top side selected.")
                            imfits_2 = savefile(fits_id, im2, instrum.get_filter(h,'C{}b'.format(cam_i)), h)
                            addtodict(dict=fits_list_dict, 
                                key=instrum.get_filter(h,'C{}b'.format(cam_i)), value='{}{}'.format(workdir, imfits_2))
                        else:
                            if m2 < vmin:
                                af.print_warn("\t* Top side rejected:\tUNDEREXPOSED.")
                            else:
                                af.print_warn("\t* Top side rejected:\tSATURATED.")
                    
                    #Not split frame            
                    else:
                        
                        # check whether median value is in specified range
                        if m > vmin and m < vmax:
                            af.print_blue("\t* Frame selected.")
                            addtodict(dict=fits_list_dict, 
                                key=instrum.get_filter(h,'C{}'.format(cam_i)), value='{}{}'.format(workdir, fits_fn))

                        else:
                            if m < vmin:
                                af.print_warn("\t* Frame rejected:\tUNDEREXPOSED.")
                            else:
                                af.print_warn("\t* Frame rejected:\tSATURATED.")

            # display image and prompt user
            else:
                
                if instrum.is_cam_split(cam_i):
                
                    if (sfrac1 < amin) or (sfrac1 > amax) or (sfrac2 < amin) or (sfrac2 > amax):
                        af.print_warn("Warning: median value outside specified range of {:.0%} - {:.0%} of saturation value in frame.  Skipping frame {}.".format(amin, amax, fits_fn))
                        continue
                        
                    if not noplot:
                        # show top frame
                        ax1 = fig.add_subplot(221)
                        plot_params_calib(ax1, im1, m1, s1, sat_pt, hist=False)
                    
                        # show pixel distribution
                        axhist = fig.add_subplot(222)
                        plot_params_calib(axhist, im1, m1, s1, sat_pt, hist=True)
                    
                        # show bottom frame
                        ax2 = fig.add_subplot(223)
                        plot_params_calib(ax2, im2, m2, s2, sat_pt, hist=False)
                    
                        # show pixel distribution
                        axhist = fig.add_subplot(224)
                        plot_params_calib(axhist, im2, m2, s2, sat_pt, hist=True)
                        fig.subplots_adjust(wspace=0.1, hspace=0.45)
                    
                else:
                    if (sfrac < amin) or (sfrac > amax):
                        af.print_warn("Warning: median value outside specified range of {:.0%} - {:.0%} of saturation value in frame.  Skipping frame {}.".format(amin, amax, fits_fn))
                        continue
                    
                    if not noplot:
                        # show frame
                        ax = fig.add_subplot(121)
                        plot_params_calib(ax, im1, m, s, sat_pt, hist=False)
                    
                        # show pixel distribution
                        axhist = fig.add_subplot(122)
                        plot_params_calib(axhist, im1, m, s, sat_pt, hist=True)
                    
                if not noplot: fig.canvas.draw()
                        
                # query user until valid response is provided
                valid_entry = False
                while not valid_entry:

                    user = raw_input("\nType Y for YES, N for NO, Q for QUIT: ")
                            
                    if user.lower() == 'y':
                                
                        if instrum.is_cam_split(cam_i) == True:
                            imfits_1 = savefile(fits_id, im1,
                                instrum.get_filter(h,'C{}a'.format(cam_i)), h)
                            addtodict(dict=fits_list_dict, 
                                key=instrum.get_filter(h,'C{}a'.format(cam_i)),
                                value='{}{}'.format(workdir, imfits_1))
                                
                            imfits_2 = savefile(fits_id, im2,
                                instrum.get_filter(h,'C{}b'.format(cam_i)), h)
                            addtodict(dict=fits_list_dict, 
                                key=instrum.get_filter(h,'C{}b'.format(cam_i)),
                                value='{}{}'.format(workdir, imfits_2))

                        else:
                            if ftype is instrum.flatname:
                                fl_key = instrum.get_filter(h,'C{}'.format(cam_i))
                            else:
                                fl_key = 'C{}'.format(cam_i)
                            addtodict(dict=fits_list_dict, key=fl_key, value='{}{}'.format(workdir, fits_fn))                            
                                
                        valid_entry = True
                            
                    elif user.lower() == 'q': # exit function
                        af.print_bold("Exiting...")
                        os.chdir(start_dir) # move back to starting directory
                        pl.close('all') # close image to free memory
                        return
                            
                    elif user.lower() == 'n': # 'N' selected, skip
                        valid_entry = True
                            
                    else: # invalid case
                        af.print_warn("'{}' is not a valid entry.".format(user))

            if not auto: fig.clear() # clear image
            hdulist.close() # close FITs file

    if auto:
        if not noplot: 
            af.print_head("\nDisplaying automatically selected {} frames:".format(ftype))
            af.show_list(fits_list_dict)
    else:
        if not noplot: pl.close('all') # close image to free memory

    if save_select:
        dt = datetime.datetime.now()
        fnout = '{}_'.format(ftype)+dt.isoformat().split('.')[0].replace('-','').replace(':','')+'.p' # python pickle extension
        af.print_head("\nSaving selection dictionary to {}".format(fnout))
        pickle.dump( fits_list_dict, open( fnout, 'wb' ) ) # save dictionary to pickle

    os.chdir(start_dir) # move back to starting directory

    return fits_list_dict
    
def image_summary(im, sat_pt, cam_i, instrum, split=False):
    """
    NAME:
        image_summary
    PURPOSE:
        Calculate median, robust scatter, and fraction of saturation point for slice.  If
        split array then will output information for both sides of array
    INPUTS:
        im      - data 
        sat_pt  - saturation point
        cam_i   - camera that is being used
        instrum - module that contains instrument specific information about cameras and split
        split   - boolean that tells if camera has split filters
    """
    
    if split == True:
        im1 = im[instrum.slice('C{}a'.format(cam_i))]
        m1  = np.median(im1)
        s1  = af.robust_sigma(im1)
        sfrac1 = float(m1)/sat_pt
        im2 = im[instrum.slice('C{}b'.format(cam_i))]
        m2  = np.median(im2)
        s2  = af.robust_sigma(im2)
        sfrac2 = float(m2)/sat_pt
        print '\t* Median of left side is {} counts ({:.0%} of saturation level).'.format(m1, sfrac1)
        print '\t* Median of right side is {} counts ({:.0%} of saturation level).'.format(m2, sfrac2)
        
        return [[im1,im2], [m1,m2], [s1,s2],[sfrac1,sfrac2]]
        
    else:
        im1 = im[instrum.slice('C{}'.format(cam_i))]
        m  = np.median(im1)
        s  = af.robust_sigma(im1)
        sfrac = float(m)/sat_pt
        print '\t* Median is {} counts ({:.0%} of saturation level).'.format(m, sfrac)

        return [im1,m,s,sfrac]
        
def addtodict(dict=None, key=None, value=None):
    """
    NAME:
        addtodict
    PURPOSE:
        Adds (key, value) to dictionary by initializing or appending
    INPUTS:
        dict  - dictionary to add to
        key   - key to add to dictionary
        value - value of key to add to dictionary
    EXAMPLE:
        addtodict(dict=dictionary, key='test', value='testing')
    """
    
    try:
        dict[key].append(value)
    except:
        dict[key] = [value]
        
def savefile(file, im, filter, h):
    """
    NAME:
        savefile
    PURPOSE:
        Adds filter keyword and saves file (usually for split filter cameras)
    INPUT:
        file   - filename without extension
        im     - data
        filter - filter name to add
        h      - header
    EXAMPLE:
        savefile('test', data, 'Z', header)
    """
    newfile = '{}_{}.fits'.format(file, filter)
    h['FILTER'] = filter
    if os.path.exists(newfile): os.remove(newfile) # delete old copy
    pf.writeto(newfile, im, header=h, clobber=True) # save object frame
    return newfile
    
def plot_params_calib(ax, im, m, s, sat_pt, hist=False):
    """
    NAME:
        plot_params_calib
    PURPOSE:
        Plots calibration files (image and histogram) for user selection
    INPUTS:
        ax     - plot reference that we will be using to plot images
        im     - data to plot
        m      - median
        s      - robust scatter
        sat_pt - saturation point
    EXAMPLE:
        plot_params_calib(ax, im, m, s, sat_pt, hist=False)
    NOTE:
        Will not plot unless you have a show() or something to display    
    """
    
    z1, z2 = af.zscale(im)
    if z2 <= z1:
        z1 = m1 - s1; z2 = m1 + s1
          
    if hist:
        ax.hist(im.flatten(), bins=50, normed=True, log=True, range=(0, sat_pt))
        ax.set_xlim((0, sat_pt))
        ax.set_xticks([0, 0.5*sat_pt, sat_pt])
        ax.set_xticklabels(['0%', '50%', '100%'])
        ax.set_yticks([])
        ax.grid()
        ax.set_title("Pixel distribution")
        ax.set_ylabel("Log Scale")
    else:
        ax.imshow(im, vmin=z1, vmax=z2, origin='lower', cmap=pl.cm.gray, interpolation='none')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(r"Median = {}, $\sigma$ = {:.1f}".format(int(m), s))
        ax.set_xlabel(r"Median is {:.0%} of saturation level.".format(float(m)/sat_pt))    
        
def plot_params_science(ax, disp_im, filter, h, central=False, window_zoom=4):

    """
    NAME:
        plot_params_science
    PURPOSE:
        Plots science files (image and histogram) for user selection
    INPUTS:
        ax          - plot reference that we will be using to plot images
        disp_im     - data to plot
        filter      - filter for labeling plot
        h           - header
        central     - boolean to center image to zoomed window
        window_zoom - zoom level for closer look
    NOTE:
        Will not plot unless you have a show() or something to display    
    """

    z1, z2 = af.zscale(disp_im)
    ax.set_xticks([])
    ax.set_yticks([])
    xm, ym = np.array(disp_im.shape, dtype=float)/2
    xr = xm/float(window_zoom); yr = ym/float(window_zoom)
    
    if central:    
        ax.imshow(disp_im[xm-xr:xm+xr,ym-yr:ym+yr], vmin=z1, vmax=z2, origin='lower', 
            cmap=pl.cm.gray, interpolation='none')
        ax.contour(disp_im[xm-xr:xm+xr,ym-yr:ym+yr], levels=[z2], 
            origin='lower', colors='r')
        ax.set_title("Zoomed region")
        
    else:
        ax.imshow(disp_im, vmin=z1, vmax=z2, origin='lower', 
            cmap=pl.cm.gray, interpolation='none')
        ax.contour(disp_im, levels=[z2], origin='lower', colors='r')
        ax.add_patch(Rectangle((ym-yr, xm-xr), 2*yr, 2*xr, ec='b', fc='none', lw=2))
        ax.set_title(r"{} band".format(filter))
        
def choose_science(instrument, workdir='.', targetdir='.', cams=[0,1,2,3], auto=False, save_select=True, 
                    figsize=(10,10), window_zoom=4, calibrate=False, noplot=False):

    """
    PURPOSE:
        Display science images for verification by user
    INPUT:
        instrument  - instrument name defined in instrument_dict (ex. 'ratir')
        workdir     - directory where function is to be executed
        targetdir   - directory where selected frames and lists are output
        cams        - camera numbers to process data, all by default
        auto        - select all science frames
        save_select - save dictionary of selected frames to python pickle file
        figsize     - dimensions of figure used to display frames for selection
        window_zoom - zoom level for closer look
    EXAMPLE:
        file_dict = choose_science('ratir', workdir = 'path/to/data/', 
            targetdir = 'path/to/processeddata/',cams = [#,#,...], calibrate=True)
    """

    instrum = instrument_dict[instrument]

    # check for non-list camera argument
    if type(cams) is not list:
        cams = [cams] # convert to list

    # check for non-integer camera designators
    if type(cams[0]) is not int:
        af.print_err("Error: cameras must be specified by an integer. Exiting...")
        return
    
    pl.ion() # pylab in interactive mode

    # move to working directory
    start_dir = os.getcwd()
    os.chdir(workdir)
    
    d = os.getcwd().split('/')[-1] # name of current directory
    if not auto:
        af.print_head("\nDisplaying science frames in {} for selection:".format(d))
    else:
        af.print_head("\nSelecting all science frames in {}:".format(d))

    # dictionary to store selected fits files by camera or filter
    fits_list_dict = {}
    
    # remove tailing / from target directory name if present
    if targetdir[-1] == '/':
        targetdir = targetdir[:-1]

    # make target directory if it does not exist
    if not os.path.exists(targetdir):
        af.print_blue("Creating target directory: {}".format(targetdir))
        os.makedirs(targetdir)
    # warn user if previous files may be overwritten
    else:
        af.print_warn("Warning: Target directory exists. Existing files will be overwritten.")
        resp = raw_input("Proceed? (y/n): ")
        if resp.lower() != 'y':
            af.print_bold("Exiting...")
            os.chdir(start_dir) # move back to starting directory
            return
        else:
            shutil.rmtree(targetdir)
            os.makedirs(targetdir)

    fits_check = glob('????????T??????C??.fits')

    if len(fits_check) == 0:
        files = instrum.change_file_names(glob(instrum.original_file_format()))
        
        fits_check_2 = glob('????????T??????C??.fits')
    
        if len(fits_check_2) == 0:
            af.print_err("Error: no files with correct format after file name changes")
            return

    # open figure for images if not auto
    if not auto and not noplot:
        fig = pl.figure(figsize=figsize)

    # work on FITs files for specified cameras
    for cam_i in cams:
        
        # print current camera number
        af.print_under("\n{:^50}".format('CAMERA {}'.format(cam_i)))
        
        if calibrate:
            # get master dark and bias frames for current camera if required
            if instrum.has_cam_bias(cam_i) == True:
                # if instrum.has_cam_bias(cam_i) == True:
                mbias_fn = '{}_C{}.fits'.format(instrum.biasname, cam_i)
                
                if not os.path.exists(mbias_fn):
                    af.print_err('Error: {} not found.  Move master bias file to working directory to proceed.'.format(mbias_fn))
                    continue
                else:
                    mbias_data = pf.getdata(mbias_fn)
                    
            if instrum.has_cam_dark(cam_i) == True:
                mdark_fn = '{}_C{}.fits'.format(instrum.darkname, cam_i)

                if not os.path.exists(mdark_fn):
                    af.print_err('Error: {} not found.  Move master dark file to working directory to proceed.'.format(mdark_fn))
                    continue
                else:
                    mdark_data = pf.getdata(mdark_fn)        
                    
        # find raw files of selected type for this camera
        fits_list = glob('????????T??????C{}{}.fits'.format(cam_i, instrum.ftype_post[instrum.objname]))
        if len(fits_list) == 0:
            af.print_warn("Warning: no fits files found.  Skipping camera {}.".format(cam_i))
            continue

        # look at FITs data sequentially
        for fits_fn in fits_list:

            fits_id = fits_fn.split('.')[0] # fits file name with extension removed
            print '{}'.format(fits_fn)
                
            # open data
            hdulist = pf.open(fits_fn)
            im = hdulist[0].data
            h = hdulist[0].header

            if calibrate:
                # get master flat frame for current filter
                if instrum.is_cam_split(cam_i) == True:
                    mflat_fn1 = '{}_{}.fits'.format(instrum.flatname, instrum.get_filter(h,'C{}a'.format(cam_i)))
                    if not os.path.exists(mflat_fn1):
                        af.print_err('Error: {} not found.  Move master flat file to working directory to proceed.'.format(mflat_fn1))
                        continue
                    else:
                        mflat_data1 = pf.getdata(mflat_fn1)
                    mflat_fn2 = '{}_{}.fits'.format(instrum.flatname, instrum.get_filter(h,'C{}b'.format(cam_i)))
                    if not os.path.exists(mflat_fn2):
                        af.print_err('Error: {} not found.  Move master flat file to working directory to proceed.'.format(mflat_fn2))
                        continue
                    else:
                        mflat_data2 = pf.getdata(mflat_fn2)
                
                else:
                    mflat_fn = '{}_{}.fits'.format(instrum.flatname, instrum.get_filter(h,'C{}'.format(cam_i)))
                    if not os.path.exists(mflat_fn):
                        af.print_err('Error: {} not found.  Move master flat file to working directory to proceed.'.format(mflat_fn))
                        continue
                    else:
                        mflat_data = pf.getdata(mflat_fn)
                                
            # get image statistics
            if instrum.is_cam_split(cam_i) == True:            
                im1 = im[instrum.slice('C{}a'.format(cam_i))]
                im2 = im[instrum.slice('C{}b'.format(cam_i))]
            else:
                im1 = im[instrum.slice('C{}'.format(cam_i))]
                    
            # display image and prompt user
            if not auto:
            
                if instrum.is_cam_split(cam_i) == True:
                    
                    disp_im1 = np.copy(im1)
                    disp_im2 = np.copy(im2)
                    
                    if calibrate:
                        if instrum.has_cam_bias(cam_i):
                            disp_im1 -= mbias_data
                            disp_im2 -= mbias_data
                        
                        if instrum.has_cam_dark(cam_i):
                            disp_im1 -= mbias_data*instrum.get_exptime(h)
                            disp_im2 -= mbias_data*instrum.get_exptime(h)
                        
                        disp_im1 = np.divide(disp_im1, mflat_data1)
                        disp_im2 = np.divide(disp_im2, mflat_data2)
                    
                    if not noplot:
                        # display top
                        ax1 = fig.add_subplot(221)
                        plot_params_science(ax1, disp_im1, instrum.get_filter(h,'C{}a'.format(cam_i)),
                            h, central=False)

                        # and central subregion
                        ax1s = fig.add_subplot(222)
                        plot_params_science(ax1s, disp_im1, instrum.get_filter(h,'C{}a'.format(cam_i)),
                            h, central=True)

                        # display bottom
                        ax2 = fig.add_subplot(223)
                        plot_params_science(ax2, disp_im2, instrum.get_filter(h,'C{}b'.format(cam_i)),
                            h, central=False)

                        # and central subregion
                        ax2s = fig.add_subplot(224)
                        plot_params_science(ax2s, disp_im2, instrum.get_filter(h, 'C{}b'.format(cam_i)),
                            h, central=True)
                
                else:
                    disp_im1 = np.copy(im1)
                    
                    if calibrate:
                        if instrum.has_cam_bias(cam_i):
                            disp_im1 -= mbias_data
                        
                        if instrum.has_cam_dark(cam_i):
                            disp_im1 -= mdark_data*instrum.get_exptime(h)
                        
                        disp_im1 = np.divide(disp_im1, mflat_data)
                    
                    if not noplot:
                        ax = fig.add_subplot(121)
                        plot_params_science(ax, disp_im1, instrum.get_filter(h, 'C{}'.format(cam_i)),
                            h, central=False)

                        # and central subregion
                        axs = fig.add_subplot(122)
                        plot_params_science(axs, disp_im1, instrum.get_filter(h, 'C{}'.format(cam_i)),
                            h, central=True)                    
                
                if not noplot:
                    fig.set_tight_layout(True)
                    fig.canvas.draw()
            
            if instrum.is_cam_split(cam_i) == True:
                if instrum.get_centered_filter(h, cam_i).count(instrum.get_filter(h, 'C{}a'.format(cam_i))) != 0:
                    print "\t* The target is focused on the {} filter.".format(instrum.get_filter(h, 'C{}a'.format(cam_i)))
                elif instrum.get_centered_filter(h, cam_i).count(instrum.get_filter(h, 'C{}b'.format(cam_i))) != 0:
                    print "\t* The target is focused on the {} filter.".format(instrum.get_filter(h, 'C{}b'.format(cam_i)))
                else:
                    af.print_warn("\t* Warning: The target is NOT focused on an split filter. The target is focused on the {} filter.".format(instrum.get_centered_filter(h, cam_i)))
            else:
                # print filter name
                print '\t* Filter used: {}'.format(instrum.get_filter(h, 'C{}'.format(cam_i)))

            # query user until valid response is provided
            valid_entry = False
            while not valid_entry:
                    
                # either select all if auto, or have user select
                if auto:
                    user = 'y'
                else:
                    user = raw_input("\nType Y for YES, N for NO, Q for QUIT: ")
                
                if user.lower() == 'y' and instrum.is_cam_split(cam_i): 
                    if instrum.get_centered_filter(h, cam_i).count(instrum.get_filter(h, 'C{}b'.format(cam_i))) != 0:
                        direction = 't'
                    elif instrum.get_centered_filter(h, cam_i).count(instrum.get_filter(h, 'C{}a'.format(cam_i))) != 0:
                        direction = 'b'
                    else:
                        af.print_warn("\t* Warning: Skipping frame not centered on split filter.")
                        user = 'n'
                        direction = ''
               
                if user.lower() == 'y':

                    h_c = h.copy()
                    
                    if instrum.is_cam_split(cam_i) == True:
                        
                        if direction.lower() == 'b':
                            f_img_post = 'a'
                            f_sky_post = 'b'
                        else:
                            f_img_post = 'b'
                            f_sky_post = 'a'
                        
                        imfits = '{}/{}_{}_{}.fits'.format(targetdir, fits_id, instrum.objname, 
                            instrum.get_filter(h_c,'C{}{}'.format(cam_i, f_img_post)))
                        im_img = im[instrum.slice('C{}{}'.format(cam_i, f_img_post))]
                        hnew = instrum.change_header_keywords(h_c, 'C{}{}'.format(cam_i, f_img_post))
                        pf.writeto(imfits, im_img, header=hnew, clobber=True) # save object frame
                        if fits_list_dict.has_key(instrum.get_filter(h_c,'C{}{}'.format(cam_i, f_img_post))):
                            fits_list_dict[instrum.get_filter(h_c,'C{}{}'.format(cam_i, f_img_post))].append(imfits)
                        else:
                            fits_list_dict[instrum.get_filter(h_c,'C{}{}'.format(cam_i, f_img_post))] = [imfits]
                        
                        # filter side with sky, now saved as object, but different list to keep track
                        skyfits = '{}/{}_{}_{}.fits'.format(targetdir, fits_id, instrum.objname, 
                            instrum.get_filter(h_c,'C{}{}'.format(cam_i, f_sky_post)))
                        im_sky = im[instrum.slice('C{}{}'.format(cam_i, f_sky_post))]
                        hnew = instrum.change_header_keywords(h_c, 'C{}{}'.format(cam_i, f_sky_post))
                        pf.writeto(skyfits, im_sky, header=hnew, clobber=True) # save sky frame
                        if fits_list_dict.has_key(instrum.get_filter(h_c,'C{}{}'.format(cam_i, f_sky_post))):
                            fits_list_dict[instrum.get_filter(h_c,'C{}{}'.format(cam_i, f_sky_post))].append(skyfits)
                        else:
                            fits_list_dict[instrum.get_filter(h_c,'C{}{}'.format(cam_i, f_sky_post))] = [skyfits]

                        valid_entry = True

                    else:
                        
                        imfits = '{}/{}_{}_{}.fits'.format(targetdir, fits_id, instrum.objname, cam_i)
                        im_img = im[instrum.slice('C{}'.format(cam_i))]
                        hnew = instrum.change_header_keywords(h_c, 'C{}'.format(cam_i))
                        pf.writeto(imfits, im_img, header=hnew, clobber=True)
                        if fits_list_dict.has_key(instrum.get_filter(h_c,'C{}'.format(cam_i))):
                            fits_list_dict[instrum.get_filter(h_c,'C{}'.format(cam_i))].append(imfits)
                        else:
                            fits_list_dict[instrum.get_filter(h_c,'C{}'.format(cam_i))] = [imfits]

                        valid_entry = True                      
                    
                elif user.lower() == 'q': # exit function
                    af.print_bold("Exiting...")
                    os.chdir(start_dir) # move back to starting directory
                    pl.close('all') # close image to free memory
                    return
                
                elif user.lower() == 'n': # 'N' selected, skip
                    valid_entry = True
                        
                else: # invalid case
                    af.print_warn("'{}' is not a valid entry.".format(user))

            if not auto:
                fig.clear() # clear image
            hdulist.close() # close FITs file

    if auto:
        if not noplot:
            af.print_head("\nDisplaying automatically selected science frames:")
            af.show_list(fits_list_dict)
    else:
        if not noplot:
            pl.close('all') # close image to free memory
            
    if save_select:
        dt = datetime.datetime.now()
        fnout = 'object_'+dt.isoformat().split('.')[0].replace('-','').replace(':','')+'.p' # python pickle extension
        af.print_head("\nSaving selection dictionary to {}".format(fnout))
        pickle.dump( fits_list_dict, open( fnout, 'wb' ) ) # save dictionary to pickle

    os.chdir(start_dir) # move back to starting directory

    return fits_list_dict
    

def mkmaster(instrument, fn_dict, mtype, fmin=5):

    """
    PURPOSE:        
        Make master calibration frames (bias, dark, flat)
        * currently no outlier rejection other than median combine
    INPUT:
        instrument  - instrument name defined in instrument_dict (ex. 'ratir')
        fn_dict     - dictionary output by choose_calib() containing organized 
                      fits file names.  can also provide file name of pickled dictionary.
        mtype       - type of master frame. should be either 'flat', 'dark' or 'bias'
        fmin        - minimum number of files needed to make a master frame
    EXAMPLE:
        mkmaster('ratir', fn_dict=output from choose_calib(), mtype = bias, dark or flat name)
    FUTURE IMPROVEMENTS:
        - Better outlier rejection
    """
    instrum = instrument_dict[instrument]

    # check if input is a file name
    if type(fn_dict) is str:
        if fn_dict.split('.')[-1] == 'p':
            af.print_bold("Loading pickled dictionary from file.")
            fn_dict = pickle.load(open(fn_dict, 'rb'))
        else:
            af.print_err("Invalid pickle file extension detected. Exiting...")
            return

    # check for valid mtype
    if mtype not in [instrum.flatname, instrum.biasname, instrum.darkname]:
        af.print_err("Error: valid arguments for mtype are {}, {} and {}. Exiting...".format(instrum.flatname, instrum.biasname, instrum.darkname))
        return

    bands = fn_dict.keys()

    sorttype = 'BAND'

    if mtype in [instrum.biasname, instrum.darkname]:
        sorttype = 'CAMERA' 
        
    d = os.getcwd().split('/')[-1] # name of current directory
    af.print_head("\nMaking master {} frame in {}:".format(mtype, d))

    # work on FITs files for specified photometric bands
    for band in bands:

        # print current band
        af.print_under("\n{:^50}".format('{} {}'.format(band, sorttype)))

        first_file = fn_dict[band][0]
        ind_C = first_file.index('C')
        cam = first_file[ind_C:ind_C + 2]
        cam_i = int(cam[1])
  
        # check if required files are present
        if instrum.has_cam_bias(cam_i):
            mbias_fn = '{}_{}.fits'.format(instrum.biasname, cam)
        else:
            mbias_fn = None
            
        if instrum.has_cam_dark(cam_i):
            mdark_fn = '{}_{}.fits'.format(instrum.darkname, cam)
        else:
            mdark_fn = None
            
        if mtype is not instrum.biasname:
            if mbias_fn is not None:
                if not os.path.exists(mbias_fn):
                    af.print_err('Error: {} not found.  Move master bias file to working directory to proceed.'.format(mbias_fn))
                    continue
            if mdark_fn is not None:
                if (mtype is instrum.flatname) and (not os.path.exists(mdark_fn)):
                    af.print_err('Error: {} not found.  Move master dark file to working directory to proceed.'.format(mdark_fn))
                    continue

        # check dictionary entries
        fns = fn_dict[band]
        if len(fns) < fmin:
            if len(fns) == 0:
                af.print_err('Error: no frames available to make master {} for {} {}.'.format(mtype, band, sorttype.lower()))
                continue
            else:
                temp = raw_input(af.bcolors.WARNING+"Only {} frames available to make master {} for {} {}.  Continue? (y/n): ".format(len(fns), mtype, band, sorttype.lower())+af.bcolors.ENDC)
                if temp.lower() != 'y' and temp.lower() != 'yes':
                    af.print_warn("Skipping {}...".format(band))
                    continue

        # load calibration data
        hdu = pf.PrimaryHDU()
        filter_arr = [] # to check that all frames used the same filter
        exptime_arr = [] # to check that all frames have the same exposure time (where relevant)
        data_arr = []
        i = 0
        for fn in fns:
            print fn
            hdu.header[FITS_IN_KEY(i)] = fn # add flat fn to master flat header
            hdulist = pf.open(fn)
            data_arr.append(hdulist[0].data)
            filter_arr.append(hdulist[0].header['FILTER'])
            exptime_arr.append(hdulist[0].header['EXPTIME'])
            i += 1
        data_arr = np.array(data_arr, dtype=np.float)

        # check that frames match
        for i in range(len(fns)-1):
            if (filter_arr[i+1] != filter_arr[0]) and (mtype is instrum.flatname):
                af.print_err("Error: cannot combine flat frames with different filters. Skipping {} {}...".format(band, sorttype.lower()))
                continue
            if (exptime_arr[i+1] != exptime_arr[0]) and (mtype is instrum.darkname):
                af.print_err("Error: cannot combine dark frames with different exposure times. Skipping {} {}...".format(band, sorttype.lower()))
                continue
        if instrum.flatname:
            hdu.header['FILTER'] = filter_arr[0] # add filter keyword to master frame
        if instrum.darkname:
            hdu.header['EXPTIME'] = exptime_arr[0] # add exposure time keyword to master frame
        
        # add CAMERA header keyword
        hdu.header['CAMERA'] = cam_i  # add camera keyword to master frame

        # crop bias frames
        if mtype is instrum.biasname:
            data_arr = data_arr[(np.s_[:],instrum.slice(cam)[0],instrum.slice(cam)[1])]
        # crop dark frames and perform calculations
        elif mtype is instrum.darkname:
            data_arr = data_arr[(np.s_[:],instrum.slice(cam)[0],instrum.slice(cam)[1])]
            data_arr = (data_arr - pf.getdata(mbias_fn))/hdu.header['EXPTIME'] # calculate dark current
        # crop flat frames and perform calculations
        elif mtype is instrum.flatname:
            if mbias_fn is not None:
                mbd = pf.getdata(mbias_fn)
                mbd = mbd
            if mdark_fn is not None:
                mdd = pf.getdata(mdark_fn)
                mdd = mdd
            
            if instrum.is_cam_split(cam_i):
                pass # split data is already cropped
            else:  
                data_arr = data_arr[(np.s_[:],instrum.slice(cam)[0],instrum.slice(cam)[1])]

            for i in range(len(exptime_arr)):
                if mbias_fn is not None:
                    data_arr[i] -= mbd
                if mdark_fn is not None:
                    data_arr[i] -= mdd*exptime_arr[i]
                data_arr[i] /= np.median(data_arr[i])

        # make master frame
        master = af.imcombine(data_arr, type='median').astype(np.float)

        # add master to hdu
        if mtype is instrum.flatname:
            hdu.data = master/np.median(master) # normalize master flat
        else:
            hdu.data = master

        # save master to fits
        hdulist = pf.HDUList([hdu])
        hdulist.writeto('{}_{}.fits'.format(mtype, band), clobber=True)