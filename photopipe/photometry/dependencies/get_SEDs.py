'''
Library to produce a catalog of fully-populated SEDs and calculating
zeropoints for any arbitrary location on the sky, using a combination 
of online catalogs (USNOB1, 2MASS, SDSS) and synthetic photometry.

To Do:
- add proper distutils setup

'''


############################################
# IMPORTS
############################################
import numpy as np
from subprocess import Popen, PIPE
from scipy.optimize import fmin_bfgs
from scipy.interpolate import interp1d
import scipy.spatial.kdtree
from os.path import isfile, dirname, join
from threading import Thread
from time import time, strftime
from urllib2 import urlopen
import pickle
import sys
import warnings

import multiprocessing as mp
N_CORES = mp.cpu_count()  # use all the cpus you have

# use the __file__ variable to point to the static files
#  Note: __file__ points to the location of this file,
#  so the static files below MUST be in the same folder.
try:
    MODELS = np.load( open( join(dirname(__file__), 'all_models.npy'),'r') )
    # rezero so that K=0 for all models
    for row in MODELS[1:]:
        row[1:] = row[1:] - row[-1]
except:
    raise IOError('cannot find models file')

try:
    err_dict = pickle.load( open( join(dirname(__file__), 'err_dict.p'), 'r') )
    ERR_FUNCTIONS = {}
    for mode in [0,1,2]:
        ERR_FUNCTIONS[mode] = {}
        ERR_FUNCTIONS[mode]['range'] = (err_dict[mode]['x'][0], err_dict[mode]['x'][-1])
        for band in err_dict[mode].keys():
            if band == 'x': continue
            ERR_FUNCTIONS[mode][band] = interp1d( err_dict[mode]['x'], err_dict[mode][band] )
except:
    raise IOError('cannot find error dictionary file')


ALL_FILTERS = ['u','g','r','i','z','y','B','V','R','I','J','H','K']
npALL_FILTERS = np.array(ALL_FILTERS) # it's helpful to have a version that works with numpy masks
# band: (central wavelength (AA), zeropoint (erg/s/cm^2/AA), effective width (AA), catalog index)
FILTER_PARAMS =  {'u': (3551., 8.6387e-9, 558.4, 0), 'g': (4686., 4.9607e-9, 1158.4, 1),
                  'r': (6165., 2.8660e-9, 1111.2, 2), 'i': (7481., 1.9464e-9, 1044.5, 3),
                  'z': (8931., 1.3657e-9, 1124.6, 4), 'y': (10091., 1.0696e-9,  1114.8, 5),
                  'B': (4400., 6.6000e-9, 912.8, 6), 'V': (5490., 3.6100e-9, 857.3, 7),
                  'R':(6500., 2.1900e-9, 1320.2, 8),  'I': (7885., 1.1900e-9, 1220.1, 9),
                  'J':(12350., 3.1353e-10, 1624.3, 10), 'H':(16620., 1.1121e-10, 2509.4, 11),
                  'K':(21590., 4.2909e-11, 2618.9, 12)}


############################################
# ONLINE CATALOG MANAGEMENT
############################################

class online_catalog_query():
    '''
    A class to handle all queries of remote catalogs.
    The main function, query_all(), queries all catalogs in a
     multithreaded manner, to minimize lag from I/O communications.
     
    Standard usage:
     q = online_catalog_query( ra, dec, field_size ) #ra,dec in decimal degrees, field_size in arcsec
     Mass, SDSS, USNOB1 = q.query_all()
     # OR #
     Mass = q.query_2mass() #same for all catalogs
     
     If an ignore argument (string or list: ['sdss','apass','usnob']) is included,
      will not attempt to query that catalog.
    '''
    def __init__(self, ra, dec, boxsize=10., ignore=None ):
        self.coords = (ra, dec) #decimal degrees
        self.boxsize = boxsize  #arcseconds
        self.MASS, self.SDSS, self.USNOB, self.APASS = self._query_all( ignore=ignore )
    
    def query_sdss( self ):
        return self.SDSS
    
    def query_2mass( self ):
        return self.MASS
    
    def query_usnob1( self ):
        return self.USNOB
    
    def query_apass( self ):
        return self.APASS
    
    def query_all( self ):
        return self.MASS, self.SDSS, self.USNOB, self.APASS
    
    def _parse_apass( self, s ):
        '''
        Parse an APASS web request string.
        
        s: a string of the CSV returned by a query to the APASS page.
          Example:
          s = urllib2.urlopen('http://www.aavso.org/cgi-bin/apass_download.pl?ra=0.5&dec=85.&radius=.25&outtype=1').read()
        
        Note: it is common for APASS to report -0 as the error for an observation; from their website,
         this merely means that the star was observed only once, and so their error determination
         code does not predict the error correctly.  Currently, this function adopts an error of 0.15mag
         for all of these observations. Missing observations are filled in with a zero.
        '''
        out = []
        for line in s.split('\n')[1:]:
            # discard any sources that have more than one two observations
            if line.count('NA') > 4:
                continue
            try:
                row = line.split(',')
                ra = float(row[0])
                dec = float(row[2])
                # magnitudes and errors
                #  in order: V, V_sig, B, B_sig, g, g_sig, r, r_sig, i, i_sig
                tmp = [ra, dec]
                for obs in row[5:]:
                    if (obs == '-0') or (obs == '0'):
                        # mark missing errors as 0.15
                        tmp.append( 0.15 )
                    elif obs == 'NA':
                        tmp.append( 0.0 )
                    else:
                        tmp.append( np.abs(float(obs)) )
                # reshuffle so that it goes g, r, i, B, V
                tmp = tmp[:2]+tmp[6:]+tmp[4:6]+tmp[2:4]
                out.append( tmp )
            except:
                # silently fail on sources that are not formatted properly,
                #  they are probably anomalous anyways.
                pass
        return np.array(out)
    
    
    def _query_apass( self, container=None, cont_index=3 ):
        '''
        Query apass server for sources found in a box of width self.boxsize (arcsecs)
        around self.coords.
        Returns an array (if objects found) or None (if not)
        
        container: a list into which to put the result; used in multithreaded requests
        cont_index: index of list into which to put the result
        '''
        ra,dec  = self.coords
        boxsize = self.boxsize / 3600. #in degrees
        # search for APASS objects around ra, dec, with a search radius of boxsize,
        #  and return the result in a CSV string
        request = 'http://www.aavso.org/cgi-bin/apass_download.pl?ra={}&dec={}&radius={}&outtype=1'.format(ra, dec, boxsize)
        o = urlopen( request ).read()
        apass_objects = self._parse_apass(o)
        if len(apass_objects) == 0:
            # no matches
            output = None
        else:
            output = apass_objects
        if container == None:
            return output
        else:
            container[ cont_index ] = output
    
    
    def _parse_sdss( self, s ):
        '''
        Parse findsdss8 output.
        
        s: a string as returned by findsdss8 using the flag "-e0,"
        returns: a 2-d array, with each row containing the results for one object:
          [ra, dec, u, u_sigma, g, g_sigma, r, r_sigma, i, i_sigma, z, z_sigma]
        ''' 
        out = []
        lines = [lll for lll in s.split('\n') if lll and lll[0]!='#']
        for line in lines:
            try:
                line = [ val for val in line[4:].split(' ') if val ] #drop first part to avoid inconsistency
                # RA and DEC
                if '+' in line[0]:
                    char = '+'
                else: 
                    char = '-'
                ra  = float(line[0].split(char)[0])
                dec = float(char + line[0].split(char)[1])
                # magnitudes and errors
                #  in order: u, u_sig, g, g_sig, r, r_sig, i, i_sig, z, z_sig 
                tmp = [ra, dec]
                for band in line[3:8]:
                    tmp += map( float, band.split(':') )
                out.append( tmp )
            except:
                # silently fail on sources that are not formatted properly,
                #  they are probably anomalous anyways.
                pass
        return np.array(out)
    
    
    def _query_sdss( self, container=None, cont_index=1, trim_mag=21. ):
        '''
        Query sdss8 server for sources found in a box of width self.boxsize (arcsecs)
        around self.coords.
        Returns an array (if objects found) or None (if not)
        
        container: a list into which to put the result; used in multithreaded requests
        cont_index: index of list into which to put the result
        trim_mag: do not return sources with r > trim_mag
        '''
        ra,dec  = self.coords
        boxsize = self.boxsize
        # search for SDSS objects around coordinates with
        #  defined box size, return only basic parameters, and
        #  sort by distance from coordinates, and return only those
        #  sources brighter than trim_mag
        request = 'findsdss8 -c "{} {}" -bs {} -lmr 0,{} -e0 -lc 4,6 -sr -m 1000000'.format( ra, dec, boxsize, trim_mag )
        #print "FINDSDSS8",ra,dec,boxsize,trim_mag

        out = Popen(request, shell=True, stdout=PIPE, stderr=PIPE)
        o,e = out.communicate()
        if e:
        	if e[:8] != '#...url=':
        		raise IOError('findsdss8 problem: '+e)
        # parse the response
        sdss_objects = self._parse_sdss(o)
        if len(sdss_objects) == 0:
            # no matches
            output = None
        else:
            sdss = np.array( [obj for obj in sdss_objects if obj[6] < trim_mag] )
            output = sdss
        if container == None:
            return output
        else:
            container[ cont_index ] = output
    
    def _parse_2mass( self, s ):
        '''
        parse find2mass output.
        
        s: a string as returned by find2mass using the flag "-eb"
        returns: a 2-d array, with each row containing the results for one object:
          [ra, dec, J, J_sigma, H, H_sigma, K, K_sigma]
        '''
        out = []
        lines = [lll for lll in s.split('\n') if lll and lll[0]!='#']
        for line in lines:
            try:
                line = line.split('|')
                # RA and DEC
                ra, dec = map( float, line[0].split(' ') )
                # magnitudes and errors
                #  in order: J, J_sig, H, H_sig, K, K_sig
                tmp = [ra, dec]
                for band in line[2:5]:
                    mag, err = [val for val in band.split(' ') if val]
                    mag = float(mag)
                    if '-' not in err:
                        # some objects do not have reported errors - for these,
                        #  assume a conservative error of 0.25mag
                        err = float(err)
                    else:
                        err = .25
                    tmp += [mag, err]
                out.append( tmp )
            except:
                # silently fail on sources that are not formatted properly,
                #  they are probably anomalous anyways.
                pass
        return np.array(out)
    
    
    def _query_2mass( self, container=None, cont_index=0 ):
        '''
        Query 2mass server for sources found in a box of width self.boxsize (arcsecs)
        around self.coords.
        Returns an array (if objects found) or None (if not)
        
        container: a list into which to put the result; used in multithreaded requests
        cont_index: index of list into which to put the result
        '''
        ra,dec  = self.coords
        boxsize = self.boxsize
        # search for 2Mass point sources around coordinates with
        #  defined box size, return only basic parameters, and 
        #  sort by distance from coordinates, and return a 
        #  maximum of 1000000 sources (i.e. return everything)
        request = 'find2mass -c {} {} -bs {} -eb -sr -m 1000000'.format( ra, dec, boxsize )
        out = Popen(request, shell=True, stdout=PIPE, stderr=PIPE)
        o,e = out.communicate()
        if e:
        	if e[:8] != '#...url=':
        		raise IOError('find2mass problem: '+e)
        # parse the response
        mass_objects = self._parse_2mass(o)
        if len(mass_objects) == 0:
            # no matches
            output = None
        else:
            output = mass_objects
        if container == None:
            return output
        else:
            container[ cont_index ] = output
    
    
    def _parse_usnob1( self, s ):
        '''
        Parse findusnob1 output.
        The photometric errors in USNOB1 are pretty terrible (~.3mag for each observation),
          but most sources have more than one observation, so I average all results
          and return the average magnitudes and the error of the mean.
        
        s: a string as returned by usnob1 using the flag "-eb"
        returns: a 2-d array, with each row containing the results for one object:
          [ra, dec, avg_B, B_sigma, avg_R, R_sigma, I, I_sigma]
        '''
        # parse the header to see how many of each magnitude are reported
        header = s.split('\n')[3]
        obs_count = header.count('Bmag')
        # just as a sanity check, make sure it reports the same number of B,R mags
        assert( obs_count == header.count('Rmag') )
        
        out = []
        lines = [lll for lll in s.split('\n') if lll and lll[0]!='#']
        for line in lines:
            try:
                line = line.split('|')
                # RA and DEC
                tmp = [ val for val in line[0].split(' ') if val ]
                if '+' in tmp[1]:
                    char = '+'
                else: 
                    char = '-'
                ra  = float(tmp[1].split(char)[0])
                dec = float(char + tmp[1].split(char)[1])
            
                # magnitudes and errors
                #  in order: B, B_sigma, R, R_sigma
                Bs, Rs = [], []
                for i in range(obs_count):
                    tmp = line[1 + 2*i]
                    if '-' not in tmp:
                        Bs.append( float(tmp) )
                    tmp = line[2 + 2*i]
                    if '-' not in tmp:
                        Rs.append( float(tmp) )
                tmp = line[1 + 2*obs_count]
                if '-' not in tmp:
                    # NOTE: results are consistently better without the I-band photometry
                    #I = float(tmp)
                    #I_err = 0.3
                    I = I_err = 0
                else:
                    I = I_err = 0
                        
                # ignore sources that don't have at least two observations
                # each measure has an error of about .3 mag
                if sum( [not Bs, not Rs, not I] ) > 1: continue
                if Bs:
                    B = np.mean(Bs)
                    B_err = 0.3/np.sqrt(len(Bs))
                else:
                    B = B_err = 0
                if Rs:
                    R = np.mean(Rs)
                    R_err = 0.3/np.sqrt(len(Rs))
                else:
                    R = R_err = 0
                out.append( [ra, dec, B, B_err, R, R_err, I, I_err] )
            except:
                # silently fail on sources that are not formatted properly,
                #  they are probably anomalous anyways.
                pass
        
        return np.array(out)
    
    
    def _query_usnob1( self, container=None, cont_index=2 ):
        '''
        Query usnob1 server for sources found in a box of width self.boxsize (arcsecs)
        around self.coords.
        Returns an array (if objects found) or None (if not)
        
        container: a list into which to put the result; used in multithreaded requests
        cont_index: index of list into which to put the result
        '''
        ra,dec  = self.coords
        boxsize = self.boxsize
        # search for USNOB1 point sources around coordinates with
        #  defined box size, return only basic parameters, and 
        #  sort by distance from coordinates, and return a maximum
        #  of 10000 objects (i.e. return everything)
        request = 'findusnob1 -c {} {} -bs {} -eb -sr -m 1000000'.format( ra, dec, boxsize )
        out = Popen(request, shell=True, stdout=PIPE, stderr=PIPE)
        o,e = out.communicate()

        if e:
        	if e[:8] != '#...url=':
        		raise IOError('findusnob1 problem: '+e)
        usnob1_objects = self._parse_usnob1(o)
        if len(usnob1_objects) == 0:
            # no matches
            output = None
        else:
            output = usnob1_objects
        if container == None:
            return output
        else:
            container[ cont_index ] = output
    
    
    def _query_all( self, ignore=None ):
        '''
        Query all sources, with an independent thread for each
         so that the communications happen concurrently, to save
         time.
        
        returns: [2Mass, SDSS, USNOB1, APASS]
        '''
        ra,dec  = self.coords
        boxsize = self.boxsize
        # results is a container into which the threads will put their responses
        results = [None]*4
        # only query the ones we want
        t0 = Thread(target=self._query_2mass, args=(results,))
        threads = [t0]
        if ignore==None or 'sdss' not in ignore:
            t1 = Thread(target=self._query_sdss, args=(results,))
            threads.append(t1)
        if ignore==None or 'usnob' not in ignore:
            t2 = Thread(target=self._query_usnob1, args=(results,))
            threads.append(t2)
        if ignore==None or 'apass' not in ignore:
            t3 = Thread(target=self._query_apass, args=(results,))
            threads.append(t3)
        for t in threads:
            t.start()
        for t in threads:
            # join each thread and wait until each is completed
            t.join()
        return results
    


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


def find_field( coords, extend=3. ):
    '''
    Determine the best field for a list of star coordinates,
    so we can perform only a single query.
    
    star_coords: a list/array of star coordinates, (decimal degrees)
    extend: the buffer beyond the requested coordinates to add to the field in both directions
      (in arcseconds)
    
    returns: (coordinates of center in decimal degrees), (RA_width_of_box, DEC_width_of_box) (both in arcseconds)
    '''
    coords = np.array(coords) #make sure the input coords are an array
    
    # handle RA rollovers
    if any( coords[:,0] > 359. ) and any( coords[:,0] < 1. ):
        # put all coordinates on the high side, will be rolled back at end
        mask = coords[:,0]<1.
        coords[:,0][ mask ] = coords[:,0][ mask ] + 360.
    
    # the field widths in each direction
    r_max = np.deg2rad( np.max(coords[:,0]) )
    r_min = np.deg2rad( np.min(coords[:,0]) )
    d_max = np.deg2rad( np.max(coords[:,1]) )
    d_min = np.deg2rad( np.min(coords[:,1]) )
    
    # the field center
    r_c = r_min + (r_max-r_min)/2.
    d_c = d_min + (d_max-d_min)/2.
    
    # the conversions between angles in RA,Dec and actual angles
    X_ra = lambda r,d: (np.cos(d)*np.sin(r-r_c))/(np.cos(d_c)*np.cos(d)*np.cos(r-r_c) +\
                                              np.sin(d_c)*np.sin(d))
    Y_dec = lambda r,d: (np.cos(d_c)*np.sin(d)-np.cos(d)*np.sin(d_c)*np.cos(r-r_c))/ \
                   (np.sin(d_c)*np.sin(d)+np.cos(d_c)*np.cos(d)*np.cos(r-r_c))
    
    # the actual widths in rads
    X_width = np.abs(X_ra(r_max, d_c) - X_ra(r_min, d_c))
    Y_width = np.abs(Y_dec(r_c, d_max) - Y_dec(r_c, d_min))
    
    # convert center to decimal degrees and make sure to rollback any RA rollovers (if they happened)
    out_c = np.rad2deg( [r_c, d_c] ).tolist()
    out_c[0] = out_c[0]%360.
    # convert widths to arcseconds and extend them
    widths = (np.rad2deg([X_width, Y_width])*3600 + 2*extend).tolist()
    return out_c, widths


def split_field( field_center, field_width, nsplit ):
    '''
    Split a large field (of width field_width) into nsplit smaller fields,
     tiled to fill a square of edgesize field_width.
    field_center in decimal degrees, field_width in arcseconds
    Returns an array of pointings and the distances between field centers (all in decimal degrees)
    '''
    fw = field_width/3600. # in degreees
    fw = np.deg2rad(fw)
    # assume dec is flat
    dec = np.deg2rad(field_center[1])
    decs = np.array([ (dec-fw/2 + (fw/nsplit/2)) + (fw/nsplit)*i for i in range(nsplit) ])
    # do actual trig to find ra
    ra = np.deg2rad(field_center[0])
    dr = (fw/nsplit)/np.cos(dec)
    ras = np.array([ (ra-(fw/2)/np.cos(dec) +(fw/nsplit/2)/np.cos(dec)) + dr*i for i in range(nsplit) ])
    centers = []
    for i in range(len(ras)):
        for j in range(len(decs)):
            centers.append( [ras[i],decs[j]] )
    return np.rad2deg(centers), np.rad2deg( [dr, fw/nsplit] )

    


############################################
# MODEL-FITTING FUNCTIONS
############################################

def _error_C(C, model, obs, weights):
    '''
    C: number, a constant akin to distance modulus
    model: array-like, model mags
    obs: array-like, observed mags
    weights: array-like, 1/(observation errors)
    
    returns: weighted sum squared errors
    '''
    nm = model+C
    return np.sum( weights*(nm-obs)**2 )

def choose_model( obs, mask, models=MODELS, allow_cut=False ):
    '''
    Find and return the best model for obs.
    Do this by fitting to all magnitudes weighted by error.
    
    Returns: model, temperature, quality_parameter
    
    obs: an array of the observed magnitudes and errors, in
           order as defined by the mode key (see below)
    mask: defines what colors to use in the fit, i.e. what observations exist
           in order: [u,g,r,i,z,y,B,V,R,I,J,H,K]
    models: an array of modeled SEDs, where 0th entry is temperature
             of the model, and the rest are magnitudes
    allow_cut: if true, allows one datapoint to be cut from fit
    '''
    # mask is an array used to choose which modeled magnitudes
    #  correspond to the included observations.
    #  Note: I append a leading zero to ignore the temperature of the model
    #   (saved in position 0 for all models)
    mask = np.hstack( (np.array([False]), np.array(mask).astype(bool)) )
    
    mags = obs[::2]
    Jmag = mags[-1]
    zerod_mags = mags - Jmag # recenter to compare to models
    weights = 1./obs[1::2]
    
    # Go through all models and choose the one with the most similar SED
    #  Keep track of the error, as returned by _error_C()
    sum_sqrs, Cs = [], []
    for model in models[1:]:
        res = fmin_bfgs( _error_C, 0., args=(model[mask], zerod_mags, weights), full_output=True, disp=False )
        Cs.append(res[0][0])
        sum_sqrs.append(res[1])
    
    i_best = np.argmin(sum_sqrs)
    best_model = models[1:][ i_best ]
    
    i_cut = None
    if allow_cut:
        max_diff = 1.
        # if there's one point more than <max_diff> away from best model, try again without that point
        #  (set that weight to zero, and return the index of the cut value)
        if max( np.abs(best_model[mask] - zerod_mags) ) > max_diff:
            i_cut = np.argmax(np.abs(best_model[mask] - zerod_mags))
            weights[i_cut] = 0.
            # Go through all models again, with new weights
            sum_sqrs, Cs = [], []
            for model in models[1:]:
                res = fmin_bfgs( _error_C, 0., args=(model[mask], zerod_mags, weights), full_output=True, disp=False )
                Cs.append(res[0][0])
                sum_sqrs.append(res[1])
                
            i_best = np.argmin(sum_sqrs)
            best_model = models[1:][ i_best ]
    
    # now add back in the Jmag value to get a model for the non-zeroed observations
    C = Cs[i_best] + Jmag
    
    # return all magnitudes for best model, the offset C, the index, and a quality metric for the best fit
    #  The quality metric is the reduced Chi^2 statistic (assuming model + distance are two fitting parameters).
    sum_sqr_err = sum_sqrs[i_best]
    if i_cut != None:
        metric = sum_sqr_err/(len(mags)-1 - 2)
    else:
        metric = sum_sqr_err/(len(mags) - 2)
    return (best_model[1:] + C, C, best_model[0], metric, i_cut )

def fit_sources( inn, f_err=ERR_FUNCTIONS, return_cut=False ):
    '''
    Wrapper function for choose_model to facilitate multiprocessor use in catalog.produce_catalog()
    
    If return_cut == True, any observations that are ignored in the fitting
     procedure are replaced by the modeled result.  Otherwise, observed
     photometry is always returned, whether used in model fit or not.
    '''
    # the masks show, in order, which bands are included
    #  order: u,g,r,i,z, y, B,V,R,I, J,H,K
    mode, obs = inn
    if mode == 0: # sdss+2mass
        mask = [1,1,1,1,1, 0, 0,0,0,0, 1,1,1]
        allow_cut = True
    elif mode == 1: # apass+2mass
        mask = [0,1,1,1,0, 0, 1,1,0,0, 1,1,1]
        # we allow 3 - 5 APASS observations, so make sure to handle that correctly
        for i,imask in enumerate([1,2,3,6,7]):
            if obs[2*i] == 0:
                mask[imask] = 0
        obs = obs[ obs>0 ]
        allow_cut = True
    elif mode == 2: # usnob+2mass
        mask = [0,0,0,0,0, 0, 1,0,1,1, 1,1,1]
        # we can allow for 2 or 3 usnob obs, so make sure we handle that correctly
        for i,imask in enumerate([6,8,9]):
            if obs[2*i] == 0:
                mask[imask] = 0
        obs = obs[ obs>0 ]
        allow_cut = True
        
    mask = np.array(mask).astype(bool)
    model, C, index, err, i_cut = choose_model( obs, mask, allow_cut=allow_cut )
    
    sed = np.empty(len(mask))
    full_errs = np.empty(len(mask))
    # keep the real observations
    sed[mask] = obs[::2]
    full_errs[mask] = obs[1::2]
    # fill in rest with modeled magnitudes
    sed[~mask] = model[~mask]
    # and return errors as estimated by model Chi^2
    for band in npALL_FILTERS[~mask]:
        i_band = FILTER_PARAMS[band][-1]
        # temporary fill-in for when USNOB bands are gone or get cut
        if (mode == 2) and (band in ['B','R']):  #note: assuming we always model I
            full_errs[i_band] = 0.5
        # temporary fill-in for when APASS bands are gone or get cut
        elif (mode == 1) and (band in ['B','V','g','r','i']):
            full_errs[i_band] = 0.25
        else:
            full_errs[i_band] = f_err[mode][band]( min(err, f_err[mode]['range'][1]) )
    
    # if a value was cut while fitting, return the modeled magnitude instead of the observed
    if return_cut and i_cut != None:
        # do some gymnastics to get the cut passband since it's behind a mask
        cut_band = npALL_FILTERS[mask][i_cut]
        i_cut_band = FILTER_PARAMS[cut_band][-1]
        sed[i_cut_band] = model[mask][i_cut]
        # For now simply report a (conservative) error of 0.5, though this is wrong.
        # Should run a series of tests seeing how good each observed band (in each mode)
        #  is predicted by the rest of the observed bands, and put those results into f_err as below
        #full_errs[i_cut_band] = f_err[mode][cut_band]( min(err, f_err[mode]['range'][1]) )
        full_errs[i_cut_band] = 0.5
        
    return ( sed, full_errs, err, index )


############################################
# SYNTHETIC CATALOG MANAGEMENT
############################################

class catalog():
    '''
    A class to handle all catalog creation.
    The main function produces, for any requested field,
     a catalog of cross-matched sources from
     2-MASS and (SDSS or USNOB1), with any missing
     bands filled in through modeling.
    
    Standard usage:
     c = catalog( (ra, dec), field_size )  # with ra, dec in degrees and field_size in arcseconds
     catalog_coords = c.coords
     catalog_SEDs = c.SEDs
    
    Optional arguments:
     input_coords: if given, will only attempt to produce a catalog for these sources.
     ignore: option to ignore any of the bands but 2mass. Can be list or single item, options
      are:  ["sdss", "usnob", "apass"]
    '''
    MAX_SIZE = 7200 # max size of largest single query
    ERR_CUT  = (8., 5., 2.5)   # maximum reduced chi^2 to keep a fit (SDSS, APASS, USNOB)
    
    def __init__( self, field_center, field_width, input_coords=None, ignore=None ):
        self.field_center = field_center
        self.field_width = field_width
        if input_coords != None:
            self.input_coords = np.array(input_coords)
        else:
            self.input_coords = input_coords
        self.coords = []
        self.SEDs = []
        self.full_errors = []
        self.model_errors = []
        self.models = []
        self.modes = []
        
        #VLT
        self.cmodes = []
        self.serr = []
        self.scat = []
        self.ccoords = []
        
        self.numcut = 0
        self.bands = ALL_FILTERS
        # this switch controls whether we ignore any catalogs.
        #  can be list or string in ['usnob','apass','sdss']
        self.ignore = ignore
        if field_width > self.MAX_SIZE:
            # simply don't allow queries that are too large
            raise ValueError( 'Field is too large. Max allowed: {}"'.format(self.MAX_SIZE) )
        else:
            self.produce_catalog()
    
    
    def produce_catalog( self ):
        '''
        Create a catalog of all objects found in field.
        Requires records in 2MASS + (SDSS and/or USNOB1).
        '''
        ra,dec = self.field_center
        q = online_catalog_query( ra, dec, self.field_width, ignore=self.ignore )
        mass, sdss, usnob, apass = q.query_all()
        
        object_mags = []
        modes = []
        object_coords = []
        cmode = -1
        if mass != None:
            if self.input_coords != None:
                # if input coordinates were given, ignore all other objects
                input_matches, tmp = identify_matches( mass[:,:2], self.input_coords )
                if not np.sum( np.isfinite(tmp) ):
                    raise ValueError( 'No matches to input objects found!' )
                mass = mass[ input_matches>=0 ]
                
            # match sdss, apass, usnob objects to 2mass objects
            if sdss != None:
                sdss_matches, tmp = identify_matches( mass[:,:2], sdss[:,:2] )
                if cmode == -1:
                	cmask = [0,-1,1,-1,2,-1,3,-1,4,-1]
                	cmode = 0
                	ocat = sdss              
            else:
                sdss_matches = -9999*np.ones(len(mass), dtype='int')
                
            if apass != None:
                apass_matches, tmp = identify_matches( mass[:,:2], apass[:,:2] )
                if cmode == -1:
                	cmask = [1,-1,2,-1,3,-1,6,-1,7,-1] 
                	cmode = 1
                	ocat = apass           	
            else:
                apass_matches = -9999*np.ones(len(mass), dtype='int')
            if usnob != None:
                usnob_matches, tmp = identify_matches( mass[:,:2], usnob[:,:2] ) 
                if cmode == -1:
                	cmask = [6,-1,8,-1]
                	cmode = 2
                	ocat = usnob
            else:
                usnob_matches = -9999*np.ones(len(mass), dtype='int')
            
            #VLT extra optical catalog objects
            ccords = ocat[:,:2]
            scat = np.zeros([len(ocat),13])+99
            serr = np.zeros([len(ocat),13])+9
            for p,val in enumerate(cmask):
                if val != -1:
                	scat[:,val] = ocat[:,p+2]
                	serr[:,val] = ocat[:,p+3]
        	
        	self.scat = scat
        	self.serr = serr
        	self.ccoords = ccords
        	self.cmodes = np.zeros(len(scat), dtype=int) + cmode
                	          
            # Go through 2mass objects and assemble a catalog
            #  of all objects present in 2mass and (sdss or apass or usnob)
            #  Preference ranking: 2MASS + (SDSS > APASS > USNOB)
            for i,obj in enumerate(mass):
                if (sdss_matches[i]>=0):
                    i_sdss = sdss_matches[i]
                    obs = np.hstack( (sdss[i_sdss][2:], obj[2:]) )
                    mode = 0
                elif (apass_matches[i]>=0):
                    i_apass = apass_matches[i]
                    obs = np.hstack( (apass[i_apass][2:], obj[2:]) )
                    mode = 1
                elif (usnob_matches[i]>=0):
                    i_usnob = usnob_matches[i]
                    obs = np.hstack( (usnob[i_usnob][2:], obj[2:]) )
                    mode = 2
                else:
                    continue
                object_mags.append( obs )
                modes.append( mode )
                object_coords.append( obj[:2] )
        if len(object_coords) < 1:
            raise ValueError( "No good sources in this field!" )
        

        # send all of these matches to the CPU pool to get modeled
        objects = zip( modes, object_mags )
        pool = mp.Pool( processes=N_CORES )
        results = pool.map( fit_sources, objects )
        #results = [fit_sources(obj) for obj in objects]
        pool.close()
        pool.join()
        
        # now go through results and construct the final values
        for i,row in enumerate(results):
            # each row is (sed, full_errs, model_err, index)
            if row[2] > self.ERR_CUT[modes[i]]: #apply quality-of-fit cut
                self.numcut += 1
                pass
            else:
                self.coords.append( object_coords[i] )
                self.SEDs.append( row[0] )
                self.full_errors.append( row[1] )
                self.model_errors.append( row[2] )
                self.models.append( row[3] )
                self.modes.append( modes[i] )
        # keep the multi-dimensional data in numpy arrays
        self.coords = np.array(self.coords)
        self.SEDs = np.array(self.SEDs)
        self.full_errors = np.array(self.full_errors)
        print 'cut', self.numcut, 'sources out of', len(results)
        
        gmatch, tmp = identify_matches( self.coords, self.ccoords )
        
        for ind,match in enumerate(gmatch):
        	if match > 0:
        		self.ccoords[match] = self.coords[ind]
        		self.scat[match] = self.SEDs[ind]
        		self.serr[match] = self.full_errors[ind]
        		self.cmodes[match] = self.modes[ind]
    
    def save_catalog( self, file_name ):
        save_catalog( self.coords, self.SEDs, self.full_errors, self.modes, file_name )
    
    def get_source( self, ra, dec ):
        '''
        Return a dictionary with the coordinates and photometry for the object
         at ra, dec (if that object is matched in the catalog).  Returns None
         if no match found.
        output = { 'coords':[ra, dec], 'y':[y_mag, y_err], 'z':[z_mag, z_err], ... }
        '''
        index, dist = identify_matches( [[ra,dec]], self.coords )
        index = index[0] #take the best match
        if index < 0:
            # no match found
            return None
        else:
            out_dict = { 'coords':self.coords[index] }
            for i,band in enumerate(self.bands):
                out_dict[band] = [ self.SEDs[index][i], self.full_errors[index][i] ]
            return out_dict
    

def save_catalog( coordinates, seds, errors, modes, file_name ):
    '''
    Save an output ASCII file of the catalog.
    file_name: output file to create
    '''
    
    fff = open(file_name,'w')
    fff.write('# Observed/modeled SEDs produced by get_SEDs.py \n' +
              '# Generated: {}\n'.format(strftime("%H:%M %B %d, %Y")) +
              '#  Mode = 0: -> B,V,R,I,y modeled from SDSS and 2-MASS\n' +
              '#       = 1: -> u,y,R,I modeled from APASS and 2-MASS\n' +
              '#       = 2: -> u,g,r,i,z,y,V,I modeled from USNOB-1 and 2-MASS\n' +
              "# " + "RA".ljust(10) + "DEC".ljust(12) + "".join([f.ljust(8) for f in ALL_FILTERS]) + \
              "".join([(f+"_err").ljust(8) for f in ALL_FILTERS]) + "Mode\n")
    for i,row in enumerate(seds):
        row_txt = "".join([ s.ljust(12) for s in map(lambda x: "%.6f"%x, coordinates[i]) ]) +\
                  "".join([ s.ljust(8) for s in map(lambda x: "%.3f"%x, row) ]) +\
                  "".join([ s.ljust(8) for s in map(lambda x: "%.3f"%x, errors[i]) ]) +\
                  str(modes[i])+"\n"
        fff.write( row_txt )
    fff.close()


############################################
# ZEROPOINT CALCULATION
############################################

def clip_me( inn, sig_clip=3., max_iter=5, convergence=.02 ):
    '''
    Perform iterative sigma clipping about the median of inn (a numpy array).
    
    sig_clip: iteratively trim values more than sig_clip*std away from median
    max_iter: stop the above sigma clipping after max_iter iterations
    convergence: the fractional convergence required
    '''
    for count in range(max_iter):
        in_len = len( inn )
        inn = inn[ np.abs(inn-np.median(inn)) < sig_clip*np.std(inn) ]
        if float(in_len-len(inn))/in_len < convergence:
            break
    return inn


def calc_zeropoint( input_coords, catalog_coords, input_mags, catalog_mags, clip=True, return_zps=True ):
    '''
    Calculate the zeropoint for a set of input stars and set of catalog stars.
    
    input_coords: a 2D array of [ [RA, DEC], [..., ...] ... ]
    catalog_coords: similar array for objects created with catalog()
    input_mags: a 1D array of instrumental magnitudes
    catalog_mags: a similar array of true magnitudes as created with catalog()
    clip: perform sigma clipping about median on zeropoint array
    return_zps: return all zeropoint estimates, not just the median
    
    Returns: zeropoint (mags), the median average deviation, and a list of matched indices for input and catalog sources.
    '''
    # make sure all inputs are numpy arrays
    input_coords = np.array(input_coords)
    catalog_coords = np.array(catalog_coords)
    input_mags = np.array(input_mags)
    catalog_mags = np.array(catalog_mags)
    
    matches, tmp = identify_matches( input_coords, catalog_coords )
    matched_inputs = input_mags[ matches>=0 ]
    matched_catalogs = catalog_mags[ matches[ matches>=0 ] ]
    
    zp_estimates = []
    for i,inst_mag in enumerate(matched_inputs):
        zp_estimates.append( matched_catalogs[i] - inst_mag )
    zp = np.array(zp_estimates)
    if clip:
        zp = clip_me( zp )
    mad = np.median( np.abs( zp-np.median(zp) ) )
    if return_zps:
        return np.median(zp), mad, matches, zp
    else:
        return np.median(zp), mad, matches
        


def zeropoint( input_file, band, output_file=None, usnob_thresh=15, alloptstars=False, quiet=False):
    '''
    Calculate <band> zeropoint for stars in <input_file>.
    
    Expects a space-or-tab-delimited ascii input file with the
     first column RA, the second DEC, and the third instrumental magnitude.
     Header/comments should be #-demarcated, and all non-commented rows in the
     file should be numbers only.
    If an output_file name is given, it saves the entire catalog to that file.
    usnob_thresh: the minium number of APASS+SDSS sources required before starting to use USNOB sources
    '''
    
    if quiet == 'False': quiet = False
    if alloptstars == 'False': alloptstars = False
    usnob_thresh = int(usnob_thresh)

    # load the data and produce a catalog
    in_data = np.loadtxt( input_file )
    input_coords = in_data[:, :2]
    input_mags = in_data[:, 2]
    field_center, field_width = find_field( input_coords )
    c = catalog( field_center, max(field_width), input_coords=input_coords )
    
    band_index = FILTER_PARAMS[band][-1]
    # check to see whether we need to use USNOB sources
    mask = np.array(c.modes)<2

    if sum( mask ) >= usnob_thresh:
        if quiet == False: print 'Using',sum(mask),'APASS and/or SDSS sources.'
        cat_mags = c.SEDs[:, band_index][mask]
        cat_coords = c.coords[mask]
    else:
        if quiet == False: print 'Using',sum(mask),'USNOB, APASS, and/or SDSS sources.'
        cat_mags = c.SEDs[:, band_index]
        cat_coords = c.coords
        
    if alloptstars:
    	zp, mad, matches, ze = calc_zeropoint( input_coords, c.ccoords, input_mags, c.scat[:,band_index], return_zps=True )
    	errors = c.serr
    	catmag = c.scat
    	modes  = c.cmodes
    else:
    	zp, mad, matches = calc_zeropoint( input_coords, cat_coords, input_mags, cat_mags, return_zps=False )
    	errors = c.full_errors
    	catmag = c.SEDs
    	modes  = c.modes
    	
    # save matched catalog to file
    if output_file:
    	oc, os, oe, om = [],[],[],[]
    	for i,match in enumerate(matches):
    		oc.append( input_coords[i] )
    		if match >= 0:
    			os.append( catmag[match] )
    			oe.append( errors[match] )
    			om.append( modes[match] )
    		else:
    			os.append( [99]*len(ALL_FILTERS) )
    			oe.append( [9]*len(ALL_FILTERS) )
    			om.append( -1 )
    	save_catalog( oc, os, oe, om, output_file )
    	
    return zp, mad
    
if __name__ == "__main__":
	warnings.filterwarnings('ignore')
	
	zeropoint(*sys.argv[1:])