# -*- mode: python; coding: utf-8 -*-

"""
NAME:
	printhtml
PURPOSE:
	Create photcomp.plot.png comparing magnitude and errors as well as create HTML page that has
	all of the data easily displayed for up to 9 filters.
	Plot SEDs for each source and all filters.
OUTPUT:
	photcomp.plot.png - shows magnitude vs. error for each filter
	abmag-energy.plot.png - shows flux vs. magnitude for each filter
	photom.html   - html page showing information about sources
	seds/(INDEX).plot.png - shows plot of SEDs for all filters of each source
DEPENDENCIES:
	Images and files (finalmags.txt) created from plotphotom.py

Translated from printratirhtml.pro by John Capone (jicapone@astro.umd.edu).
Modified 7/31/2014 by Vicki Toy (vtoy@astro.umd.edu)
"""

import glob
import numpy as np
#import pylab as pl
import photprocesslibrary as pplib

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl
# Disable interactive mode
pl.ioff()

from datetime import datetime

def plot_mag(filters, plotdict, colors):
	pl.figure()
	pl.xlim([15,22])
	pl.ylim([-0.01,0.2])
	
	#For each filter, plot mag vs. error for photocomp.png
	counter = 0
	for filter in filters:
		pl.plot(plotdict[filter+'mag'], plotdict[filter+'magerr'], marker='o',linestyle='None', label=filter, color=colors[counter])
		counter = counter + 1
	
	pl.xlabel('AB Magnitude')
	pl.ylabel(r"$\Delta$ Mag")
	pl.legend(loc='lower right')
	pl.savefig('photcomp.plot.png', bbox_inches='tight')
	pl.clf()

def plot_seds(filters):
    wavelengths = { 'r': 6122.3, 'i': 7439.5, 'z': 8897.1, 'y': 10289.4, 'J': 12350.0, 'K': 21590.0, 'u': 3594.9, 'g': 4640.4, 'H':16620.0 }
    file_mask = ['RA', 'DEC', 'u', 'g', 'r', 'i', 'z', 'y', 'B', 'V', 'R', 'I', 'J', 'H', 'K', 'u_err', 'g_err', 'r_err', 'i_err', 'z_err', 'y_err', 'B_err', 'V_err', 'R_err', 'I_err', 'J_err', 'H_err', 'K_err', 'Mode', 'Timestamp', 'Filter']
    file_formats = ('<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', 'U18', 'U1')

    files = glob.glob('seds/*.seds.txt')

    if len(files) == 0:
        print 'Did not find any SEDs files! Check your data directory path!'
        return    

    figure = pl.figure()
    pl.ylim([15,22])
    pl.xlim([3000, 22000])
    
    for file in files:
        source_index = (file.split(".")[0]).split('/')[1]

        # data = np.genfromtxt(file, delimiter=' ', dtype=None)
        data = np.loadtxt(file, dtype={'names': file_mask, 'formats': file_formats}, delimiter=' ', skiprows=0)
        filters_map = {}

        if data.size == 1:
            data = np.array([data])

        for i in range(data.size):
            filter = (data[i])[file_mask.index('Filter')]
            filters_map[filter] = np.array([data[i]]) if not (filter in filters_map) else np.append(filters_map[filter], np.array([data[i]]), 0)

        ncols = (2 if len(filters_map) > 1 else 1)
        nrows = (int(len(filters_map) / 2) if len(filters_map) > 1 else 1)
        fig, ax = pl.subplots(nrows=nrows, ncols=ncols, figsize=(6 * ncols, 3 * nrows + 1))
        # pl.subplots_adjust(top=0.8, hspace=0.3, wspace=0.3) with top title
        pl.subplots_adjust(hspace=0.3, wspace=0.3)

        ax_row = 0
        ax_col = 0 if len(filters_map) > 2 else 1

        for filter, seds in filters_map.iteritems():
            for i, sed in enumerate(seds):
                plot_x = []
                plot_y = []

                for f in filters:
                    plot_x.append(wavelengths[f])
                    plot_y.append(sed[file_mask.index(f)])

            col = (ax[ax_row])[ax_col] if len(filters_map) > 2 else (ax[ax_row] if len(filters_map) > 1 else ax)    

            col.set_title('Filter ' + filter)
            col.set_xlabel(r'Wavelength ($\AA$)')
            col.set_ylabel('Magnitude')
            
            plot_label = datetime.strptime(sed[file_mask.index('Timestamp')],'%Y%m%dT%H%M%S%f').strftime('%Y-%m-%d %H:%M:%S')
            col.plot(plot_x, plot_y, label=plot_label)
            
            # https://stackoverflow.com/a/4701285
            # Shrink current axis's height by 10% on the bottom
            box = col.get_position()
            col.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])
            # Put a legend below current axis
            col.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True, shadow=False, ncol=5, prop={'size': 11})

            if (ax_col + 1) % 2 == 0:
                ax_row = ax_row + 1
                ax_col = 0
            else:
                ax_col = ax_col + 1            

        # pl.suptitle('Source #' + source_index, fontsize=16)        
        pl.savefig('seds/' + source_index + '.plot.png', bbox_inches='tight')
        pl.clf()
        pl.close(fig)

    pl.close(figure)

def plot_seds_test():
    filters = ['r','i','z','y','J','H']
    plot_seds(filters)

def generate_headers_table(headers, headers_names):
	t = '<table border="1"><tr>'
	for name in headers_names:
		t = t + '<th>' + name['title'] + '</th>'
	t = t + '</tr>'

	for header in headers:
		row = '<tr>'
		for name in headers_names:
			header_value = name['format'] % header[name['key']] if 'format' in name else header[name['key']] 
			row = row + '<td>' + header_value + '</td>'
		row = row + '</tr>'

		t = t + row

	t = t + '</table>'

	return t

def printhtml(filters, colnames, omitted_colnames, headers, headers_names):
	#Reads in final magnitudes
	colgrab = np.loadtxt('./finalmags.txt', unpack=True)
	
	#Store each column into dictionary based on colnames and create header for HTML header from colnames
	plotdict = {}
	t = '<tr><th>#</th>'
	for i in np.arange(len(colnames)):
		plotdict[colnames[i]] = colgrab[i,:]
        if not (colnames[i] in omitted_colnames):
            t = t + '<th>' + colnames[i]+'</th>'
    t = t + '</tr>\n' 

	colors = ['black', 'purple', 'blue', 'aqua', 'green', 'orange', 'red', 'yellow', 'magenta']
	print 'Plotting AB Magnitude comparison'
    plot_mag(filters, plotdict, colors)
    print 'Plotting SEDs'
	plot_seds(filters)

	#Create html page that displays images made in plotphotom.py and values from finalmags.txt
	f = open( './photom.html', 'w' )
	f.write( '<!DOCTYPE HTML>\n' )
	f.write( '<HTML>\n' )
	f.write( '<HEAD>\n' )
	f.write( '<TITLE>PHOTOMETRY DATA</TITLE>\n' )
	f.write( '</HEAD>\n' )
	f.write( '<BODY BGCOLOR="#FFFFFF" TEXT="#003300">\n' )

	#Finds filter image files with green circles over sources and writes to HTML
	prefchar = 'coadd'
	zffiles = pplib.choosefiles( prefchar + '*_?.png' )

	im_wid = 400

	f.write( '<IMG SRC="./color.png" width="' + `im_wid` + '"><BR>\n' )
	for i in range(len(zffiles)):
	    f.write( '<IMG SRC="./' + zffiles[i] + '" width="' + `im_wid` + '">\n' )
	    if (i+1)%3 == 0:
	    	f.write( '<BR>\n' )

	#Write table with header info from each image
	f.write( '<BR>\n' )
	f.write(generate_headers_table(headers, headers_names))
	f.write( '<BR>\n' )

	f.write( '<BR><HR><FONT SIZE="+2" COLOR="#006600">AB System Photometry (sources within 1 arcmin):</FONT><BR>\n' )
	f.write( 'Notes: Non-zero magnitudes with uncertainty of zero are 3-sigma upper limits.  Sources with magnitudes of 0.0000 are unobserved.<BR>\n' )
	f.write( '		 Circles in images above match aperture size used in sextractor.<BR>\n' )
	f.write( '<BR>\n' )


	#Writes table with all magnitudes
	f.write( '<table border="1" width="100%">\n' )
	
	f.write( t )
	for j in np.arange(len(plotdict[colnames[0]])):
		f.write('<tr><td>{:.0f}</td>'.format(j))
		for col in colnames:
			if not (col in omitted_colnames):
				f.write('<td>{:.3f}</td>'.format(plotdict[col][j]))
		f.write('<tr>\n')
	f.write( '</table>\n' )
	
	f.write( '<BR><HR>\n' )
	f.write( '<table><tbody><tr><td><IMG SRC="photcomp.plot.png"/></td></tr></tbody></table>\n')
	f.write( '<BR><HR>\n' )

	f.write( '</BODY>\n' )
	f.write( '</HTML>\n' )
	f.close()