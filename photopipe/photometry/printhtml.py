"""
NAME:
	printhtml
PURPOSE:
	Create photcomp.plot.png comparing magnitude and errors as well as create HTML page that has
	all of the data easily displayed for up to 9 filters
OUTPUT:
	photcomp.plot.png - shows magnitude vs. error for each filter
	abmag-energy.plot.png - shows flux vs. magnitude for each filter
	photom.html   - html page showing information about sources
DEPENDENCIES:
	Images and files (finalmags.txt) created from plotphotom.py

Translated from printratirhtml.pro by John Capone (jicapone@astro.umd.edu).
Modified 7/31/2014 by Vicki Toy (vtoy@astro.umd.edu)
"""

import numpy as np
import pylab as pl
import photprocesslibrary as pplib

def plotMag(filters, plotdict, colors):
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

def plotFlux(filters, plotdict, colors):
	pl.figure()
	pl.xlim([0,22000])
	pl.ylim([15,22])

	#For each filter, plot energy vs. abmag for abmag-energy.plot.png
	counter = 0
	for filter in filters:
		pl.plot(plotdict[filter+'flux'], plotdict[filter+'mag'], marker='o',linestyle='None', label=filter, color=colors[counter])
		counter = counter + 1
	
	pl.xlabel('cnts/sec (FLUX_APER)')
	pl.ylabel('AB Magnitude')
	pl.legend(loc='lower right')
	pl.savefig('abmag-energy.plot.png', bbox_inches='tight')
	pl.clf()	

def printhtml(filters, colnames):
	#Reads in final magnitudes
	colgrab = np.loadtxt('./finalmags.txt', unpack=True)
	
	#Store each column into dictionary based on colnames and create header for HTML header from colnames
	plotdict = {}
	t = '<tr><th>#</th>'
	for i in np.arange(len(colnames)):
		plotdict[colnames[i]] = colgrab[i,:]		
		t = t + '<th>' + colnames[i]+'</th>'
	
	t = t + '<tr>\n' 

	colors = ['black', 'purple','blue','aqua','green','orange','red','yellow', 'magenta']
	plotMag(filters, plotdict, colors)
	plotFlux(filters, plotdict, colors)

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
			f.write('<td>{:.3f}</td>'.format(plotdict[col][j]))
		f.write('<tr>\n')
	f.write( '</table>\n' )
	
	f.write( '<BR><HR>\n' )
	f.write( '<table><tbody><tr><td><IMG SRC="photcomp.plot.png"/></td><td><IMG SRC="abmag-energy.plot.png"/></td></tr></tbody></table>\n')
	f.write( '<BR><HR>\n' )

	f.write( '</BODY>\n' )
	f.write( '</HTML>\n' )
	f.close()