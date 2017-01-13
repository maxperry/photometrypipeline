import numpy

def median(l):
   a = numpy.array(l)
   return numpy.median(a)

####################################
def stdev(l):
   a = numpy.array(l)
   return numpy.std(a)
   
####################################
"""
NAME:
	most
PURPOSE:
	Finds the value in the list with the most amount of values that are within
	val-vmin and val+vmax.  An approximate mode calculation
INPUT:
	list - array of values you want to search
OPTIONAL KEYWORDS:
	vmin - sets search minimum to be value-vmin (default is 1)
	vmax - sets search maximum to be value+vmax (default is 1)
EXAMPLE:
	most(numpy.array([1.2, 4.3, 3.1, 5.1, 3.9, 7.1]), vmin=1, vmax=1)
		---> 4.3
"""  
def most(list, vmin=1, vmax=1):
	counter = numpy.zeros(len(list))
	
	#Finds number of elements that are between current element +/- vmax/vmin
	for i in range(0, len(list)):
		counter[i] =((list[i]+vmax >= list) & (list[i]-vmin <= list)).sum()
	
        #Returns the element that has the most values within vmin and vmax range
        #If no element with most, returns median
        if len(set(counter)) == 1:
            return numpy.median(list)
        else:
            return list[counter.argmax()]
	
####################################
def rasex2deg(rastr):
    rastr = str(rastr).strip()
    ra=rastr.split(':')
    if len(ra) == 1: return float(rastr)
    return 15*(float(ra[0])+float(ra[1])/60.0+float(ra[2])/3600.0)
    
####################################
def decsex2deg(decstr):
    decstr = str(decstr).strip()
    dec=decstr.split(':')
    if len(dec) == 1: return float(decstr)
    sign=1
    if (decstr[0] == '-'): sign=-1
    return sign*(abs(float(dec[0]))+float(dec[1])/60.0+float(dec[2])/3600.0)

####################################
#Compare objects using magnitude.
def magcomp(obj1, obj2): #useful for sorting
    return (obj1.mag > obj2.mag) - (obj1.mag < obj2.mag)

########################################
def unique(inlist):
    lis = inlist[:] #make a copy
    lis.sort()
    llen = len(lis)
    i = 0
    while i < llen-1:
        if lis[i+1] == lis[i]: 
            del lis[i+1]
            llen = llen - 1
        else:
            i = i + 1
    return lis