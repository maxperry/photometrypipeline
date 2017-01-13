#!/usr/bin/env python
'''
Command-line tool providing an interface to the 
 get_SEDs.py library.

Usage:
 python zeropoint.py <input_filename> <passband> <output_filename>
'''

from get_SEDs import *
from sys import argv

if __name__ == '__main__':
    in_file  = argv[1]
    passband = argv[2]
    if len(argv) > 3:
        out_file = argv[3]
    else:
        out_file = None
    
    zp, mad = zeropoint( in_file, passband, output_file=out_file )
    print zp, mad