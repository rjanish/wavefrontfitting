#! /usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import subprocess

main_file = 'abr_fit_image.cpp'
implementations = [filename for filename in os.listdir(os.getcwd()) 
					if ((os.path.splitext(filename)[1] == ".cpp" or 
					   os.path.splitext(filename)[1] == ".c") and
					   filename !=  main_file)]
compiler = ["g++", "-Wall"] 
cfitsio_hdir = "/home/rjanish/cfitsio/include" # location of cfitsio library
cfitsio_ldir = "/home/rjanish/cfitsio/lib" # location of cfitsio library
header_dirs = ["-I", cfitsio_hdir] # add location of cfitsio library to list of directories in which to search for header files
link = ["-L%s" % cfitsio_ldir, "-lcfitsio"] # link to cfitsio library
if len(sys.argv) > 1:
	if sys.argv[1] == "-db":
		compiler.append("-ggdb")
	else:
		sys.exit('Invalid command line argument')
cmd = compiler + implementations + [main_file, "-o", main_file[:-4]] + link + header_dirs
compiler = subprocess.Popen(cmd, shell=False,
							stdout=subprocess.PIPE,
    					    stderr=subprocess.PIPE)
stdout, stderr = compiler.communicate()
stderr = stderr.split("\n")
for line in stderr:
	useless_messages = line.count('warning') + \
					   line.count('In function') + \
					   line.count('At global scope')
	if useless_messages < 1:
		print line
