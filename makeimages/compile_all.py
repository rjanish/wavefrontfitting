#! /usr/bin/python
# -*- coding: utf-8 -*-

import os
from glob import glob
import subprocess

base_codes = glob("nr*.c")
for driver in glob("*.cpp"):
	print driver
	target = driver[:-4]
	cmd = ["g++", "-Wall"] + base_codes + [driver, "-o", target]
	compiler = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE,
        					    stderr=subprocess.PIPE)
	stdout, stderr = compiler.communicate()
	stderr = stderr.split("\n")
	common_errors = ["defined but not used", "unused variable", "In function ‘void newt(double*, int, int*, int*, void (*)(int, double*, double*))’", "nrutil.h: At global scope"]
	for line in stderr:
		if line == "":
			continue
		fake_error = False 
		for error in common_errors:
			fake_error = fake_error or (error in line) 
		if not fake_error:
			print line
