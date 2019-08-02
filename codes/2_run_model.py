#!/usr/bin/env python
from __future__ import print_function

# Copyright 2019 Juliane Mai - juliane.mai(at)uwaterloo.ca
#
# License
# This file is part of the code library for "Model Variable Augmentation (MVA) 
# for Diagnostic Assessment of Sensitivity Analysis Results".
#
# The MVA code library is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The MVA code library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with The MVA code library.
# If not, see <https://github.com/julemai/MVA/blob/master/LICENSE>.
#
# If you use this method in a publication please cite:
#
#    Mai, J., & Tolson, B. A. ( 2019). 
#    Model Variable Augmentation (MVA) for diagnostic assessment of sensitivity analysis results. 
#    Water Resources Research, 55, 2631-2651.
#    https://doi.org/10.1029/2018WR023382


# An example calling sequence to derive model outputs for Ishigami-Homa function (or any other model):
#
# python 2_run_model.py -i parameter_sets_original.out \
#                       -o model_output_original.out

"""
Runs a model for a bunch of parameter sets and stores the (scalar) model outputs in a file.

History
-------
Written,  JM, Mar 2019
"""

# -------------------------------------------------------------------------
# Command line arguments - if script
#

# Comment|Uncomment - Begin
#if __name__ == '__main__':

# -----------------------
# add subolder scripts/lib to search path
# -----------------------
import sys
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path+'/lib')

import argparse
import numpy as np
import scipy.stats as stats
import copy
import sobol                  # in lib/

infile      = 'parameter_sets_original.out'      # name of file containing sampled parameter sets to run the model
outfile     = 'model_output_original.out'        # name of file used to save (scalar) model outputs

parser   = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                  description='''Runs a model for a bunch of parameter sets and stores the (scalar) model outputs in a file.''')
parser.add_argument('-i', '--infile', action='store',
                    default=infile, dest='infile', metavar='infile',
                    help="Name of file containing sampled parameter sets to run the model (default: 'parameter_sets_original.out').")
parser.add_argument('-o', '--outfile', action='store',
                    default=outfile, dest='outfile', metavar='outfile',
                    help="Name of file used to save (scalar) model outputs (default: 'model_output_original.out').")

args     = parser.parse_args()
infile   = args.infile
outfile  = args.outfile

del parser, args


def model_function(paraset):
    # function that takes parameter set and returns (scalar) model output
    # here: Ishigami-Homa function (Ishigami and Homma, [1990])
    #            f(x) = sin(p1) + a * sin(p2)**2 + b * p3**4 * sin(p1)
    #       with
    #            a=2.0 and b=1.0

    a = 2.0
    b = 1.0
    
    model = np.sin(paraset[0]) + a * np.sin(paraset[1])**2 + b * paraset[2]**4 * np.sin(paraset[0])
    
    return model

# read parameter sets
ff = open(infile, "r")
parasets = ff.readlines()
ff.close()

# write model outputs
ff = open(outfile, "w")

for paraset in parasets:
    paraset = map(float,paraset.strip().split())
    model = model_function(paraset)
    ff.write(str(model)+'\n')

ff.close()
print("wrote:   '"+outfile+"'")
        
