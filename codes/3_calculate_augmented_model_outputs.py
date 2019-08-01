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


# An example calling sequence to derive augmented model outputs based on original model outputs
# and the sampled, augmented parameter sets:
#
# python calculate_augmented_model_outputs.py \
#                                -i model_output_original.out \
#                                -a parameter_sets_augmented.out \
#                                -o model_output_augmented.out \
#                                -e 0.1 \
#                                -n 1000 \
#                                -m ['sobol']

"""
Derives the augmented model outputs y_MVA based on:
- the original model outputs y_Ori (saved in file specified with -i),
- the sampled augmented parameters z0, z1, and z2 (saved in file given with -a) as well as 
- the MVA parameter delta (specified with -e).

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

nsets       = 10                               # number of Sobol sequences
infile      = 'model_output_original.out'      # name of file used to save (scalar) original model outputs
augfile     = 'parameter_sets_augmented.out'   # name of file that contains sampled augmented parameters z0, z1, and z2
outfile     = 'model_output_augmented.out'     # name of file used to save (scalar) augmented model outputs
delta       = 0.1                              # MVA parameter delta
method      = ['sobol']                        # SA method that is going to be applied later:
#                                              # supported options:
#                                              #       'sobol'
#                                              #       ['pawn',Nf] where Nf is number of samples per approximated CDF

parser   = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                  description='''Derives the augmented model outputs y_MVA based on: (a) the original model outputs y_Ori (saved in file specified with -i), (b) the sampled augmented parameters z0, z1, and z2 (saved in file given with -a) as well as (c) the MVA parameter delta (specified with -d).''')
parser.add_argument('-i', '--infile', action='store',
                    default=infile, dest='infile', metavar='infile',
                    help="Name of file used to save (scalar) original model outputs (default: 'model_output_original.out').")
parser.add_argument('-a', '--augfile', action='store',
                    default=augfile, dest='augfile', metavar='augfile',
                    help="Name of file that contains sampled augmented parameters z0, z1, and z2 (default: 'parameter_sets_augmented.out').")
parser.add_argument('-o', '--outfile', action='store',
                    default=outfile, dest='outfile', metavar='outfile',
                    help="Name of file used to save (scalar) augmented model outputs (default: 'model_output_augmented.out').")
parser.add_argument('-e', '--delta', action='store',
                    default=delta, dest='delta', metavar='delta',
                    help="Parameter delta used for MVA estimates (default: 0.1).")
parser.add_argument('-n', '--nsets', action='store',
                    default=nsets, dest='nsets', metavar='nsets',
                    help='Number of sensitivity samples (default: nsets=10).')
parser.add_argument('-m', '--method', action='store',
                    default=None, dest='method', metavar='method',
                    help="SA method that is going to be applied later. Supported options are ['sobol'] and ['pawn',Nf] where Nf is number of samples per approximated CDF. (default: ['sobol']).")

args     = parser.parse_args()
infile   = args.infile
outfile  = args.outfile
delta    = np.float(args.delta)
nsets    = np.int(args.nsets)

# some arguments need some formatting
if args.method is not None:
    tmp = args.method

    # '[sobol]'     --> ['sobol']
    # '[pawn, 50]'  --> ['pawn',50]
    tmp = tmp.replace('[','').split('],')
    tmp = [ s.replace(']','').split(',') for s in tmp ][0]
    
    tmp[0] = str(tmp[0])
    if len(tmp) == 2:
        tmp[1] = np.int(tmp[1])
    method = tmp

del parser, args

# read parameter sets for augmented variables z0, z1, z2
ff = open(augfile, "r")
zz = ff.readlines()
ff.close()
zz = np.array([ map(float,izz.strip().split()) for izz in zz ])
z0 = zz[:,0]
z1 = zz[:,1]
z2 = zz[:,2]

# read original model outputs
ff = open(infile, "r")
y_ori = ff.readlines()
ff.close()
y_ori = np.array(map(float,y_ori))

npara = np.shape(y_ori)[0]/nsets - 2

# derive c used for Eq. 2 and derived in Eq. 3 in Mai & Tolson (2019)
if method[0] == 'sobol':
    
    model_ori_a = y_ori[0:nsets]
    model_ori_b = y_ori[nsets:2*nsets]
    
elif method[0] == 'pawn':

    print("ToDo: Sampling to run PAWN method")
    
else:
    print('method = ',method[0])
    raise ValueError('This method is not implemented yet! Only "sobol" and "pawn".')    

mu_f  = np.mean(np.append(model_ori_a, model_ori_b))
mu_f2 = np.mean(np.append(model_ori_a**2,model_ori_b**2))
var_f = np.var(np.append(model_ori_a,model_ori_b))
cc2   = (2.0*(mu_f2-mu_f**2)-2.0/3.0*delta**2*(var_f+mu_f**2)-var_f)/var_f
if (cc2 < -0.5):
    print("")
    print("c**2 = ", cc2, " but must be positive!!")    #### THIS SHOULD NOT HAPPEN!!!
    print("You likely need more samples or chosen delta is too large!")
    stop
elif (cc2 < 0.0):  # is between -0.1 and 0.0 --> sampling error
    cc = 0.0
else:
    cc = np.sqrt(cc2)

# write augmented model outputs
ff = open(outfile, "w")

if method[0] == 'sobol':
    
    # derive augmented model output y_MVA as described in Eq. 2 in Mai & Tolson (2019)
    # A sets
    y_MVA_a = 0.0*z0[0:nsets]       + z1[0:nsets]*model_ori_a       - z2[0:nsets]*model_ori_a       + cc*model_ori_a
    for iset in range(nsets):
        ff.write( str(y_MVA_a[iset])+'\n' )
    y_MVA_b = 0.0*z0[nsets:2*nsets] + z1[nsets:2*nsets]*model_ori_b - z2[nsets:2*nsets]*model_ori_b + cc*model_ori_b
    # B sets
    for iset in range(nsets):
        ff.write( str(y_MVA_b[iset])+'\n' )
    # Ci sets with original model variables
    for ipara in range(npara):
        model_ori_ci = y_ori[nsets*(2+ipara):nsets*(3+ipara)]
        y_MVA_ci = 0.0*z0[0:nsets] + z1[0:nsets]*model_ori_ci - z2[0:nsets]*model_ori_ci + cc*model_ori_ci
        for iset in range(nsets):
            ff.write( str(y_MVA_ci[iset])+'\n' )
    # Ci sets with augmented model variables
    y_MVA_ci_z0 = 0.0*z0[nsets:2*nsets] + z1[0:nsets]*model_ori_a       - z2[0:nsets]*model_ori_a       + cc*model_ori_a
    y_MVA_ci_z1 = 0.0*z0[0:nsets]       + z1[nsets:2*nsets]*model_ori_a - z2[0:nsets]*model_ori_a       + cc*model_ori_a
    y_MVA_ci_z2 = 0.0*z0[0:nsets]       + z1[0:nsets]*model_ori_a       - z2[nsets:2*nsets]*model_ori_a + cc*model_ori_a
    for iset in range(nsets):
        ff.write( str(y_MVA_ci_z0[iset])+'\n' )
    for iset in range(nsets):
        ff.write( str(y_MVA_ci_z1[iset])+'\n' )
    for iset in range(nsets):
        ff.write( str(y_MVA_ci_z2[iset])+'\n' )

elif method[0] == 'pawn':

    print("ToDo: Sampling to run PAWN method")
    
else:
    print('method = ',method[0])
    raise ValueError('This method is not implemented yet! Only "sobol" and "pawn".')

ff.close()
print("wrote:   '"+outfile+"'")


