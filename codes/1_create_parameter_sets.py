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


# An example calling sequence to create 1000 parameter sets with one uniform
# and one Gaussian distributed parameter to perform a PAWN sensitivity analysis with 50 samples for each CDF
# and a delta for MVA of 0.1 is:
#
# python create_parameter_sets.py -n 1000 \
#                                 -d [[uniform,0.1,0.4],[gaussian,0.0,2.0]] \
#                                 -e 0.1 \
#                                 -m [sobol] \
#                                 -o parameter_sets
#
# or a calling sequence for a Sobol' sensitivity analysis with three parameters that are all uniformly distributed
# between -pi and pi and a delta for MVA of 0.1 is:
#
# python create_parameter_sets.py -n 1000 \
#                                 -d [['uniform',-3.14159,3.14159],['uniform',-3.14159,3.14159],['uniform',-3.14159,3.14159]] \
#                                 -e 0.1 \
#                                 -m ['sobol'] \
#                                 -o parameter_sets

"""
Sample parameter sets using stratified sampling of Sobol' sequences. The parameters can be specified to be 
either uniform or Gaussian distributed. The parameter sets can be specified to be either used for a Sobol' 
or PAWN sensitivity analysis. The parameter sets are stored in a text file.

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

nsets       = 10                                             # number of Sobol sequences
dists       = [['uniform',0.1,0.4],['gaussian',0.0,2.0]]     # distribution and distribution parameters for each model parameter
delta       = 0.1                                            # MVA parameter delta
method      = ['sobol']                                      # SA method that is going to be applied later:
#                                                            # supported options:
#                                                            #       'sobol'
#                                                            #       ['pawn',Nf] where Nf is number of samples per approximated CDF
outfile     = 'parameter_sets'                               # name of file to store sampled parameter sets

parser   = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                  description='''Sample parameter sets using stratified sampling of Sobol' sequences. The parameters can be specified to be 
er uniform or Gaussian distributed. The parameter sets can be specified to be either used for a Sobol' 
AWN sensitivity analysis. The parameter sets are stored in a text file.''')
parser.add_argument('-n', '--nsets', action='store',
                    default=nsets, dest='nsets', metavar='nsets',
                    help='Number of sensitivity samples (default: nsets=10).')
parser.add_argument('-d', '--dists', action='store',
                    default=None, dest='dists', metavar='dists',
                    help="Distribution and distribution parameters for each model parameter (default: [['uniform',0.1,0.4],['gaussian',0.0,2.0]]).")
parser.add_argument('-e', '--delta', action='store',
                    default=delta, dest='delta', metavar='delta',
                    help="Parameter delta used for MVA estimates (default: 0.1).")
parser.add_argument('-m', '--method', action='store',
                    default=None, dest='method', metavar='method',
                    help="SA method that is going to be applied later. Supported options are ['sobol'] and ['pawn',Nf] where Nf is number of samples per approximated CDF. (default: ['sobol']).")
parser.add_argument('-o', '--outfile', action='store',
                    default=outfile, dest='outfile', metavar='outfile',
                    help="Baseame of file to store sampled parameter sets. Two files will be created <outfile>'_original.out' containing parametersets to run the model and <outfile>'_augmented.out' containing the sampled MVA variables (default: 'parameter_sets').")

args     = parser.parse_args()
nsets    = np.int(args.nsets)
outfile  = args.outfile
delta    = np.float(args.delta)

# some arguments need some formatting
if args.dists is not None:
    tmp = args.dists

    # '[[uniform,0.1,0.4],[gaussian,0.0,2.0]]' --> [['uniform',0.1,0.4],['gaussian',0.0,2.0]]
    tmp = tmp.replace('[','').split('],')
    tmp = [ s.replace(']','').split(',') for s in tmp]
    
    for ii,idist in enumerate(tmp):
        tmp[ii][0] = str(tmp[ii][0])
        tmp[ii][1] = np.float(tmp[ii][1])
        tmp[ii][2] = np.float(tmp[ii][2])
    dists = tmp

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

print('dists:  ',dists)
print('method: ',method)
print("")


npara = len(dists)
skip  = 30000

if method[0] == 'sobol':
    sobol_sets = sobol.i4_sobol_generate(2*(npara+3),nsets,skip)    # add three MVA parameters z0, z1, and z2

    # scale sobol sets to
    for ipara in range(npara):
        if dists[ipara][0] == 'uniform':
            sobol_sets[ipara,:]         = sobol_sets[ipara,:]         * (dists[ipara][2]-dists[ipara][1]) + dists[ipara][1]
            sobol_sets[ipara+npara+3,:] = sobol_sets[ipara+npara+3,:] * (dists[ipara][2]-dists[ipara][1]) + dists[ipara][1]
        elif dists[ipara][0] == 'gaussian':
            # N[0,1]; cut at 1% and 99% percentiles
            sobol_sets[ipara,:]         = 0.01 + 0.98 * sobol_sets[ipara,:]
            sobol_sets[ipara+npara+3,:] = 0.01 + 0.98 * sobol_sets[ipara+npara+3,:]
            sobol_sets[ipara,:]         = stats.norm.ppf(sobol_sets[ipara,:],         loc=dists[ipara][1], scale=dists[ipara][2])
            sobol_sets[ipara+npara+3,:] = stats.norm.ppf(sobol_sets[ipara+npara+3,:], loc=dists[ipara][1], scale=dists[ipara][2])
        else:
            print('dist = ',dists[ipara])
            raise ValueError('This distribution type is not implemented yet! Only "uniform" and "gaussian".')

    # ---------------------
    # sets of variables to run the model
    # ---------------------
    ff = open(outfile+'_original.out', "w")

    # A sets: first half of Sobol' sets  (without three MVA parameters: sobol_sets[0:npara,:])   
    a_sets = copy.deepcopy(np.transpose(sobol_sets[0:npara,:]))
    for iset in a_sets:
        ff.write(' '.join([ str(ii) for ii in iset])+'\n')
        
    # B sets: second half of Sobol' sets (without three MVA parameters: sobol_sets[npara+3:2*npara+3,:])
    b_sets = copy.deepcopy(np.transpose(sobol_sets[npara+3:2*npara+3,:]))
    for iset in b_sets:
        ff.write(' '.join([ str(ii) for ii in iset])+'\n')

    # Ci sets: everything from A except parameter i is from B
    for ipara in range(npara):
        c_sets          = copy.deepcopy(a_sets)
        c_sets[:,ipara] = b_sets[:,ipara]
        for iset in c_sets:
            ff.write(' '.join([ str(ii) for ii in iset])+'\n')

    ff.close()
    print("wrote:   '"+outfile+'_original.out'+"'")

    # ---------------------
    # sets of MVA variables
    # ---------------------
    ff = open(outfile+'_augmented.out', "w")

    # A sets: first half of Sobol' sets  (only three MVA parameters: sobol_sets[npara:npara+3,:]) 
    a_sets = copy.deepcopy(np.transpose(sobol_sets[npara:npara+3,:]))
    a_sets[:,1] = a_sets[:,1] * 2.0 * delta + (1.-delta)     # scale to range: [1-delta, 1+delta]
    a_sets[:,2] = a_sets[:,2] * 2.0 * delta + (1.-delta)     # scale to range: [1-delta, 1+delta]
    for iset in a_sets:
        ff.write(' '.join([ str(ii) for ii in iset])+'\n')
        
    # B sets: second half of Sobol' sets (without three MVA parameters: sobol_sets[2*npara+3:2*npara+6,:])
    b_sets = copy.deepcopy(np.transpose(sobol_sets[2*npara+3:2*npara+6,:]))
    b_sets[:,1] = b_sets[:,1] * 2.0 * delta + (1.-delta)     # scale to range: [1-delta, 1+delta]
    b_sets[:,2] = b_sets[:,2] * 2.0 * delta + (1.-delta)     # scale to range: [1-delta, 1+delta]
    for iset in b_sets:
        ff.write(' '.join([ str(ii) for ii in iset])+'\n')

    ff.close()
    print("wrote:   '"+outfile+'_augmented.out'+"'")

elif method[0] == 'pawn':

    print("ToDo: Sampling to run PAWN method")
    
else:
    print('method = ',method[0])
    raise ValueError('This method is not implemented yet! Only "sobol" and "pawn".')    
