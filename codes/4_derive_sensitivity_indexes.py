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


# An example calling sequence to calculate the sensitivity indexes for the previously derived model outputs. Technically this script works also with the original (non-augmented) model outputs. It will then derive the sensitivity indexes for the original model parameters only:
#
# python 4_derive_sensitivity_indexes.py \
#                  -i model_output_augmented.out \
#                  -o sensitivity_indexes_augmented.out \
#                  -n 1000 \
#                  -m ['sobol']

"""
Calculate the sensitivity indexes for the previously derived model outputs stored in the file specified by -i. Technically this script works also with the original (non-augmented) model outputs. It will then derive the sensitivity indexes for the original model parameters only.

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
import sobol_index              # in lib/
import pawn_index               # in lib/
from autostring import astr     # in lib/

nsets       = 10                                      # number of Sobol sequences
infile      = 'model_output_augmented.out'            # name of file used to save (scalar) original model outputs
outfile     = 'sensitivity_indexes_augmented.out'     # name of file used to save (scalar) augmented model outputs
method      = ['sobol']                               # SA method that is going to be applied later:
#                                                     # supported options:
#                                                     #       'sobol'
#                                                     #       ['pawn',Nf,stat,alpha] where Nf is number of conditioning values in PAWN method
#                                                     #                   Nf = parameter 'n' in Pianosi & Wagener (2015)
#                                                     #                   stat = statistic used in PAWN
#                                                     #                   alpha = confidence level of Kolmogorov-Smirnov test

parser   = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                  description='''Calculate the sensitivity indexes for the previously derived model outputs stored in the file specified by -i. Technically this script works also with the original (non-augmented) model outputs. It will then derive the sensitivity indexes for the original model parameters only.''')
parser.add_argument('-i', '--infile', action='store',
                    default=infile, dest='infile', metavar='infile',
                    help="Name of file used to save (scalar) model outputs (default: 'model_output_augmented.out').")
parser.add_argument('-o', '--outfile', action='store',
                    default=outfile, dest='outfile', metavar='outfile',
                    help="Name of file where sensitivity index estimates will be stored (default: 'sensitivity_indexes_augmented.out').")
parser.add_argument('-n', '--nsets', action='store',
                    default=nsets, dest='nsets', metavar='nsets',
                    help='Number of sensitivity samples (default: nsets=10).')
parser.add_argument('-m', '--method', action='store',
                    default=None, dest='method', metavar='method',
                    help="SA method that is applied. Supported options are ['sobol'] and ['pawn',Nf,stat,alpha] where Nf is number of conditioning values used for PAWN method (parameter 'n' in Pianosi & Wagener, 2015), stat is the statistic used (e.g., mean, median, or max), and alpha is the confidence level for the Kolmogorov-Smirnov test. (default: ['sobol']).")

args     = parser.parse_args()
infile   = args.infile
outfile  = args.outfile
nsets    = int(args.nsets)

# some arguments need some formatting
if args.method is not None:
    tmp = args.method

    # '[sobol]'     --> ['sobol']
    # '[pawn, 50]'  --> ['pawn',50]
    tmp = tmp.replace('[','').split('],')
    tmp = [ s.replace(']','').split(',') for s in tmp ][0]

    tmp[0] = str(tmp[0])
    if tmp[0] == 'pawn':
        tmp[1] = int(tmp[1])       # number of conditioning values
        tmp[2] = str(tmp[2])          # statistic: 'max', 'mean', or 'median'
        tmp[3] = float(tmp[3])     # alpha
    method = tmp

del parser, args

def array_to_string(array_to_write,dim):

    if dim == 0:
        string = str(array_to_write)+'\n'
    elif dim == 1:
        string = ' '.join(astr(array_to_write,prec=8))+'\n'
    else:
        raise ValueError('3_calculate_augmented_model_outputs.py: Only scalar and 1-dimensional model outputs are supported yet!')

    return string


# read original model outputs
ff = open(infile, "r")
y_MVA = ff.readlines()
ff.close()
y_MVA = np.array([ list(map(float,ii.strip().split())) for ii in y_MVA ])

# derive c used for Eq. 2 and derived in Eq. 3 in Mai & Tolson (2019)
if method[0] == 'sobol':

    npara = int( np.shape(y_MVA)[0]/nsets - 2 )
    ntime = int( np.shape(y_MVA)[1] )

    if ntime == 1:
        model_a = y_MVA[0:nsets,0]
        model_b = y_MVA[nsets:2*nsets,0]
        model_c = np.reshape(y_MVA[2*nsets:],[npara,nsets])
    else:
        model_a = y_MVA[0:nsets]
        model_b = y_MVA[nsets:2*nsets]
        model_c = np.reshape(y_MVA[2*nsets:],[npara,nsets,ntime])

        # time axis is expected to be first dimension
        model_a   = np.swapaxes(model_a,0,1)      # [nsets,ntime] --> [ntime,nsets]
        model_b   = np.swapaxes(model_b,0,1)      # [nsets,ntime] --> [ntime,nsets]
        model_c   = np.swapaxes(model_c,0,2)      # [npara,nsets,ntime] --> [ntime,nsets,npara]
        model_c   = np.swapaxes(model_c,1,2)      # [ntime,nsets,npara] --> [ntime,npara,nsets]

elif method[0] == 'pawn':

    nrepl      = method[1]   # number of conditioning values used in PAWN method (parameter "n" in Pianosi & Wagener, 2015)
    pawn_stat  = method[2]
    alpha      = method[3]
    npara      = int( (np.shape(y_MVA)[0]-nsets) / (nrepl*nsets) )
    ntime      = int(  np.shape(y_MVA)[1] )

    uncond = y_MVA[0:nsets]
    cond   = np.reshape(y_MVA[nsets:nsets+nrepl*nsets*(npara+3)],[npara,nrepl,nsets,ntime])   # shape is    [npara,nrepl,nsets,ntime]
    #                                                                                         # needs to be [nrepl,npara,nsets,ntime]
    cond   = np.swapaxes(cond,0,1)

else:
    print('method = ',method[0])
    raise ValueError('This method is not implemented yet! Only "sobol" and "pawn".')

# derive sensitivity index estimates with specified method
ff = open(outfile, "w")

if method[0] == 'sobol':

    si_MVA, sti_MVA = sobol_index.sobol_index(
                ya=model_a,
                yb=model_b,
                yc=model_c,
                si=True, sti=True,
                method='Mai1999')

    print("")
    print("si  = ",list(map(str,si_MVA)))
    print("sti = ",list(map(str,sti_MVA)))
    print("")

    ff.write("# Si(1:"+str(ntime)+")    STi(1:"+str(ntime)+") \n")
    for ipara in range(npara):
        if ntime == 1:
            ff.write(array_to_string(np.append(si_MVA[ipara],sti_MVA[ipara]),1))
        else:
            ff.write(array_to_string(np.append(si_MVA[:,ipara],sti_MVA[:,ipara]),1))

elif method[0] == 'pawn':

    pawn        = [[] for itime in range(ntime)]
    influential = [[] for itime in range(ntime)]
    for itime in range(ntime):
        pawn[itime], influential[itime] = pawn_index.pawn_index(
            uncond[:,itime],
            cond[:,:,:,itime],
            pawn_stat=pawn_stat,
            alpha=alpha)
    pawn        = np.array(pawn)
    influential = np.array(influential)

    print("")
    print("pawn index  = ",list(map(str,pawn)))
    print("influential = ",list(map(str,influential)))
    print("")

    ff.write("# PAWN(1:"+str(ntime)+")    Influential(1:"+str(ntime)+") \n")
    for ipara in range(npara):
        ff.write(array_to_string(np.append(pawn[:,ipara],influential[:,ipara]),1))

else:
    print('method = ',method[0])
    raise ValueError('This method is not implemented yet! Only "sobol" and "pawn".')

ff.close()
print("wrote:   '"+outfile+"'")
