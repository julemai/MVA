#!/usr/bin/env python
from   __future__ import print_function
import numpy as np
from   statsmodels.distributions.empirical_distribution import ECDF
# statsmodels seems to have an issue with scipy under Python 3.x:
# ticket:
#     https://github.com/statsmodels/statsmodels/issues/5747
# try:
#     pip install statsmodels==0.10.0rc2 --pre

def pawn_index( uncond, cond,
                alpha=0.05,
                pawn_stat='median'):
    """
        Estimation of parameter sensitivity using the PAWN method introduced by
              Pianosi, F. & Wagener T., 2015
              A simple and efficient method for global sensitivity analysis based on cumulative distribution functions.
              Environmental Modelling & Software 67, 1-11.


        Definition
        ----------
        def pawn_index( uncond, cond,
                        alpha=0.05,
                        pawn_stat='median'):

        Input
        ----------
        uncond     unconditional sets of model output
                   dim_1 = nu
        cond       conditional sets of model output
                   dim_1 = nn
                   dim_2 = npara
                   dim_3 = nc


        Optional Input
        --------------
        alpha           significance level of KS test
                        default: 0.05
        pawn_stat       statistic used for PAWN index; e.g. 'max', 'mean', 'median'
                        default: 'median'


        Output
        ------
        T_stats[0:npara]        PAWN sensitivity index per parameter (close to 1.0 very sensitive, close to 0.0 not sensitive)
        informative[0:npara]    certainty of parameter being sensitive (close to 1.0 very certain, close to 0.0 uncertain)


        Restrictions
        ------------
        1. PAWN statistics are limited to 'max', 'mean', 'median' at the moment


        References
        ----------
        Pianosi, F. & Wagener T., 2015
              A simple and efficient method for global sensitivity analysis based on cumulative distribution functions.
              Environmental Modelling & Software 67, 1-11.


        Examples
        --------

        >>> import numpy as np
        >>> import sobol
        >>> import copy
        >>> from autostring import astr

        >>> # number of unconditional parameter sets
        >>> nu = 100
        >>> # number of parameters
        >>> npara = 3
        >>> # number of conditional sets 
        >>> nc = 100
        >>> # number of replicates
        >>> nn = 50
        >>> # statistic used for PAWN index
        >>> pawn_stat = 'median'
        >>> # significance level of KS test
        >>> alpha = 0.05
        >>>
        >>> # sample base
        >>> # --> for unconditional exactly these sets are used
        >>> # --> for conditional ith column is replaced with constant value
        >>> sobol_sets = sobol.i4_sobol_generate(npara,nu,40000)    
        >>> fixed_vals = sobol.i4_sobol_generate(npara,nn,40000+nu) 
        >>> 
        >>> # for Ishigami-Homma sets need to be Uniform[-Pi,Pi]
        >>> sobol_sets[0:npara] = -np.pi + 2.* np.pi * sobol_sets[0:npara]
        >>> fixed_vals[0:npara] = -np.pi + 2.* np.pi * fixed_vals[0:npara]
        >>> 
        >>> # Ishigami-Homma has two constants in the model
        >>> a = 2.0
        >>> b = 1.0
        >>> 
        >>> # model outputs of unconditional sets
        >>> uncond = np.sin(sobol_sets[0,:])  + a*(np.sin(sobol_sets[1,:]))**2  + b*sobol_sets[2,:]**4  * np.sin(sobol_sets[0,:])
        >>> 
        >>> # model outputs of conditional sets
        >>> cond      = [ [ [] for ipara in range(npara) ] for irepl in range(nn) ]  
        >>> for ipara in range(npara):
        ...     for irepl in range(nn):
        ...         sobol_sets_cond          = copy.deepcopy(sobol_sets)
        ...         sobol_sets_cond[ipara,:] = fixed_vals[ipara,irepl]
        ...         # cond[irepl][ipara] = functions.ishigami_homma( sobol_sets_cond[0:npara,:], 2.0, 1.0)
        ...         cond[irepl][ipara] = np.sin(sobol_sets_cond[0,:])  + a*(np.sin(sobol_sets_cond[1,:]))**2  + b*sobol_sets_cond[2,:]**4  * np.sin(sobol_sets_cond[0,:])
        >>>
        >>> T_stat, influential = pawn_index(uncond, cond, pawn_stat=pawn_stat, alpha=alpha)
        >>>
        >>> print('PAWN_index :: T_stat       =',astr(T_stat,3,pp=True))
        PAWN_index :: T_stat       = ['0.460' '0.110' '0.270']
        >>> print('certainty  :: influential  =',astr(influential,3,pp=True))
        certainty  :: influential  = ['1.000' '1.000' '1.000']

        License
        -------
        This file is part of the JAMS Python package.

        It is NOT released under the GNU Lesser General Public License, yet.

        If you use this routine, please contact Juliane Mai.

        Copyright 2017 Juliane Mai


        History
        -------
        Written,  JM, Dec 2017
    """

    uncond = np.array(uncond)
    cond   = np.array(cond)
    
    nn    = np.shape(cond)[0]      # number replicates
    npara = np.shape(cond)[1]      # number parameters
    nc    = np.shape(cond)[2]      # number conditional sets
    nu    = np.shape(uncond)[0]    # number unconditional sets
    
    # estimated CDF of unconditional sets
    ecdf_uncond = ECDF(uncond)

    # estimated CDF of conditional sets
    ecdf_cond = [ [ ECDF(cond[irepl][ipara]) for ipara in range(npara) ] for irepl in range(nn) ]

    # get CDF values
    lb  = min(np.min(uncond),np.min(cond))
    ub  = max(np.max(uncond),np.max(cond))
    lb  = min(np.percentile(uncond, 0.5),np.percentile(cond, 0.5))
    ub  = max(np.percentile(uncond,99.5),np.percentile(cond,99.5))
    lb1 = lb - (ub-lb)*0.15
    ub1 = ub + (ub-lb)*0.15
    lb  = lb1
    ub  = ub1
    x_grid = np.arange(lb,ub+(ub-lb)/4000.,(ub-lb)/4000.)

    uncond_cdf = ecdf_uncond(x_grid)

    T_stats     = []   # PAWN sensitivity index
    informative = []   # gives certainty of parameter being informative (1.0=informative, 0.0=non-informative)
    cond_cdf    = [ [ [] for ipara in range(npara) ] for irepl in range(nn) ]
    for ipara in range(npara):
        KS_stat              = []
        count_influential    = 0.0
        count_noninformative = 0.0
        for irepl in range(nn):

            cond_cdf[irepl][ipara] = ecdf_cond[irepl][ipara](x_grid)
                
            # Eq. 10 in Zadeh et al. (2017), EMS
            #           "Comparison of variance-based and moment-independent global
            #           sensitivity approaches by application to the SWAT model"
            stat = np.max(np.abs(cond_cdf[irepl][ipara]-uncond_cdf))
            # if ipara==1:
                # print("stat = ",stat)
                # print("cond_cdf[",irepl,"][",ipara,"] = ",cond_cdf[irepl][ipara])
                # print("cond[irepl][ipara] = ", cond[irepl][ipara])
                # print("x_grid             = ", x_grid)
            KS_stat.append(stat)

            # null-hypothesis: two samples come from common distribution
            # rejected null hypothesis means that parameter is influential
            # reject null: p-value is equal or smaller than significance level alpha=0.05
            # found on:
            #       https://en.m.wikipedia.org/wiki/Kolmogorov-Smirnov_test
            # and:
            #       Massey (1951) and Marsaglia et al. (2003)
            ks_crit = np.sqrt(-0.5*np.log(alpha/2.0)) * np.sqrt((nu + nc) / (nu * nc))
            if stat > ks_crit:
                count_influential    += 1.0
            else:
                count_noninformative += 1.0

        if (pawn_stat == 'max'):
            T_stat = np.max(KS_stat)
        elif (pawn_stat == 'median'):
            T_stat = np.median(KS_stat)
        elif (pawn_stat == 'mean'):
            T_stat = np.mean(KS_stat)
        else:
            raise ValueError("Statistic not implemented yet!")

        informative.append(count_influential*1./nn)
        T_stats.append(T_stat)

    T_stats     = np.array(T_stats)
    informative = np.array(informative)    

    return T_stats, informative

if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
