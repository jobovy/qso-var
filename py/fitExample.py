import sys
import os, os.path
import scipy
import numpy as nu
from astrometry.util import pyfits_utils as fu
from gp.trainGP import trainGP, read_train_out
from gp.trainingSet import trainingSet
#from gp.powerlawSF import covarFunc
from gp.OU_ARD import covarFunc
def fitExample(N=None,datafilename='../data/SDSSJ203817.37+003029.8.fits',
               band='g'):
    """
    NAME:
       fitExample
    PURPOSE:
       Fit an example of a power-law structure function GP covariance function
       to an SDSS quasar
    INPUT:
       N - use a random subsample of N points
       datafilename
       band - band to fit
    OUTPUT:
    HISTORY:
       2010-02-21 - Written - Bovy (NYU)
    """
    nu.random.seed(1)
    #Get data
    file= fu.table_fields(datafilename)
    if band == 'u':
        mjd_g= nu.array(file.mjd_u)/365.25
        g= nu.array(file.u)
        err_g= nu.array(file.err_u)
    elif band == 'g':
        mjd_g= nu.array(file.mjd_g)/365.25
        g= nu.array(file.g)
        err_g= nu.array(file.err_g)
    elif band == 'r':
        mjd_g= nu.array(file.mjd_r)/365.25
        g= nu.array(file.r)
        err_g= nu.array(file.err_r)
    elif band == 'i':
        mjd_g= nu.array(file.mjd_i)/365.25
        g= nu.array(file.i)
        err_g= nu.array(file.err_i)
    elif band == 'z':
        mjd_g= nu.array(file.mjd_z)/365.25
        g= nu.array(file.z)
        err_g= nu.array(file.err_z)

    mask= (mjd_g != 0)*(g < 20.6)*(g > 19.7)
    g= g[mask]
    g-= nu.mean(g)
    err_g= err_g[mask]
    mjd_g= mjd_g[mask]
    mjd_g-= nu.amin(mjd_g)
    meanErr_g= nu.mean(err_g)

    trainSet= trainingSet(listx=mjd_g,listy=g,noise=err_g)
    SF= covarFunc(a=.0001,l=.001)
    #SF= covarFunc(gamma=.2,A=0.000524)
    #SF= covarFunc(a=0.18,l=0.005)
    
    #Train
    print "Training ..."
    outcovarFunc= trainGP(trainSet,SF,useDerivs=False)

    print "Best-fit:"
    print outcovarFunc._dict

    return None

if __name__ == '__main__':
    if len(sys.argv) > 1:
        fitExample(band=sys.argv[1])
    else:
        fitExample()
