import sys
import os, os.path
import scipy
import numpy as nu
from astrometry.util import pyfits_utils as fu
from gp.trainGP import trainGP, read_train_out
from gp.trainingSet import trainingSet
from drwcross import covarFunc
def fitExampleDRW(datafilename='../data/SDSSJ203817.37+003029.8.fits'):
    """
    NAME:
       fitExampleDRW
    PURPOSE:
       Fit an example of a DRW structure function GP covariance function
       to an SDSS quasar in g and r
    INPUT:
       datafilename
    OUTPUT:
    HISTORY:
       2010-02-21 - Written - Bovy (NYU)
    """
    nu.random.seed(1)
    #Get data
    file= fu.table_fields(datafilename)

    mjd_g= nu.array(file.mjd_g)/365.25
    g= nu.array(file.g)
    err_g= nu.array(file.err_g)
    mjd_r= nu.array(file.mjd_r)/365.25
    r= nu.array(file.r)
    err_r= nu.array(file.err_r)

    mask= (mjd_g != 0)*(g < 20.6)*(g > 19.7)#Adjust for non-g
    g= g[mask]
    g-= nu.mean(g)
    err_g= err_g[mask]
    mjd_g= mjd_g[mask]
    mjd_g-= nu.amin(mjd_g)
    meanErr_g= nu.mean(err_g)

    r= r[mask]
    r-= nu.mean(r)
    err_r= err_r[mask]
    mjd_r= mjd_r[mask]
    mjd_r-= nu.amin(mjd_r)
    meanErr_r= nu.mean(err_r)

    listx= [(t,'g') for t in mjd_g]
    listx.extend([(t,'r') for t in mjd_r])
    listy= [m for m in g]
    listy.extend([m for m in r])
    listy= nu.array(listy)
    noise= [m for m in err_g]
    noise.extend([m for m in err_r])
    noise= nu.array(noise)
    trainSet= trainingSet(listx=listx,listy=listy,noise=noise)
    
    SF= covarFunc(S=1.,tau=.2,B=.1,C=.005,gammagr=.5,Gammagr=.01)

    #Train
    print "Training ..."
    outcovarFunc= trainGP(trainSet,SF,useDerivs=False)

    print "Best-fit:"
    print outcovarFunc._dict

    return None

if __name__ == '__main__':
    fitExampleDRW()
