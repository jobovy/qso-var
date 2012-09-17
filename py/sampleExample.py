import sys
import os, os.path
import cPickle as pickle
import scipy
import numpy as nu
from numpy import array
import galpy.util.bovy_plot as bovy_plot
from astrometry.util import pyfits_utils as fu
from gp.trainGP import trainGP
from gp.sampleGP import sampleGP
from gp.trainingSet import trainingSet
#from gp.OU_ARD import covarFunc
def singleBandSampleExample(datafilename='../data/SDSSJ203817.37+003029.8.fits',
                  band='g',nsamples=1000,
                  basefilename='SDSSJ203817.37+003029.8'):
    """
    NAME:
       singleBandSampleExample
    PURPOSE:
       Sample an example of a power-law structure function GP covariance 
       function to an SDSS quasar
    INPUT:
       datafilename
       band - band to fit
       nsamples - number of samples to use
       basefilename
    OUTPUT:
    HISTORY:
       2010-08-08 - Written - Bovy (NYU)
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
    mjd_r= nu.array(file.mjd_r)/365.25
    r= nu.array(file.r)
    err_r= nu.array(file.err_r)
    mjd_i= nu.array(file.mjd_i)/365.25
    i= nu.array(file.i)
    err_i= nu.array(file.err_i)

    mask= (mjd_g != 0)*(g < 20.6)*(g > 19.7)
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

    i= i[mask]
    i-= nu.mean(i)
    err_i= err_i[mask]
    mjd_i= mjd_i[mask]
    mjd_i-= nu.amin(mjd_i)
    meanErr_i= nu.mean(err_i)

    savefilename= basefilename+'.sav'
    if os.path.exists(savefilename):
        savefile= open(savefilename,'rb')
        gammas= pickle.load(savefile)
        if 'gr' in basefilename:
            logGammas= pickle.load(savefile)
            gammagrs= pickle.load(savefile)
            logGammagrs= pickle.load(savefile)
        else:
            logAs= pickle.load(savefile)
        out= pickle.load(savefile)
        savefile.close()
    else:
        if 'gri' in basefilename:
            raise NotImplementedError
            from powerlawSFgr import covarFunc
            params= {'logGamma': array([-7.79009776]), 
                     'logGammagr': array([-28.0487848]), 
                     'gamma': array([ 0.45918053]), 
                     'gammagr': array([ 0.21333858])}
            SF= covarFunc(**params)
            listx= [(t,'g') for t in mjd_g]
            listx.extend([(t,'r') for t in mjd_r])
            listx.extend([(t,'i') for t in mjd_i])
            listy= [m for m in g]
            listy.extend([m for m in r])
            listy.extend([m for m in i])
            listy= nu.array(listy)
            noise= [m for m in err_g]
            noise.extend([m for m in err_r])
            noise.extend([m for m in err_i])
            noise= nu.array(noise)
            trainSet= trainingSet(listx=listx,listy=listy,noise=noise)
        elif 'gr' in basefilename:
            raise NotImplementedError
            from powerlawSFgr import covarFunc
            params= {'logGamma': array([-7.79009776]), 
                     'logGammagr': array([-28.0487848]), 
                     'gamma': array([ 0.45918053]), 
                     'gammagr': array([ 0.21333858])}
            SF= covarFunc(**params)
            listx= [(t,'g') for t in mjd_g]
            listx.extend([(t,'r') for t in mjd_r])
            listy= [m for m in g]
            listy.extend([m for m in r])
            listy= nu.array(listy)
            noise= [m for m in err_g]
            noise.extend([m for m in err_r])
            noise= nu.array(noise)
            trainSet= trainingSet(listx=listx,listy=listy,noise=noise)
        else:
            from gp.powerlawSF import covarFunc
            trainSet= trainingSet(listx=mjd_g,listy=g,noise=err_g)
            params= [0.,0.]
            SF= covarFunc(gamma=.2,A=0.000524) #Power
            #SF= covarFunc(a=.0001,l=.001) #OU
    
        #Train
        print "Training ..."
        outcovarFunc= trainGP(trainSet,SF,useDerivs=False)
        
        print "Best-fit:"
        print outcovarFunc._dict
        
        print "Sampling ..."
        out= sampleGP(trainSet,outcovarFunc,nsamples=nsamples,
                      step=[0.2 for ii in range(len(params))])
        
        gammas= nu.array([o._dict['gamma'] for o in out]).reshape(len(out))
        if 'gr' in basefilename:
            gammagrs= nu.array([o._dict['gammagr'] for o in out]).reshape(len(out))
            logGammas= nu.array([o._dict['logGamma'] for o in out]).reshape(len(out))
            logGammagrs= nu.array([o._dict['logGammagr'] for o in out]).reshape(len(out))
        else:
            logAs= nu.array([o._dict['logA'] for o in out]).reshape(len(out))
            print nu.mean(logAs),nu.sqrt(nu.var(logAs))
        print nu.mean(gammas),nu.sqrt(nu.var(gammas))
        
        #Save
        savefile= open(savefilename,'wb')
        pickle.dump(gammas,savefile)
        if 'gr' in basefilename:
            pickle.dump(logGammas,savefile)
            pickle.dump(gammagrs,savefile)
            pickle.dump(logGammagrs,savefile)
        else:
            pickle.dump(logAs,savefile)
        pickle.dump(out,savefile)
        savefile.close()

    #if not 'gr' in basefilename:
    #    logGammas= logAs/gammas

    #Plot
    bovy_plot.bovy_print()
    bovy_plot.scatterplot(logAs/2.,gammas,'k,',
                          ylabel=r'$\gamma$',xlabel=r'\log A',
                          xrange=[-9.21/2.,0.],
                          yrange=[0.,1.25],
                          bins=50,onedhists=True)
    bovy_plot.bovy_end_print(basefilename+'_sample_2d.png')

if __name__ == '__main__':
    if len(sys.argv) > 2:
        singleBandSampleExample(covarType=sys.argv[2],
                                basefilename=sys.argv[1])
    elif len(sys.argv) > 1:
        singleBandSampleExample(basefilename=sys.argv[1])
    else:
        singleBandSampleExample()
