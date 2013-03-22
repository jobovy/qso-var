import sys
import re
import os, os.path
import cPickle as pickle
import numpy
from numpy import linalg
from scipy import maxentropy
import fitsio
from optparse import OptionParser
from skewQSO import _ERASESTR
_SQRTTWOPI= -0.5*numpy.log(2.*numpy.pi)
def classQSO(parser):
    (options,args)= parser.parse_args()
    if len(args) == 0:
        parser.print_help()
        return
    if os.path.exists(options.outfile):
        print options.outfile+" exists"
        print "Remove this file before running ..."
        print "Returning ..."
        return None
    #Load fit params: Quasars
    if os.path.exists(options.qsomodel):
        qsofile= open(options.qsomodel,'rb')
        try:
            xamp_qso= pickle.load(qsofile)
            xmean_qso= pickle.load(qsofile)
            xcovar_qso= pickle.load(qsofile)
        finally:
            qsofile.close()
    else:
        print "Input to 'qsomodel' not recognized ..."
        print "Returning ..."
        return
    #Stars
    if os.path.exists(options.starmodel):
        starfile= open(options.starmodel,'rb')
        try:
            xamp_star= pickle.load(starfile)
            xmean_star= pickle.load(starfile)
            xcovar_star= pickle.load(starfile)
        finally:
            starfile.close()
    else:
        print "Input to 'starmodel' not recognized ..."
        print "Returning ..."
        return
    ##RR Lyrae
    if os.path.exists(options.rrlyraemodel):
        rrlyraefile= open(options.rrlyraemodel,'rb')
        try:
            xamp_rrlyrae= pickle.load(rrlyraefile)
            xmean_rrlyrae= pickle.load(rrlyraefile)
            xcovar_rrlyrae= pickle.load(rrlyraefile)
        finally:
            rrlyraefile.close()
    else:
        print "Input to 'rrlyraemodel' not recognized ..."
        print "Returning ..."
        return
    #Restore samples
    savefilename= args[0]
    print "Reading data ..."
    if os.path.exists(savefilename):
        savefile= open(savefilename,'rb')
        samples= pickle.load(savefile)
        type= pickle.load(savefile)
        band= pickle.load(savefile)
        mean= pickle.load(savefile)
        savefile.close()
    else:
        print "Input file does not exist ..."
        print "Returning ..."
        return
    #Load the overall data, to later match back to ra and dec
    if 'nuvx' in args[0].lower():
        sources= fitsio.read('../data/nUVX_woname.fit')
    elif 'uvx' in args[0].lower():
        sources= fitsio.read('../data/uvx_woname.fit')
    sourcesDict= {}
    for ii in range(len(sources)):
        sourcesDict[sources[ii]['ONAME'].strip().replace(' ', '')+'.fit']= ii
    #Classify each source
    ndata= len(samples)
    print ndata
    logpxagamma_qso= numpy.zeros(ndata)
    logpxagamma_star= numpy.zeros(ndata)
    logpxagamma_rrlyrae= numpy.zeros(ndata)
    ras= numpy.zeros(ndata)
    decs= numpy.zeros(ndata)
    for ii, key in enumerate(samples.keys()):
        if ii > 20: break
        sys.stdout.write('\r'+_ERASESTR+'\r')
        sys.stdout.flush()
        sys.stdout.write('\rWorking on %i / %i\r' % (ii+1,ndata))
        sys.stdout.flush()
        if type == 'powerlawSF':
           #Stack as A,g,Ac,gc
            loggammas= []
            logAs= []
            for sample in samples[key]:
                loggammas.append(numpy.log(sample['gamma'][0]))
                logAs.append(sample['logA'][0]) #RITABAN
            loggammas= numpy.array(loggammas)
            logAs= numpy.array(logAs)
            weights= -loggammas #the 1/gamma to get a flat prior in log gamma, but expressed as log(1/gamma)
            weights-= maxentropy.logsumexp(weights) #sum weights = 1
            #Stack the data
            thisydata= numpy.reshape(loggammas,
                                     (len(loggammas),1))
            thisydata2= numpy.reshape(logAs,(len(logAs),1))
            thisydata=numpy.column_stack( [ thisydata, thisydata2 ] )
            #Evaluate quasar/star/RR lyrae distributions
            logpxagamma_qso[ii]= maxentropy.logsumexp(weights+_eval_sumgaussians(thisydata,
                                                                 xamp_qso,
                                                                 xmean_qso,
                                                                 xcovar_qso))
            logpxagamma_star[ii]= maxentropy.logsumexp(weights+_eval_sumgaussians(thisydata,
                                                                 xamp_star,
                                                                 xmean_star,
                                                                 xcovar_star))
            logpxagamma_rrlyrae[ii]= maxentropy.logsumexp(weights+_eval_sumgaussians(thisydata,
                                                                 xamp_rrlyrae,
                                                                 xmean_rrlyrae,
                                                                 xcovar_rrlyrae))
            #Find RA and Dec
            try:
                ratmp= sources[sourcesDict[key]]['RA']
                dectmp= sources[sourcesDict[key]]['DEC']
            except KeyError:
                print "Failed to match for RA and Dec ..."
                continue
            else:
                ras[ii]= ratmp
                decs[ii]= dectmp
    sys.stdout.write('\r'+_ERASESTR+'\r')
    sys.stdout.flush()
    #Save
    saveClass(logpxagamma_qso,
              logpxagamma_star,
              logpxagamma_rrlyrae,
              ras,decs,
              options.outfile)
    return None

def _eval_sumgaussians(x,xamp,xmean,xcovar):
    """x array [ndata,ndim], return log"""
    ndata= x.shape[0]
    da= x.shape[1]
    out= numpy.zeros(ndata)
    ngauss= len(xamp)
    loglike= numpy.zeros(ngauss)
    for ii in range(ndata):
        for kk in range(ngauss):
            if xamp[kk] == 0.:
                loglike[kk]= numpy.finfo(numpy.dtype(numpy.float64)).min
                continue
            tinv= linalg.inv(xcovar[kk,:,:])
            delta= x[ii,:]-xmean[kk,:]
            loglike[kk]= numpy.log(xamp[kk])+0.5*numpy.log(linalg.det(tinv))\
                -0.5*numpy.dot(delta,numpy.dot(tinv,delta))+\
                da*_SQRTTWOPI
        out[ii]= maxentropy.logsumexp(loglike)
    return out

def saveClass(q,s,r,ra,dec,filename):
    """
    NAME:
       saveClass
    PURPOSE:
       save the classifications
    INPUT:
       q, s, r, ra, dec
       filename - name of the file that the output will be saved to
    OUTPUT:
       (none)
    HISTORY:
       2011-01-30 - Written - Bovy (NYU)
    """
    #Create recarray
    ndata= len(q)
    out= numpy.recarray((ndata,),
                        dtype=[('ra','f8'),
                               ('dec','f8'),
                               ('logpx_qso','f8'),
                               ('logpx_star','f8'),
                               ('logpx_rrlyrae','f8')])
    out.logpx_qso= q
    out.logpx_star= s
    out.logpx_rrlyrae= r
    out.ra= ra
    out.dec= dec
    #Now write to fits
    fitsio.write(filename,out,clobber=True)
    return None

class VarClass:
    """empty class to hold results"""
    pass

def logsum(array):
    """
    NAME:
       logsum
    PURPOSE:
       calculate the logarithm of the sum of an array of numbers,
       given as a set of logs
    INPUT:
       array - logarithms of the numbers to be summed
    OUTPUT:
       logarithm of the sum of the exp of the numbers in array
    REVISION HISTORY:
       2009-09-29 -Written - Bovy (NYU)
    """

    #For now Press' log-sum-exp because I am too lazy to implement 
    #my own algorithm for this
    array= numpy.array(array)
    c= numpy.amax(array)
    return numpy.log(numpy.nansum(numpy.exp(numpy.add(array,-c))))+c

def get_options():
    usage = "usage: %prog [options] <savefilename>\n\nsavefilename= name of the file that holds samples"
    parser = OptionParser(usage=usage)
    parser.add_option("-o","--outfile",dest='outfile',default=None,
                      help="Name of the output fits file")
    parser.add_option("--qsomodel",dest='qsomodel',default='test',
                      help="Model for the quasars ('test', or filename)")
    parser.add_option("--starmodel",dest='starmodel',default='test',
                      help="Model for the stars ('test', or filename)")
    parser.add_option("--rrlyraemodel",dest='rrlyraemodel',default='test',
                      help="Model for the rrlyrae ('test', or filename)")
    return parser

if __name__ == '__main__':
    classQSO(get_options())
