import re
import os, os.path
import cPickle as pickle
import numpy as nu
import pyfits
from optparse import OptionParser
from varqso import VarQso
from fitQSO import QSOfilenames
def classQSO(parser):
    (options,args)= parser.parse_args()
    if len(args) == 0:
        parser.print_help()
        return
    if os.path.exists(args[0]):
        print filename+" exists"
        print "Remove this file before running ..."
        print "Returning ..."
        return None
    #Load fit params: Quasars
    if options.qsomodel == 'test':
        if len(options.band) == 1:
            qsoparams= [{'gamma':0.3,'logA':-2.}]
        else:#Multi-band
            qsoparams= [{'gamma':0.2,'logA':-2.,
                         'gammagr':0.0001,'logAgr':-2.}]
        qsoweights= [1.]
    elif options.qsomodel == 'zero':
        qsoparams= [{}]
        qsoweights= [1.]
    elif os.path.exists(options.qsomodel):
        qsofile= open(options.qsomodel,'rb')
        qsoparams= pickle.load(qsofile)
        qsoweights= nu.array(pickle.load(qsofile),dtype='float64')
        qsofile.close()
    else:
        print "Input to 'qsomodel' not recognized ..."
        print "Returning ..."
        return
    #Stars
    if options.starmodel == 'test':
        if len(options.band) == 1:
            starparams= [{'gamma':0.0001,'logA':-3.5}]
        else:#Multi-band
            starparams= [{'gamma':0.1,'logA':-3.,
                          'gammagr':0.0001,'logAgr':-2.}]
        starweights= [1.]
    elif options.starmodel == 'zero':
        starparams= [{}]
        starweights= [1.]
    elif os.path.exists(options.starmodel):
        starfile= open(options.starmodel,'rb')
        starparams= pickle.load(starfile)
        starweights= nu.array(pickle.load(starfile),dtype='float64')
        starfile.close()
    else:
        print "Input to 'starmodel' not recognized ..."
        print "Returning ..."
        return
    #RR Lyrae
    if options.rrlyraemodel == 'test':
        if len(options.band) == 1:
            rrlyraeparams= [{'gamma':0.0001,'logA':-2.}]
        else:#Multi-band
            rrlyraeparams= [{'gamma':0.1,'logA':-2.,
                          'gammagr':0.0001,'logAgr':-2.}]
        rrlyraeweights= [1.]
    elif options.rrlyraemodel == 'zero':
        rrlyraeparams= [{}]
        rrlyraeweights= [1.]
    elif os.path.exists(options.rrlyraemodel):
        rrlyraefile= open(options.rrlyraemodel,'rb')
        rrlyraeparams= pickle.load(rrlyraefile)
        rrlyraeweights= nu.array(pickle.load(rrlyraefile),dtype='float64')
        rrlyraefile.close()
    else:
        print "Input to 'rrlyraemodel' not recognized ..."
        print "Returning ..."
        return
    #normalize weights
    qsoweights/= nu.sum(qsoweights)
    starweights/= nu.sum(starweights)
    rrlyraeweights/= nu.sum(rrlyraeweights)
    #Load location of the data
    if options.resampled:
        if os.path.exists(options.sample):
            samplefile= open(options.sample,'rb')
            objs= pickle.load(samplefile)
            samplefile.close()
        else:
            print "'--resampled' is set, but --sample= filename does not exist ..."
            print "Returning ..."
            return None
    else:
        if options.sample == 'nuvx':
            dir= '../data/nuvx/'
        if options.sample == 'nuvxall':
            dir= '../data/nuvx_all/'
        if options.sample == 'uvx':
            dir= '../data/uvx/'
        objs= QSOfilenames(dir=dir)
    #Classify each source
    out= []
    allcount, count= 0, 0
    for obj in objs:
        allcount+= 1
        if options.resampled:
            key= obj[0]
        else:
            key= os.path.basename(obj)
        #if key != 'SDSSJ013306.18-004523.8.fit':
        #    continue
        print "Working on "+str(count)+"(%i/%i): " % (allcount,len(objs))+key
        if options.resampled:
            v= obj[1]
        else:
            v= VarQso(obj)
        if v.nepochs(options.band) < options.minepochs:
            print "This object does not have enough epochs ..."
            continue
        varout= VarClass()
        varout.key= key
        #quasar likelihoods
        qsolike= []
        for ii in range(len(qsoparams)):
            qsolike.append(v.loglike(band=options.band,type=options.type,
                                     params=qsoparams[ii])
                           +nu.log(qsoweights[ii]))
        qsolike= logsum(qsolike)
        varout.qsologlike= qsolike
        #star likelihoods
        starlike= []
        for ii in range(len(starparams)):
            starlike.append(v.loglike(band=options.band,type=options.type,
                                     params=starparams[ii])
                           +nu.log(starweights[ii]))
        starlike= logsum(starlike)
        varout.starloglike= starlike
        #RR Lyrae likelihoods
        rrlyraelike= []
        for ii in range(len(rrlyraeparams)):
            rrlyraelike.append(v.loglike(band=options.band,type=options.type,
                                         params=rrlyraeparams[ii])
                               +nu.log(rrlyraeweights[ii]))
        rrlyraelike= logsum(rrlyraelike)
        varout.rrlyraeloglike= rrlyraelike
        #print qsolike, starlike
        if qsolike > starlike and qsolike > rrlyraelike:
            print qsolike, starlike, rrlyraelike
        out.append(varout)
        count+= 1
        #if count > 500: break
    #Save
    for jj in range(len(out)):
        if out[jj].qsologlike > out[jj].starloglike and out[jj].qsologlike > out[jj].rrlyraeloglike:
            print out[jj].key
    saveClass(out,args[0])
    return None

def saveClass(out,filename):
    """
    NAME:
       saveClass
    PURPOSE:
       save the classifications
    INPUT:
       out - listof varClass objects key, qsologlike, starloglike)
       filename - name of the file that the output will be saved to
    OUTPUT:
       (none)
    HISTORY:
       2011-01-30 - Written - Bovy (NYU)
    """
    key= [re.split(r'.fit',o.key)[0] for o in out]
    lenkey= [len(k) for k in key]
    #Prepare columns
    cols= []
    colkey= pyfits.Column(name='key',format=str(max(lenkey))+'A',array=key)
    cols.append(colkey)
    colqso= pyfits.Column(name='qsologlike',format='E',
                          array=[o.qsologlike for o in out])
    cols.append(colqso)
    colstar= pyfits.Column(name='starloglike',format='E',
                           array=[o.starloglike for o in out])
    cols.append(colstar)
    colrrlyrae= pyfits.Column(name='rrlyraeloglike',format='E',
                              array=[o.rrlyraeloglike for o in out])
    cols.append(colrrlyrae)
    #Save
    columns= pyfits.ColDefs(cols)  
    tbhdu= pyfits.new_table(columns)
    tbhdu.writeto(filename)

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
    array= nu.array(array)
    c= nu.amax(array)
    return nu.log(nu.nansum(nu.exp(nu.add(array,-c))))+c

def get_options():
    usage = "usage: %prog [options] <savefilename>\n\nsavefilename= name of the file that the classifications will be saved to"
    parser = OptionParser(usage=usage)
    parser.add_option("-b","--band",dest='band',default='r',
                      help="band(s) to use")
    parser.add_option("-t","--type",dest='type',default='powerlawSF',
                      help="Type of model to use (powerlawSF or DRW)")
    parser.add_option("--minepochs",dest='minepochs',default=2,type='int',
                      help="Minimum number of epochs for an object to be considered")
    parser.add_option("--sample",dest='sample',default='nuvx',
                      help="sample to classify ('nuvx', 'uvx', or resampled-filename)")
    parser.add_option("--qsomodel",dest='qsomodel',default='test',
                      help="Model for the quasars ('test', or filename)")
    parser.add_option("--starmodel",dest='starmodel',default='test',
                      help="Model for the stars ('test', or filename)")
    parser.add_option("--rrlyraemodel",dest='rrlyraemodel',default='test',
                      help="Model for the rrlyrae ('test', or filename)")
    parser.add_option("--resampled",action="store_true", dest="resampled",
                      default=False,
                      help="Sample to classify is a resampled sample")
    return parser

if __name__ == '__main__':
    classQSO(get_options())
