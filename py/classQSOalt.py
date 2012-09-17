_MEANS={'logA':-2.,'gamma':0.2,'s':0.75}
#_STDS={'logA':.5,'gamma':.2,'s':1.}
_STDS={'logA':.000005,'gamma':.00002,'s':1.}
import re
import os, os.path
import copy
import cPickle as pickle
import numpy as nu
from scipy import spatial
import pyfits
from optparse import OptionParser
from varqso import VarQso
from fitQSO import QSOfilenames
from classQSO import logsum
def classQSOalt(parser):
    (options,args)= parser.parse_args()
    if len(args) == 0:
        parser.print_help()
        return
    if options.nearest < 10000:
        print "Sparse nearest neighbor sampling not yet supported ..."
        print "Returning ..."
        return None
    if os.path.exists(args[0]):
        print filename+" exists"
        print "Remove this file before running ..."
        print "Returning ..."
        return None
    #Load fit params: Quasars
    if os.path.exists(options.qsomodel):
        qsofile= open(options.qsomodel,'rb')
        qsoparams= pickle.load(qsofile)
        qsofile.close()
    else:
        print "Input to 'qsomodel' not recognized ..."
        print "Returning ..."
        return
    #Stars
    if os.path.exists(options.starmodel):
        starfile= open(options.starmodel,'rb')
        starparams= pickle.load(starfile)
        starfile.close()
    else:
        print "Input to 'starmodel' not recognized ..."
        print "Returning ..."
        return
    #RR Lyrae
    if os.path.exists(options.rrlyraemodel):
        rrlyraefile= open(options.rrlyraemodel,'rb')
        rrlyraeparams= pickle.load(rrlyraefile)
        rrlyraefile.close()
    else:
        print "Input to 'rrlyraemodel' not recognized ..."
        print "Returning ..."
        return
    #Set up KD trees for nearest neighbors search
    print "Setting up KD tree ..."
    qsokd= spatial.KDTree(_prepare_KD(qsoparams,poploglike=True))
    starkd= spatial.KDTree(_prepare_KD(starparams,poploglike=True))
    rrlyraekd= spatial.KDTree(_prepare_KD(rrlyraeparams,poploglike=True))
    #Load location of the data
    if os.path.exists(options.input):
        savefile= open(options.input,'rb')
        fitparams= pickle.load(savefile)
        type= pickle.load(savefile)
        band= pickle.load(savefile)
        try:
            mean= pickle.load(savefile)
        except EOFError:
            pass
        finally:
            savefile.close()
    else:
        print options.input+" does not exist ..."
        print "Returning ..."
        return None
    print "Setting up objects to classify ..."
    objs= fitparams.keys()
    nobjs= len(objs)
    nparams= len(fitparams.values()[0].keys())-1
    #Prepare objs for KD Tree
    objsKDin= nu.zeros((nobjs,nparams))
    for ii in range(nobjs):
        for jj in range(nparams):
            key= fitparams.values()[ii].keys()[jj]
            if key == 'loglike': continue
            objsKDin[ii,jj]= (fitparams.values()[ii][key]-_MEANS[key])/_STDS[key]
    #Classify each source
    out= []
    count= 0
    print "Classifying sources ..."
    for ii in range(nobjs):
        key= objs[ii]
        print "Working on "+str(count)+"(%i/%i): " % (ii+1,len(objs))+key
        varout= VarClass()
        varout.key= key
        #quasar likelihoods
        #Find nearest neighbors
        (qsonn,mems)= qsokd.query(objsKDin[ii,:],k=nu.amin([options.nearest,qsokd.data.shape[0]]))
        qsolike= []
        for jj in range(nu.amin([options.nearest,qsokd.data.shape[0]])):
            qsolike.append(-0.5*qsonn[jj]**2.)
        qsolike= logsum(qsolike)-nu.log(nu.amin([options.nearest,qsokd.data.shape[0]]))\
            +fitparams.values()[ii]['loglike']
        varout.qsologlike= qsolike
        #star likelihoods
        #Find nearest neighbors
        (starnn,mems)= starkd.query(objsKDin[ii,:],k=nu.amin([options.nearest,starkd.data.shape[0]]))
        starlike= []
        for jj in range(nu.amin([options.nearest,starkd.data.shape[0]])):
            starlike.append(-0.5*starnn[jj]**2.)
        starlike= logsum(starlike)-nu.log(nu.amin([options.nearest,starkd.data.shape[0]]))\
            +fitparams.values()[ii]['loglike']
        varout.starloglike= starlike
        #RR Lyrae likelihoods
        (rrlyraenn,mems)= rrlyraekd.query(objsKDin[ii,:],k=nu.amin([options.nearest,rrlyraekd.data.shape[0]]))
        rrlyraelike= []
        for jj in range(nu.amin([options.nearest,rrlyraekd.data.shape[0]])):
            rrlyraelike.append(-0.5*rrlyraenn[jj]**2.)
        rrlyraelike= logsum(rrlyraelike)-nu.log(nu.amin([options.nearest,rrlyraekd.data.shape[0]]))\
            +fitparams.values()[ii]['loglike']
        varout.rrlyraeloglike= rrlyraelike
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

def _prepare_KD(params,poploglike=True):
    """Internal function to whiten input parameters"""
    if poploglike:
        nparams= len(params.values()[0].keys())-1
    else:
        nparams= len(params.values()[0].keys())
    nobjs= len(params)
    kdin= nu.zeros((nobjs,nparams))
    for ii in range(nobjs):
        for jj in range(nparams):
            key= params.values()[ii].keys()[jj]
            if poploglike and key == 'loglike': continue  
            kdin[ii,jj]= (params.values()[ii][key]-_MEANS[key])/_STDS[key]
    return kdin

def get_options():
    usage = "usage: %prog [options] <savefilename>\n\nsavefilename= name of the file that the classifications will be saved to"
    parser = OptionParser(usage=usage)
    parser.add_option("-b","--band",dest='band',default='r',
                      help="band(s) to use")
    parser.add_option("-t","--type",dest='type',default='powerlawSF',
                      help="Type of model to use (powerlawSF or DRW)")
    parser.add_option("-n","--nearest",dest='nearest',default=10000,
                      type='int',
                      help="Number of nearest neighbors to use")
    parser.add_option("--minepochs",dest='minepochs',default=2,type='int',
                      help="Minimum number of epochs for an object to be considered")
    parser.add_option("-i","--input",dest='input',default=None,
                      help="sample to classify, fits")
    parser.add_option("--qsomodel",dest='qsomodel',default='test',
                      help="Model for the quasars ('test', or filename)")
    parser.add_option("--starmodel",dest='starmodel',default='test',
                      help="Model for the stars ('test', or filename)")
    parser.add_option("--rrlyraemodel",dest='rrlyraemodel',default='test',
                      help="Model for the rrlyrae ('test', or filename)")
    #parser.add_option("--resampled",action="store_true", dest="resampled",
    #                  default=False,
    #                  help="Sample to classify is a resampled sample")
    return parser

if __name__ == '__main__':
    classQSOalt(get_options())
