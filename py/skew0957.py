import sys
import os, os.path
import numpy as nu
import cPickle as pickle
from optparse import OptionParser
from varqso import VarQso, LCmodel
from fitQSO import QSOfilenames
from galpy.util import save_pickles
from plotFits import open_qsos
def skew0957(parser):
    (options,args)= parser.parse_args()
    if len(args) == 0:
        parser.print_help()
        return
    savefilename= args[0]
    #Read data
    if options.vanderriest:
        vA= VarQso('../data/0957-A.fits',band=options.band)
        vB= VarQso('../data/0957-B.fits',band=options.band)
    else:
        vA= VarQso('../data/L0957-A_%s.fits' % options.band,band=options.band)
        vB= VarQso('../data/L0957-B_%s.fits' % options.band,band=options.band)
    if not options.fitsfile is None and os.path.exists(options.fitsfile):
        fitsfile= open(options.fitsfile,'rb')
        paramsA= pickle.load(fitsfile)
        paramsB= pickle.load(fitsfile)
        paramsAB= pickle.load(fitsfile)
        fitsfile.close()
        vA.LCparams= paramsA
        vB.LCparams= paramsB
        vA.LCtype='powerlawSF'
        vB.LCtype='powerlawSF'
        vA.LCmean= 'const'
        vB.LCmean= 'const'
        vA.fitband= options.band
        vB.fitband= options.band
    else:
        #Fit for means
        print "Fitting SF for both images ..."
        vA.fit(options.band,mean='const')
        vB.fit(options.band,mean='const')
    #Load into single new VarQso
    newm= list(vA.m[options.band]-vA.LCparams['m'])
    newerrm= list(vA.err_m[options.band])
    newmjd= list(vA.mjd[options.band])
    newm.extend(list(vB.m[options.band]+nu.mean(vA.m[options.band])
                                 -nu.mean(vB.m[options.band])-vB.LCparams['m']))
    newerrm.extend(list(vB.err_m[options.band]))
    newmjd.extend(list(vB.mjd[options.band]-417./365.25))#shift lagged B
    v= VarQso(newmjd,newm,newerrm,band=options.band,medianize=False)
    if not options.fitsfile is None and os.path.exists(options.fitsfile):
        v.LCparams= paramsAB
        v.LCtype='powerlawSF'
        v.LCmean= 'zero'
        v.fitband= options.band
        v.LC= LCmodel(trainSet=v._build_trainset(options.band),
                      type=v.LCtype,mean=v.LCmean,
                      init_params=paramsAB)
    else:
        v.fit(options.band)
    if not options.fitsfile is None and not os.path.exists(options.fitsfile):
        save_pickles(options.fitsfile,vA.LCparams,vB.LCparams,v.LCparams)
    taus= nu.arange(1.,201.,1.)/365.25
    thisskew= v.skew(taus,options.band,duration=0.7)
    thisgaussskews= nu.zeros((options.nsamples,len(taus)))
    print "Calculating simulated skews ..."
    for ii in range(options.nsamples):
        #First re-sample
        redshift= 1.41
        o= v.resample(v.mjd[options.band],band=options.band,noconstraints=True,
                      wedge=options.wedge,
                      wedgerate=options.wedgerate*365.25/(1.+redshift),
                      wedgetau=(1.+redshift)) #1yr
        o.LCparams= v.LCparams
        o.LC= v.LC
        o.fitband= v.fitband
        o.LCtype= v.LCtype
        o.LCmean= v.LCmean
        if options.wedge:
            o.LCparams['gamma']= 1.
            o.LCparams['logA']= o.LCparams['logA']\
                +nu.log(0.05**v.LCparams['gamma']/0.05)
            o.LCmean= 'zero' #bc we remove the mean when resampling wedge
            #Set up LC with correct params
            o.LC= LCmodel(trainSet=o._build_trainset(options.band),
                          type=o.LCtype,mean=o.LCmean,
                          init_params=o.LCparams)
        thisgaussskews[ii,:]= o.skew(taus,options.band,duration=0.7)
    skews= {}
    gaussskews= {}
    skews['0957']= thisskew
    gaussskews['0957']= thisgaussskews
    save_pickles(savefilename,skews,gaussskews,
                 None,options.band,None,taus)
    return None

def get_options():
    usage = "usage: %prog [options] <savefilename>\n\nsavefilename= name of the file that the skews will be saved to"
    parser = OptionParser(usage=usage)
    parser.add_option("-b","--band",dest='band',default='r',
                      help="band(s) to sample")
    parser.add_option("-n","--nsamples",dest='nsamples',
                      default=100,type='int',
                      help="Number of samples to take")
    parser.add_option("-f","--fitsfile",dest='fitsfile',
                      default=None,
                      help="File that holds the best-fits")
    parser.add_option("--dtau",dest='dtau',
                      default=1,type='float',
                      help="lag spacing")
    parser.add_option("--taumax",dest='taumax',
                      default=40.,type='float',
                      help="lag spacing")
    parser.add_option("--wedge",action="store_true", dest="wedge",
                      default=False,
                      help="Use wedge model")
    parser.add_option("--wedgerate",dest='wedgerate',
                      default=0.1,type='float',
                      help="wedge rate (rest-frame; /days)")
    parser.add_option("--vanderriest",action="store_true", dest="vanderriest",
                      default=False,
                      help="Use Vanderriest data")
    return parser

if __name__ == '__main__':
    skew0957(get_options())
