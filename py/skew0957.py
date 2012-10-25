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
    #Fit for means
    vA.fit(options.band,mean='const')
    vB.fit(options.band,mean='const')
    #Load into single new VarQso
    newm= list(vA.m[options.band]-vA.Lcparams['m'])
    newerrm= list(vA.err_m[options.band])
    newmjd= list(vA.mjd[options.band])
    newm.extend(newm.extend(list(vB.m['g']+nu.mean(vA.m[options.band])
                                 -nu.mean(vB.m[options.band])-vB.LCparams['m'])))
    newerrm.extend(list(vB.err_m[options.band]))
    newmjd.extend(list(vB.mjd[options.band]-417./365.25))#shift lagged B
    v= varqso.VarQso(newmjd,newm,newerrm,band='g',medianize=False)
    v.fit(options.band)
    taus= nu.arange(1.,201.,1.)/365.25
    thisskew= v.skew(taus,options.band)
    thisgaussskews= nu.zeros((len(taus)))
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
            thisgaussskews[ii,:]= o.skew(taus,options.band)
    save_pickles(savefilename,skews,gaussskews,None,options.band,None,taus)
    return None

def get_options():
    usage = "usage: %prog [options] <savefilename>\n\nsavefilename= name of the file that the skews will be saved to"
    parser = OptionParser(usage=usage)
    parser.add_option("-b","--band",dest='band',default='r',
                      help="band(s) to sample")
    parser.add_option("-n","--nsamples",dest='nsamples',
                      default=100,type='int',
                      help="Number of samples to take")
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
