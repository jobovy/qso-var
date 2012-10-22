import sys
import os, os.path
import numpy as nu
import cPickle as pickle
from optparse import OptionParser
from varqso import VarQso, LCmodel
from fitQSO import QSOfilenames
from galpy.util import save_pickles
from plotFits import open_qsos
_DEBUG=True
_ERASESTR= "                                                                                "
def skewQSO(parser):
    (options,args)= parser.parse_args()
    if len(args) == 0:
        parser.print_help()
        return
    savefilename= args[0]
    if os.path.exists(savefilename):
        savefile= open(savefilename,'rb')
        skews= pickle.load(savefile)
        gaussskews= pickle.load(savefile)
        fittype= pickle.load(savefile)
        band= pickle.load(savefile)
        mean= pickle.load(savefile)
        taus= pickle.load(savefile)      
        savefile.close()
    else:
        skews= {}
        gaussskews= {}
        fittype= options.type
        mean= options.mean
        band= options.band
        taus= nu.arange(options.dtau,options.taumax,options.dtau)/365.25
    if os.path.exists(options.fitsfile):
        fitsfile= open(options.fitsfile,'rb')
        params= pickle.load(fitsfile)
        fitsfile.close()
    else:
        raise IOError("--fitsfile (or -f) has to be set to the file holding the best-fits")
    if options.star:
        dir= '../data/star/'
    elif options.nuvx:
        dir= '../data/nuvx/'
    elif options.nuvxall:
        dir= '../data/nuvx_all/'
    elif options.uvx:
        dir= '../data/uvx/'
    elif options.rrlyrae:
        dir= '../data/rrlyrae/'
    else:
        dir= '../data/s82qsos/'
    qsos= QSOfilenames(dir=dir)
    if not options.split is None:
        splitarr= nu.arange(len(qsos)) / int(nu.ceil(len(qsos)/float(options.split)))
        splitDict= {}
        for ii, qso in enumerate(qsos):
            key= os.path.basename(qso)
            splitDict[key]= splitarr[ii]
        print "Running bin %i ..." % options.rah
    savecount= 0
    count= len(skews)
    if not options.star:
        #Read master file for redshifts
        dataqsos= open_qsos()
        qsoDict= {}
        ii=0
        for qso in dataqsos:
            qsoDict[qso.oname.strip().replace(' ', '')+'.fit']= ii
            ii+= 1
    for qso in qsos:
        key= os.path.basename(qso)
        if skews.has_key(key):
            continue
        if not options.split is None:
            if splitDict[key] != options.rah:
                continue
        else:
            try:
                if int(key[5:7]) != options.rah and options.rah != -1:
                    continue
            except ValueError:
                if options.rah == -2 or options.rah == -1:
                    pass
                else:
                    print "Skipping ValueError "+key
                    continue
        sys.stdout.write('\r'+_ERASESTR+'\r')
        sys.stdout.flush()
        sys.stdout.write('\rWorking on %s: %s\r' % (str(count),key))
        sys.stdout.flush()
        v= VarQso(qso)
        if v.nepochs(band) < 20:
            #print "This object does not have enough epochs ..."
            continue
        #Set best-fit
        v.LCparams= params[key]
        v.LC= LCmodel(trainSet=v._build_trainset(band),type=fittype,mean=mean,
                      init_params=params[key])
        v.LCtype= fittype
        v.LCmean= mean
        v.fitband= band
        #Now compute skew and Gaussian samples
        try:
            thisskew= v.skew(taus,band)
        except RuntimeError:
            continue
        thisgaussskews= nu.zeros((options.nsamples,len(taus)))
        for ii in range(options.nsamples):
            #First re-sample
            if options.star:
                redshift= 0.
            else:
                redshift= dataqsos[qsoDict[key]].z
            try:
                o= v.resample(v.mjd[band],band=band,noconstraints=True,
                              wedge=options.wedge,
                              wedgerate=options.wedgerate*365.25/(1.+redshift),
                              wedgetau=(1.+redshift)) #1yr
            except nu.linalg.linalg.LinAlgError:
                if params[key]['gamma'] > 1.5:
                    continue #re-sampling fails because of bad gamma
                else:
                    print key, params[key]
            o.LCparams= v.LCparams
            o.LC= v.LC
            o.fitband= v.fitband
            o.LCtype= v.LCtype
            o.LCmean= v.LCmean
            if options.wedge:
                o.LCparams['gamma']= 1.
                o.LCparams['logA']= o.LCparams['logA']\
                    +nu.log(0.05**v.LCparams['gamma']/0.05)
            thisgaussskews[ii,:]= o.skew(taus,band)
        skews[key]= thisskew
        gaussskews[key]= thisgaussskews
        savecount+= 1
        if savecount == options.saveevery:
            sys.stdout.write('\r'+_ERASESTR+'\r')
            sys.stdout.flush()
            sys.stdout.write('\rSaving ...\r')
            sys.stdout.flush()
            save_pickles(savefilename,skews,gaussskews,fittype,band,mean,taus)
            savecount= 0
        count+= 1
    sys.stdout.write('\r'+_ERASESTR+'\r')
    sys.stdout.flush()
    save_pickles(savefilename,skews,gaussskews,fittype,band,mean,taus)
    print "All done"

def get_options():
    usage = "usage: %prog [options] <savefilename>\n\nsavefilename= name of the file that the skews will be saved to"
    parser = OptionParser(usage=usage)
    parser.add_option("-b","--band",dest='band',default='r',
                      help="band(s) to sample")
    parser.add_option("-t","--type",dest='type',default='powerlawSF',
                      help="Type of model to sample (powerlawSF, powerlawSFratios, or DRW)")
    parser.add_option("--mean",dest='mean',default='zero',
                      help="Type of mean to sample (zero, const)")
    parser.add_option("-f","--fitsfile",dest='fitsfile',
                      default=None,
                      help="File that holds the best-fits")
    parser.add_option("-n","--nsamples",dest='nsamples',
                      default=100,type='int',
                      help="Number of samples to take")
    parser.add_option("--saveevery",dest='saveevery',type='int',
                      default=10,
                      help="Save every --saveevery iterations")
    parser.add_option("--rah",dest='rah',type='int',
                      default=-1,
                      help="RA hour to consider (-1: all, -2: ValueError), if combined with split, bin to use")
    parser.add_option("--split",dest='split',type='int',
                      default=None,
                      help="split the sample in this number of bins")
    parser.add_option("--star",action="store_true", dest="star",
                      default=False,
                      help="Sample stars")
    parser.add_option("--nuvx",action="store_true", dest="nuvx",
                      default=False,
                      help="Sample nUVX sample")
    parser.add_option("--uvx",action="store_true", dest="uvx",
                      default=False,
                      help="Sample UVX sample")
    parser.add_option("--nuvxall",action="store_true", dest="nuvxall",
                      default=False,
                      help="Sample nUVX_all sample")
    parser.add_option("--rrlyrae",action="store_true", dest="rrlyrae",
                      default=False,
                      help="Sample RR Lyrae sample")
    parser.add_option("--resampled",action="store_true", dest="resampled",
                      default=False,
                      help="Objects are 'resampled': stored in sav-file")
    parser.add_option("-i","--infile",dest='infile',
                      default=None,
                      help="Input file if --resampled")
    parser.add_option("--dtau",dest='dtau',
                      default=3,type='float',
                      help="lag spacing")
    parser.add_option("--taumax",dest='taumax',
                      default=150.,type='float',
                      help="lag spacing")
    parser.add_option("--wedge",action="store_true", dest="wedge",
                      default=False,
                      help="Use wedge model")
    parser.add_option("--wedgerate",dest='wedgerate',
                      default=0.1,type='float',
                      help="wedge rate (rest-frame; /days)")
    return parser

if __name__ == '__main__':
    skewQSO(get_options())
