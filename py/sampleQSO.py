import os, os.path
import numpy as nu
import cPickle as pickle
from optparse import OptionParser
import signal
from varqso import VarQso, LCmodel
from fitQSO import save_pickles, QSOfilenames
_DEBUG=True
#For time-out reasons
def handler(signum, frame):
    raise Exception("Sampling timed out")
def sampleQSO(parser):
    (options,args)= parser.parse_args()
    if len(args) == 0:
        parser.print_help()
        return
    if len(args) == 2:
        othersavefilename= args[1]
        othersavefile= open(othersavefilename,'rb')
        othersamples= pickle.load(othersavefile)
        othersavefile.close()
    else:
        othersamples= {}
    savefilename= args[0]
    if os.path.exists(savefilename):
        savefile= open(savefilename,'rb')
        samples= pickle.load(savefile)
        type= pickle.load(savefile)
        band= pickle.load(savefile)
        mean= pickle.load(savefile)
        savefile.close()
    else:
        samples= {}
        type= options.type
        mean= options.mean
        band= options.band
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
    if options.resampled:
        raise NotImplementedError("resampled not implemented yet")
        if os.path.exists(options.infile):
            samplefile= open(options.infile,'rb')
            qsos= pickle.load(samplefile)
            samplefile.close()
        else:
            print "'--resampled' is set, but -i filename does not exist ..."
            print "Returning ..."
            return None
    else:
        qsos= QSOfilenames(dir=dir)
    #Register time-out handler
    signal.signal(signal.SIGALRM, handler)
    savecount= 0
    count= len(samples)
    for qso in qsos:
        if options.resampled:
            key= qso[0]
        else:
            key= os.path.basename(qso)
        if samples.has_key(key) or othersamples.has_key(key):
            continue
        try:
            if int(key[5:7]) != options.rah and options.rah != -1:
                continue
        except ValueError:
            if options.rah == -2 or options.rah == -1:
                pass
            else:
                print "Skipping ValueError "+key
                continue
        print "Working on "+str(count)+": "+key
        if options.resampled:
            v= qso[1]
        else:
            v= VarQso(qso)
        if options.lownepochs:
            if v.nepochs(band) >= 20:
                print "This object has too many epochs ..."
                continue
            elif v.nepochs(band) < 3:
                print "This object has < 3 epochs ..."
                continue
        elif not options.lownepochs and v.nepochs(band) < 20:
            print "This object does not have enough epochs ..."
            continue
        #Set best-fit
        v.LCparams= params[key]
        v.LC= LCmodel(trainSet=v._build_trainset(band),type=type,mean=mean,
                      init_params=params[key])
        v.LCtype= type
        v.LCmean= mean
        v.fitband= band
        #Now sample
        signal.alarm(options.timeout)
        try:
            v.sampleGP(nsamples=options.nsamples,metropolis=options.metropolis,
                       markovpy=options.markovpy,
                       burnin=int(nu.floor(0.2*options.nsamples)))
        except Exception, exc:
            if str(exc) == "Sampling timed out":
                print exc
                continue
            else:
                raise
        signal.alarm(0)
        samples[key]= v.get_sampleGP()
        if _DEBUG:
            _print_diagnostics(samples[key])
            #print samples[key][options.nsamples]
        savecount+= 1
        if savecount == options.saveevery:
            print "Saving ..."
            save_pickles(samples,type,band,mean,savefilename)
            savecount= 0
        count+= 1
    save_pickles(samples,type,band,mean,savefilename)
    print "All done"

def _print_diagnostics(sample):
    """Internal function to print some diagnostics"""
    keys= sample[0].keys()
    for key in sample[0].keys():
        xs= []
        for ii in range(len(sample)):
            xs.append(sample[ii][key][0])
        xs= nu.array(xs)
        print key, nu.mean(xs), nu.std(xs)

def get_options():
    usage = "usage: %prog [options] <savefilename>\n\nsavefilename= name of the file that the samples will be saved to"
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
                      help="RA hour to consider (-1: all, -2: ValueError)")
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
    parser.add_option("--timeout",dest='timeout',
                      default=0,type='int',
                      help="Time out for individual object sampling (in sec)")
    parser.add_option("--metropolis",action="store_true", dest="metropolis",
                      default=False,
                      help="Use Metropolis sampling")
    parser.add_option("--markovpy",action="store_true", dest="markovpy",
                      default=False,
                      help="Use markovpy sampling")
    parser.add_option("--lownepochs",action="store_true", dest="lownepochs",
                      default=False,
                      help="Fit sources with a small number of epochs rather than a large number of epochs (< 20 rather than > 20)")
    return parser

if __name__ == '__main__':
    sampleQSO(get_options())
