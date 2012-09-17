import varqso
import os, os.path
import cPickle as pickle
from optparse import OptionParser
import numpy as nu
from galpy.util import bovy_plot
def plotGPSamples(parser):
    (options,args)= parser.parse_args()
    if len(args) == 0:
        parser.print_help()
        return
    #Restore fit
    if os.path.exists(args[0]):
        fitsfile= open(args[0],'rb')
        params= pickle.load(fitsfile)
        type= pickle.load(fitsfile)
        band= pickle.load(fitsfile)
        try:
            mean= pickle.load(fitsfile)
        except EOFError:
            mean= 'zero'
        fitsfile.close()
    else:
        raise IOError("You need to specify the file that holds the fits")
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
    #Load object
    v= varqso.VarQso(os.path.join(dir,options.key))
    v.LCparams= params[options.key]
    v.LC= varqso.LCmodel(trainSet=v._build_trainset(band),type=type,mean=mean)
    v.LCtype= type
    v.LCmean= mean
    v.fitband= band
    if os.path.exists(options.savefilename):
        print "Reusing samples from previous sampling"
        savefile= open(options.savefilename,'rb')
        samples= pickle.load(savefile)
        savefile.close()
        v.set_sampleGP(samples)
    else:
        v.sampleGP(nsamples=options.nsamples)
        samples= v.get_sampleGP()
        #Save
        savefile= open(options.savefilename,'wb')
        pickle.dump(samples,savefile)
        savefile.close()
    #Now plot
    bovy_plot.bovy_print()
    v.plot_sampleGP(d1=options.d1,d2=options.d2,whitemax=options.whitemax)
    bovy_plot.bovy_end_print(options.plotfilename)
    return None
    
def get_options():
    usage = "usage: %prog [options] <savefilename>\n\nsavefilename= name of the file that holds the fits"
    parser = OptionParser(usage=usage)
    parser.add_option("-n","--nsamples",dest='nsamples',default=1000,
                      type='int',
                      help="Number of samples to plot")
    parser.add_option("--key",dest='key',default=None,
                      help="Key of the object to plot")
    parser.add_option("--d1",dest='d1',default='logA',
                      help="parameter to plot on x")
    parser.add_option("--d2",dest='d2',default='gamma',
                      help="parameter to plot on y")
    parser.add_option("-o",dest='plotfilename',default=None,
                      help="Name for the file that will hold the plot")
    parser.add_option("-s",dest='savefilename',default=None,
                      help="Name for the file that will save the samples")
    parser.add_option("--xmin",dest='xmin',default=-9.21/2.,
                      type='float',
                      help="min of xrange")
    parser.add_option("--xmax",dest='xmax',default=0.,
                      type='float',
                      help="max of xrange")
    parser.add_option("--ymin",dest='ymin',default=0.,
                      type='float',
                      help="min of yrange")
    parser.add_option("--ymax",dest='ymax',default=1.25,
                      type='float',
                      help="max of yrange")
    parser.add_option("--star",action="store_true", dest="star",
                      default=False,
                      help="Sample star sample")
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
    parser.add_option("--whitemax",action="store_true", dest="whitemax",
                      default=False,
                      help="Use a white cross to indicate the maximum")
    return parser

if __name__ == '__main__':
    plotGPSamples(get_options())
    
