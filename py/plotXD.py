import os, os.path
import cPickle as pickle
from optparse import OptionParser
import numpy as nu
import galpy.util.bovy_plot as bovy_plot
import xdtarget
def plotXD(parser):
    (options,args)= parser.parse_args()
    if len(args) == 0:
        parser.print_help()
        return
    if os.path.exists(args[0]):
        savefile= open(args[0],'rb')  
        xamp= pickle.load(savefile)
        xmean= pickle.load(savefile)
        xcovar= pickle.load(savefile)
        savefile.close()
    else:
        print args[0]+" does not exist ..."
        print "Returning ..."
        return
    #Load XD object in xdtarget
    xdt= xdtarget.xdtarget(amp=xamp,mean=xmean,covar=xcovar)
    out= xdt.sample(nsample=options.nsamples)
    #Prepare for plotting
    if options.expd1: xs= nu.exp(out[:,options.d1])
    elif not options.divided1 is None: xs= out[:,options.d1]/options.divided1
    else: xs= out[:,options.d1]
    if options.expd2: ys= nu.exp(out[:,options.d2])
    elif not options.divided2 is None: ys= out[:,options.d2]/options.divided2
    else: ys= out[:,options.d2]
    if options.type == 'DRW':
        #plot logA, logA = 0
        if options.d1 == 0 and options.d2 == 1: 
            #Convert to logA
            xs= (nu.log(2.)+xs+nu.log(1.-nu.exp(-1./nu.exp(ys))))/2.
        elif options.d1 ==1 and options.d2 == 0:
            #Convert to logA
            ys= (nu.log(2.)+ys+nu.log(1.-nu.exp(-1./nu.exp(xs))))/2.
        else:
            print "d1 and d2 have to be 0 or 1 (and not the same!) ..."
            print "Returning ..."
            return
    bovy_plot.bovy_print()
    bovy_plot.scatterplot(xs,ys,'k,',onedhists=True,
                          xrange=[options.xmin,options.xmax],
                          yrange=[options.ymin,options.ymax],
                          xlabel=options.xlabel,
                          ylabel=options.ylabel)
    bovy_plot.bovy_end_print(options.plotfilename)


def get_options():
    usage = "usage: %prog [options] <savefilename>\n\nsavefilename= name of the file that holds the XD fit"
    parser = OptionParser(usage=usage)
    parser.add_option("-o",dest='plotfilename',default=None,
                      help="Name for the file that will hold the plot")
    parser.add_option("-t","--type",dest='type',default='powerlawSF',
                      help="Type of model that we are considering")
    parser.add_option("-n","--nsamples",dest='nsamples',default=10000,
                      type='int',
                      help="Number of samples to plot")
    parser.add_option("--d1",dest='d1',default=0,
                      type='int',
                      help="index of first dimension")
    parser.add_option("--d2",dest='d2',default=1,
                      type='int',
                      help="index of second dimension")
    parser.add_option("--expd1",action="store_true",
                      dest="expd1",
                      default=False,
                      help="Plot exp(d1)")
    parser.add_option("--expd2",action="store_true",
                      dest="expd2",
                      default=False,
                      help="Plot exp(d2)")
    parser.add_option("--divided1",type='float',
                      dest="divided1",
                      default=None,
                      help="Divide d1 by this number (not compatible with --expd1)")
    parser.add_option("--divided2",type='float',
                      dest="divided2",
                      default=None,
                      help="Divide d2 by this number (not compatible with --expd2)")
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
    parser.add_option("--xlabel",dest='xlabel',default=r'$\log A$',
                      help="xlabel")
    parser.add_option("--ylabel",dest='ylabel',default=r'$\gamma$',
                      help="ylabel")
    return parser


if __name__ == '__main__':
    plotXD(get_options())
