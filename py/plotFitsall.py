import os, os.path
import cPickle as pickle
from optparse import OptionParser
import numpy as nu
import scipy as sc
import matplotlib
from galpy.util import bovy_plot
import xdtarget
def plotFitsall(parser):
    nu.random.seed(1)
    (options,args)= parser.parse_args()
    if len(args) == 0:
        parser.print_help()
        return
    if os.path.exists(args[0]):
        savefile= open(args[0],'rb')  
        params= pickle.load(savefile)
        savefile.close()
    else:
        print args[0]+" does not exist ..."
        print "Returning ..."
        return
    if os.path.exists(args[1]):
        savefile= open(args[1],'rb')  
        starparams= pickle.load(savefile)
        savefile.close()
    else:
        print args[1]+" does not exist ..."
        print "Returning ..."
        return
    if os.path.exists(args[2]):
        savefile= open(args[2],'rb')  
        rrparams= pickle.load(savefile)
        savefile.close()
    else:
        print args[2]+" does not exist ..."
        print "Returning ..."
        return
    #Prepare for plotting
    xs= nu.array([p[options.d1] for p in params.values()]).reshape(len(params))
    ys= nu.array([p[options.d2] for p in params.values()]).reshape(len(params))
    if options.expd1: xs= nu.exp(xs)
    elif not options.divided1 is None: xs/= options.divided1
    if options.expd2: ys= nu.exp(ys)
    elif not options.divided2 is None: ys/= options.divided2
    if options.type == 'DRW':
        #plot logA, logA = 0
        if options.d1 == 'loga2' and options.d2 == 'logl':
            #Convert to logA
            xs= (nu.log(2.)+xs+nu.log(1.-nu.exp(-1./nu.exp(ys))))/2.
        elif options.d2 == 'loga2' and options.d1 == 'logl':           
            #Convert to logA
            ys= (nu.log(2.)+ys+nu.log(1.-nu.exp(-1./nu.exp(xs))))/2.
        else:
            print "d1 and d2 have to be 0 or 1 (and not the same!) ..."
            print "Returning ..."
            return
    #stars
    starxs= nu.array([p[options.d1] for p in starparams.values()]).reshape(len(starparams))
    starys= nu.array([p[options.d2] for p in starparams.values()]).reshape(len(starparams))
    if options.expd1: starxs= nu.exp(starxs)
    elif not options.divided1 is None: starxs/= options.divided1
    if options.expd2: starys= nu.exp(starys)
    elif not options.divided2 is None: starys/= options.divided2
    if options.type == 'DRW':
        #plot logA, logA = 0
        if options.d1 == 'loga2' and options.d2 == 'logl':
            #Convert to logA
            starxs= (nu.log(2.)+starxs+nu.log(1.-nu.exp(-1./nu.exp(starys))))/2.
        elif options.d2 == 'loga2' and options.d1 == 'logl':           
            #Convert to logA
            starys= (nu.log(2.)+starys+nu.log(1.-nu.exp(-1./nu.exp(starxs))))/2.
        else:
            print "d1 and d2 have to be 0 or 1 (and not the same!) ..."
            print "Returning ..."
            return
    #RR Lyrae
    rrxs= nu.array([p[options.d1] for p in rrparams.values()]).reshape(len(rrparams))
    rrys= nu.array([p[options.d2] for p in rrparams.values()]).reshape(len(rrparams))
    if options.expd1: rrxs= nu.exp(rrxs)
    elif not options.divided1 is None: rrxs/= options.divided1
    if options.expd2: rrys= nu.exp(rrys)
    elif not options.divided2 is None: rrys/= options.divided2
    if options.type == 'DRW':
        #plot logA, logA = 0
        if options.d1 == 'loga2' and options.d2 == 'logl':
            #Convert to logA
            rrxs= (nu.log(2.)+rrxs+nu.log(1.-nu.exp(-1./nu.exp(rrys))))/2.
        elif options.d2 == 'loga2' and options.d1 == 'logl':           
            #Convert to logA
            rrys= (nu.log(2.)+rrys+nu.log(1.-nu.exp(-1./nu.exp(rrxs))))/2.
        else:
            print "d1 and d2 have to be 0 or 1 (and not the same!) ..."
            print "Returning ..."
            return
    #Plot
    xrange=[options.xmin,options.xmax]
    yrange=[options.ymin,options.ymax]
    bins= int(round(0.3*sc.sqrt(len(xs))))
    if options.type == 'powerlawSF':
        onedhistyweights=nu.ones(len(ys))/40.
    else:
        onedhistyweights=nu.ones(len(ys))/4000.
    bovy_plot.bovy_print()
    bovy_plot.bovy_plot(xs,ys,'b,',onedhists=True,
                        bins=bins,
                        onedhistynormed=False,
                        onedhistyweights=onedhistyweights,
                        xrange=xrange,
                        yrange=yrange,
                        onedhistec='b',
                        xlabel=options.xlabel,
                        ylabel=options.ylabel)
    #Stars
    bovy_plot.bovy_plot(starxs,starys,'k,',onedhists=True,
                        bins=bins,
                        xrange=xrange,
                        yrange=yrange,
                        overplot=True)
    #RR Lyrae
    bovy_plot.bovy_plot(rrxs,rrys,'r,',onedhists=True,
                        bins=bins,
                        onedhistec='r',
                        xrange=xrange,
                        yrange=yrange,
                        overplot=True)
    #Label
    if options.type == 'powerlawSF':
        bovy_plot.bovy_text(-4.4,1.15,r'$\mathrm{F/G\ stars}$',color='k')
        bovy_plot.bovy_text(-4.4,1.05,r'$\mathrm{QSOs}$',color='b')
        bovy_plot.bovy_text(-4.4,0.95,r'$\mathrm{RR\ Lyrae}$',color='r')
    elif options.type == 'DRW':
        bovy_plot.bovy_text(-4.4,2.88,r'$\mathrm{F/G\ stars}$',color='k')
        bovy_plot.bovy_text(-4.4,1.76,r'$\mathrm{QSOs}$',color='b')
        bovy_plot.bovy_text(-4.4,.64,r'$\mathrm{RR\ Lyrae}$',color='r')
    bovy_plot.bovy_end_print(options.plotfilename)

def get_options():
    usage = "usage: %prog [options] <savefilename>x3\n\nsavefilename= name of the file that holds the fits for 1) quasars, 2) stars, 3) RR Lyrae"
    parser = OptionParser(usage=usage)
    parser.add_option("-o",dest='plotfilename',default=None,
                      help="Name for the file that will hold the plot")
    parser.add_option("-t","--type",dest='type',default='powerlawSF',
                      help="Model that was fit ('powerlawSF' or 'DRW')")
    parser.add_option("--d1",dest='d1',default='logA',
                      help="index of first dimension")
    parser.add_option("--d2",dest='d2',default='gamma',
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
    parser.add_option("--xlabel",dest='xlabel',default=r'$\log A\ \mathrm{(amplitude\ at\ 1\ yr)}$',
                      help="xlabel")
    parser.add_option("--ylabel",dest='ylabel',default=r'$\gamma\ \mathrm{(power-law\ exponent)}$',
                      help="ylabel")
    return parser


if __name__ == '__main__':
    plotFitsall(get_options())
