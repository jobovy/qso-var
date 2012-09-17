import os, os.path
import cPickle as pickle
from optparse import OptionParser
import numpy as nu
import scipy as sc
import matplotlib
from galpy.util import bovy_plot
import xdtarget
def plotXDall(parser):
    nu.random.seed(1)
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
    if os.path.exists(args[1]):
        savefile= open(args[1],'rb')  
        starxamp= pickle.load(savefile)
        starxmean= pickle.load(savefile)
        starxcovar= pickle.load(savefile)
        savefile.close()
    else:
        print args[1]+" does not exist ..."
        print "Returning ..."
        return
    if os.path.exists(args[2]):
        savefile= open(args[2],'rb')  
        rrxamp= pickle.load(savefile)
        rrxmean= pickle.load(savefile)
        rrxcovar= pickle.load(savefile)
        savefile.close()
    else:
        print args[2]+" does not exist ..."
        print "Returning ..."
        return
    if options.nsamplesstar is None: options.nsamplesstar= options.nsamples
    if options.nsamplesrrlyrae is None: options.nsamplesrrlyrae= options.nsamples
    #Load XD object in xdtarget
    xdt= xdtarget.xdtarget(amp=xamp,mean=xmean,covar=xcovar)
    out= xdt.sample(nsample=options.nsamples)
    xdt= xdtarget.xdtarget(amp=starxamp,mean=starxmean,covar=starxcovar)
    starout= xdt.sample(nsample=options.nsamplesstar)
    xdt= xdtarget.xdtarget(amp=rrxamp,mean=rrxmean,covar=rrxcovar)
    rrout= xdt.sample(nsample=options.nsamplesrrlyrae)
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
    #stars
    if options.expd1: starxs= nu.exp(starout[:,options.d1])
    elif not options.divided1 is None: starxs= starout[:,options.d1]/options.divided1
    else: starxs= starout[:,options.d1]
    if options.expd2: starys= nu.exp(starout[:,options.d2])
    elif not options.divided2 is None: starys= starout[:,options.d2]/options.divided2
    else: starys= starout[:,options.d2]
    if options.type == 'DRW':
        #plot logA, logA = 0
        if options.d1 == 0 and options.d2 == 1: 
            #Convert to logA
            starxs= (nu.log(2.)+starxs+nu.log(1.-nu.exp(-1./nu.exp(starys))))/2.
        elif options.d1 ==1 and options.d2 == 0:
            #Convert to logA
            starys= (nu.log(2.)+starys+nu.log(1.-nu.exp(-1./nu.exp(starxs))))/2.
        else:
            print "d1 and d2 have to be 0 or 1 (and not the same!) ..."
            print "Returning ..."
            return
    #RR Lyrae
    if options.expd1: rrxs= nu.exp(rrout[:,options.d1])
    elif not options.divided1 is None: rrxs= rrout[:,options.d1]/options.divided1
    else: rrxs= rrout[:,options.d1]
    if options.expd2: rrys= nu.exp(rrout[:,options.d2])
    elif not options.divided2 is None: rrys= rrout[:,options.d2]/options.divided2
    else: rrys= rrout[:,options.d2]
    if options.type == 'DRW':
        #plot logA, logA = 0
        if options.d1 == 0 and options.d2 == 1: 
            #Convert to logA
            rrxs= (nu.log(2.)+rrxs+nu.log(1.-nu.exp(-1./nu.exp(rrys))))/2.
        elif options.d1 ==1 and options.d2 == 0:
            #Convert to logA
            rrys= (nu.log(2.)+rrys+nu.log(1.-nu.exp(-1./nu.exp(rrxs))))/2.
        else:
            print "d1 and d2 have to be 0 or 1 (and not the same!) ..."
            print "Returning ..."
            return
    #Plot
    xrange=[options.xmin,options.xmax]
    yrange=[options.ymin,options.ymax]
    data= sc.array([xs,ys]).T
    bins= int(round(0.3*sc.sqrt(options.nsamples)))
    hist, edges= sc.histogramdd(data,bins=bins,range=[xrange,yrange])
    #Censor hist ASSUMES gamma=[0.,1.2], logA=[-9.21/2.,0.] for powerlawSF
    x= nu.zeros((bins,bins))
    y= nu.zeros((bins,bins))
    for bb in range(bins):
        x[:,bb]= nu.linspace(options.xmin,options.xmax,bins)
        y[bb,:]= nu.linspace(options.ymin,options.ymax,bins)  
    #mask
    if options.type =='powerlawSF':
        hist[(y < 0.1) * (x > -3.) * ( x < -1.5)]= nu.nan
        hist[(x < -3.)]= nu.nan
        hist[(x > -2.) * (y < (0.25*x+0.6))]= nu.nan
        onedhistyweights=nu.ones(len(ys))/100.
    elif options.type == 'DRW':
        hist[(y < (-4.223*(x+2)-10))]= nu.nan
        hist[(y < -4.153)*(y < (58.47*(x+2.1)-10.))]= nu.nan
        hist[(y > -4.153)*(y < (2.93*x+2.)-1.)]= nu.nan
        onedhistyweights=nu.ones(len(ys))/2500.
    bovy_plot.bovy_print()
    #First just plot contours
    cdict = {'red': ((.0, 1.0, 1.0),
                     (1.0, 1.0, 1.0)),
             'green': ((.0, 1.0, 1.0),
                       (1.0, 1.0, 1.0)),
             'blue': ((.0, 1.0, 1.0),
                      (1.0, 1.0, 1.0))}
    allwhite = matplotlib.colors.LinearSegmentedColormap('allwhite',cdict,256)
    bovy_plot.scatterplot(xs,ys,'b,',onedhists=True,
                          bins=bins,
                          cmap=allwhite,
                          onedhistynormed=False,
                          onedhistyweights=onedhistyweights,
                          xrange=xrange,
                          yrange=yrange,
                          onedhistec='b',
                          xlabel=options.xlabel,
                          ylabel=options.ylabel)
    bovy_plot.scatterplot(starxs,starys,'k,',onedhists=True,
                          bins=bins,
                          cmap=allwhite,
                          xrange=xrange,
                          yrange=yrange,
                          overplot=True)
    bovy_plot.scatterplot(rrxs,rrys,'r,',onedhists=True,
                          bins=bins,
                          cmap=allwhite,
                          onedhistec='r',
                          xrange=xrange,
                          yrange=yrange,
                          overplot=True)
    hist/= nu.nansum(hist)
    #Custom colormap
    cdict = {'red': ((.0, 1.0, 1.0),
                     (1.0, .0, .0)),
             'green': ((0.0, 1.0, 1.0),
                       (1.0, .0, .0)),
             'blue': ((0.0, 1.0, 1.0),
                      (1.0, 1.0, 1.0))}
    my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    bovy_plot.scatterplot(xs,ys,'b,',onedhists=False,contours=False,
                          levels=[1.01],
                          bins=bins,
                          cmap=my_cmap,
                          hist=hist,edges=edges,
                          onedhistynormed=False,
                          onedhistyweights=onedhistyweights,
                          xrange=xrange,
                          yrange=yrange,
                          overplot=True)
    #Stars
    data= sc.array([starxs,starys]).T
    hist, edges= sc.histogramdd(data,bins=bins,range=[xrange,yrange])
    if options.type == 'powerlawSF':
        hist[(x > -2.5)]= nu.nan
        hist[(x < -2.5)*(y > (-.19*(x+2.5)))]= nu.nan
    elif options.type == 'DRW':
        hist[(y >= (-4.223*(x+2)-10))]= nu.nan
    hist/= nu.nansum(hist)
    bovy_plot.scatterplot(starxs,starys,'k,',onedhists=True,contours=False,
                          levels=[1.01],#HACK such that outliers aren't plotted
                          bins=bins,
                          hist=hist,edges=edges,
                          xrange=xrange,
                          yrange=yrange,
                          overplot=True)
    #RR Lyrae
    data= sc.array([rrxs,rrys]).T
    hist, edges= sc.histogramdd(data,bins=bins,range=[xrange,yrange])
    if options.type == 'powerlawSF':
        hist[(x < -2.5)]= nu.nan
        hist[(x > -2.5)*(y > ((x+2.5)/1.5)**9.*1.15+0.1)]= nu.nan
        #hist[(x > -2.5)*(y > (.25*x+.6))]= nu.nan
    elif options.type == 'DRW':
        hist[(y < -4.153)*(y >= (58.47*(x+2.1)-10.))*(y > -7.5)]= nu.nan
        hist[(y < -7.5)*(x < -2.1)]= nu.nan
        hist[(y > -4.153)*(y >= (2.93*x+2.-1.5))]= nu.nan
    #Custom colormap
    cdict = {'red': ((.0, 1.0, 1.0),
                     (1.0, 1.0, 1.0)),
             'green': ((0.0, 1.0, 1.0),
                       (1.0, .0, .0)),
             'blue': ((0.0, 1.0, 1.0),
                      (1.0, .0, .0))}
    my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    hist/= nu.nansum(hist)
    bovy_plot.scatterplot(rrxs,rrys,'r,',onedhists=False,contours=False,
                          levels=[1.01],#HACK such that outliers aren't plotted
                          bins=bins,
                          cmap=my_cmap,
                          hist=hist,edges=edges,
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
    usage = "usage: %prog [options] <savefilename>x3\n\nsavefilename= name of the file that holds the XD fit for 1) quasars, 2) stars, 3) RR Lyrae"
    parser = OptionParser(usage=usage)
    parser.add_option("-o",dest='plotfilename',default=None,
                      help="Name for the file that will hold the plot")
    parser.add_option("-t","--type",dest='type',default='powerlawSF',
                      help="Type of model that was fit (powerlawSF or DRW)")
    parser.add_option("-n","--nsamples",dest='nsamples',default=10000,
                      type='int',
                      help="Number of samples to plot")
    parser.add_option("--nstar",dest='nsamplesstar',default=None,
                      type='int',
                      help="Number of samples to plot for stars (default=same as quasars)")
    parser.add_option("--nrrlyrae",dest='nsamplesrrlyrae',default=None,
                      type='int',
                      help="Number of samples to plot for RR Lyrae (default: same as quasars)")
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
    parser.add_option("--xlabel",dest='xlabel',default=r'$\log A\ \mathrm{(amplitude\ at\ 1\ yr)}$',
                      help="xlabel")
    parser.add_option("--ylabel",dest='ylabel',default=r'$\gamma\ \mathrm{(power-law\ exponent)}$',
                      help="ylabel")
    return parser


if __name__ == '__main__':
    plotXDall(get_options())
