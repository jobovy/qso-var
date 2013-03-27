import os, os.path
import cPickle as pickle
from optparse import OptionParser
import numpy
from galpy.util import bovy_plot
from matplotlib import pyplot
def plotSkew(parser):
    (options,args)= parser.parse_args()
    if len(args) == 0:
        parser.print_help()
        return
    savefilename= args[0]
    if os.path.exists(savefilename):
        savefile= open(savefilename,'rb')
        skews= pickle.load(savefile)
        gaussskews= pickle.load(savefile)
        type= pickle.load(savefile)
        band= pickle.load(savefile)
        mean= pickle.load(savefile)
        taus= pickle.load(savefile)      
        savefile.close()
    else:
        parser.print_help()
        return
    #Accumulate
    keys= skews.keys()
    allskews= numpy.zeros((len(keys),len(taus)))
    allgaussskews= numpy.zeros((len(keys),gaussskews[keys[0]].shape[0],
                                len(taus)))
    for ii, key in enumerate(keys):
        allskews[ii,:]= skews[key]
        allgaussskews[ii,:]= gaussskews[key] #Don't flip to get correct wedges
    #Statistic
    q= 0.99
    statistic= numpy.median
    if not options.indx is None:
        sigma= quantile(allgaussskews[options.indx,:,:],q=q)
        bovy_plot.bovy_print(fig_width=7.)
        bovy_plot.bovy_plot(taus*365.25,allskews[options.indx,:],'k-',
                            xlabel=r'$\mathrm{lag}\ \tau\ [\mathrm{days}]$',
                            ylabel=r'$\mathrm{skew}(\tau)$',
                            zorder=5,
                            yrange=[-1.,1.])
        pyplot.fill_between(taus*365.25,sigma[0,:],sigma[1,:],color='0.75',zorder=0)
        bovy_plot.bovy_plot(taus*365.25,statistic(allgaussskews[options.indx,:,:],
                                                   axis=0),'-',
                            overplot=True,color='0.5')
        bovy_plot.bovy_end_print(options.plotfilename)
        return None
    indx= numpy.all(numpy.isnan(allskews),axis=1)
    allskews= allskews[True-indx,:]
    allgaussskews= allgaussskews[True-indx,:,:]
    #Median
    medianskew= statistic(allskews,axis=0)
    mediangaussskew= statistic(allgaussskews,axis=0)
    #Determine 1-sigma
    sigma= quantile(mediangaussskew,q=q)
    #Plot
    bovy_plot.bovy_print(fig_width=7.)
    bovy_plot.bovy_plot(taus*365.25,medianskew,'k-',
                        xlabel=r'$\mathrm{lag}\ \tau\ [\mathrm{days}]$',
                        ylabel=r'$\mathrm{skew}(\tau)$',
                        zorder=5,
                        yrange=[-1.,1.])
    pyplot.fill_between(taus*365.25,sigma[0,:],sigma[1,:],color='0.75',zorder=0)
    bovy_plot.bovy_plot(taus*365.25,statistic(mediangaussskew,axis=0),'-',
                        overplot=True,color='0.5')
    bovy_plot.bovy_end_print(options.plotfilename)
    return None

def quantile(t,q=0.68):
    """calculate the quantile of a distribution for this problem,
    t= ngauss,nt, quantile compute for each t"""
    out= numpy.zeros((2,t.shape[1]))
    for ii in range(out.shape[1]):
        sortedt= sorted(t[:,ii])
        low= sortedt[int(numpy.floor((0.5-q/2.)*len(t)))]
        high= sortedt[int(numpy.floor((0.5+q/2.)*len(t)))]
        out[0,ii]= low
        out[1,ii]= high
    return out

def get_options():
    usage = "usage: %prog [options] <savefilename>\n\nsavefilename= name of the file that holds the skews"
    parser = OptionParser(usage=usage)
    parser.add_option("-o",dest='plotfilename',default=None,
                      help="Name for plotfile")
    parser.add_option("-i","--indx",dest='indx',default=None,type='float',
                      help="Just plot this index")
    parser.add_option("-b","--band",dest='band',default='r',
                      help="band(s) to sample")
    parser.add_option("-t","--type",dest='type',default='powerlawSF',
                      help="Type of model to sample (powerlawSF, powerlawSFratios, or DRW)")
    parser.add_option("--mean",dest='mean',default='zero',
                      help="Type of mean to sample (zero, const)")
    parser.add_option("-f","--fitsfile",dest='fitsfile',
                      default=None,
                      help="File that holds the best-fits")
    return parser
    
if __name__ == '__main__':
    plotSkew(get_options())
