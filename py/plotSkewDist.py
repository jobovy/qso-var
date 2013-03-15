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
        allskews[ii,:]= -skews[key] #go to regular definition
        allgaussskews[ii,:]= -gaussskews[key]
    #Statistic
    q= 0.99
    statistic= numpy.median
    if not options.indx is None:
        print "indx option not allowed"
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
    tauindx= 5
    print "Showing tau %f" % (taus[tauindx]*365.25)
    bovy_plot.bovy_print(fig_width=7.)
    bovy_plot.bovy_hist(allskews[:,tauindx],bins=31,
                        color='k',normed=True,histtype='step',
                        xlabel=r'$\mathrm{skew}(\tau = %i\ \mathrm{days})$' % (int(taus[tauindx]*365.25)))
    bovy_plot.bovy_plot([numpy.median(allskews[:,tauindx]),numpy.median(allskews[:,tauindx])],
                        [0.,10.],'k-',overplot=True)
    bovy_plot.bovy_plot([numpy.mean(allskews[:,tauindx]),numpy.mean(allskews[:,tauindx])],
                        [0.,10.],'k--',overplot=True)
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
