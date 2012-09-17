import os, os.path
import cPickle as pickle
from optparse import OptionParser
import numpy as nu
import astrometry.util.pyfits_utils as pyfits_utils
import galpy.util.bovy_plot as bovy_plot
from varqso import VarQso
def compareMagFluxFits(parser):
    (options,args)= parser.parse_args()
    if len(args) == 0:
        parser.print_help()
        return
    params= []
    for filename in args:
        if os.path.exists(filename):
            savefile= open(filename,'rb')  
            params.append(pickle.load(savefile))
            savefile.close()
        else:
            print filename+" does not exist ..."
            print "Returning ..."
            return
    if options.plottype == 'AA':
        ys= []
        xs= []
        for key in params[1].keys():
            try:
                ys.append(params[0][key]['logA']/2.)
                xs.append(params[1][key]['logA']/2.)
            except KeyError:
                continue
        xs= nu.array(xs).reshape(len(xs))
        ys= nu.array(ys).reshape(len(xs))
        ys= xs-ys-nu.log(nu.log(10.)/2.5)
        xrange=[-9.21/2.,0.]
        yrange=[-.25,.25]
        xlabel= r'$\log A^{\mathrm{flux}}_'+options.band+r'\ \mathrm{(amplitude\ at\ 1\ yr)}$'
        ylabel= r'$\log A^{\mathrm{flux}}_'+options.band+r'-\log A^{\mathrm{mag}}_'+options.band+r'- \log\left(\frac{\log 10}{2.5}\right)$'
    elif options.plottype == 'gg':
        ys= []
        xs= []
        for key in params[1].keys():
            try:
                ys.append(params[0][key]['gamma'])
                xs.append(params[1][key]['gamma'])
            except KeyError:
                continue
        xs= nu.array(xs).reshape(len(xs))
        ys= nu.array(ys).reshape(len(xs))
        print len(xs)
        ys= xs-ys
        xrange=[0.,1.2]
        yrange=[-.25,.25]
        xlabel= r'$\gamma^{\mathrm{flux}}_'+options.band+r'\ \mathrm{(power-law\ exponent)}$'
        ylabel= r'$\gamma^{\mathrm{flux}}_'+options.band+r'- \gamma^{\mathrm{mag}}_'+options.band+'$'
    elif options.plottype == 'loglike2':
        ys= []
        xs= []
        nepochs= []
        cnt= 0
        for key in params[1].keys():
            #cnt+= 1
            #print cnt
            #if cnt > 10: break
            #Get the number of epochs
            v= VarQso(os.path.join('../data/s82qsos/',key))
            try:
                ys.append(params[0][key]['loglike']
                          -v.nepochs(options.band)*nu.log(nu.log(10.)/2.5))
                nepochs.append(v.nepochs(options.band))
                xs.append(params[1][key]['loglike'])
            except KeyError:
                continue
        xs= -nu.array(xs).reshape(len(xs))
        ys= -nu.array(ys).reshape(len(xs))
        nepochs= nu.array(nepochs).reshape(len(nepochs))
        ys/= nepochs
        xs/= nepochs
        ys= xs-ys
        xrange=[-3.,0.]
        yrange=[-.1,.1]
        xlabel= r'$\log \mathcal{L}^{\mathrm{flux}}_{'+options.band+r',\mathrm{red}}\ \mathrm{(amplitude\ at\ 1\ yr)}$'
        ylabel= r'$\log \mathcal{L}^{\mathrm{flux}}_{'+options.band+r',\mathrm{red}}- \log \mathcal{L}^{\mathrm{mag}}_{'+options.band+',\mathrm{red}}'+r'- \log\left(\frac{\log 10}{2.5}\right)$'
    bovy_plot.bovy_print()
    bovy_plot.scatterplot(xs,ys,'k,',onedhists=True,
                          yrange=yrange,
                          xrange=xrange,bins=31,
                          xlabel=xlabel,
                          ylabel=ylabel)
    bovy_plot.bovy_plot(nu.array(xrange),[0.,0.],'0.5',
                        overplot=True)
    bovy_plot.bovy_end_print(options.plotfilename)
    return None

def get_options():
    usage = "usage: %prog [options] <savefilename_mags> <savefilename_flux>\n\nsavefilename_mags= name of the file that holds the fits for magnitudes\nsavefilename_flux= name of the file that holds the fits for relative fluxes"
    parser = OptionParser(usage=usage)
    parser.add_option("-b","--band",dest='band',default='r',
                      help="band(s) to fit")
    parser.add_option("-t","--type",dest='type',default='powerlawSF',
                      help="Type of model (powerlawSF or DRW, DRW not implemented yet)")
    parser.add_option("-o",dest='plotfilename',default=None,
                      help="Name for the file that will hold the plot")
    parser.add_option("--plottype",dest='plottype',default='AA',
                      help="Type of plot (AA, gg, loglike2")
    return parser

if __name__ == '__main__':
    compareMagFluxFits(get_options())
