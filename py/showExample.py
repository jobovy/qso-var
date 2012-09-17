import sys
import numpy as nu
from numpy import *
from astrometry.util import pyfits_utils as fu
import pyfits
import galpy.util.bovy_plot as plot
from matplotlib import pyplot
from gp.eval_gp import *

def mean(x,params):
    return 0

def covar(x,y,params):
    """covar: covariance function"""
    return params.evaluate(x,y)

def covar(x,y,params):
    """covar: covariance function"""
    #sig, A, gamma= .1,.24,0.38
    sig, A, gamma= nu.exp(-2.1762761),nu.exp(-3.02045715),0.5868293
    return sig*(1.-A/sig/2.*sc.fabs(x-y)**gamma)

def showExample(filename='../data/SDSSJ203817.37+003029.8.fits',band='g',
                constraints=None,basefilename='SDSSJ203817.37+003029.8'):
    """
    NAME:
       showExample
    PURPOSE:
       show an example of a power-law structure function GP covariance function
       fit to an SDSS quasar
    INPUT:
       filename - filename with the data
       band - band to consider
       constraints - if None, use all constraints, if [] use no constraints!
       basefilename - basefilename for plots
    OUTPUT:
       writes several plots
    HISTORY:
       2010-06-20 - Written - Bovy (NYU)
    """
    file= fu.table_fields(filename)
    if band == 'u':
        mjd_g= nu.array(file.mjd_u)/365.25
        g= nu.array(file.u)
        err_g= nu.array(file.err_u)
    elif band == 'g':
        mjd_g= nu.array(file.mjd_g)/365.25
        g= nu.array(file.g)
        err_g= nu.array(file.err_g)
    elif band == 'r':
        mjd_g= nu.array(file.mjd_r)/365.25
        g= nu.array(file.r)
        err_g= nu.array(file.err_r)
    elif band == 'i':
        mjd_g= nu.array(file.mjd_i)/365.25
        g= nu.array(file.i)
        err_g= nu.array(file.err_i)
    elif band == 'z':
        mjd_g= nu.array(file.mjd_z)/365.25
        g= nu.array(file.z)
        err_g= nu.array(file.err_z)

    mask= (mjd_g != 0)*(g < 20.6)*(g > 19.7)#Adjust for non-g
    g= g[mask]
    g-= nu.mean(g)
    err_g= err_g[mask]
    mjd_g= mjd_g[mask]
    mjd_g-= nu.amin(mjd_g)
    meanErr_g= nu.mean(err_g)


    nu.random.seed(4)
    nGP=5
    nx=201
    params_mean= ()
    from powerlawSFgr import covarFunc
    params= {'logA': array([ 1.63714183]), 'gamma': array([ 0.60216581])}
    params= {'gamma': array([ 0.49500723]), 'logA': array([-3.36044037])}
    cf= covarFunc(**params)
    params_covar= (cf)
    ndata= len(g)
    if constraints is None:
        listx= mjd_g
        listy= g
        noise= err_g
        trainSet= trainingSet(listx=listx,listy=listy,noise=noise)
        constraints= trainSet
    else:
        constraints= None

    useconstraints= constraints
    xs= nu.linspace(-0.1,6.5,nx)
    GPsamples= eval_gp(xs,mean,covar,(),params_covar,nGP=nGP,constraints=useconstraints)
    thismean= calc_constrained_mean(xs,mean,params_mean,covar,params_covar,useconstraints)
    thiscovar= calc_constrained_covar(xs,covar,params_covar,useconstraints)
    plot.bovy_print()
    pyplot.plot(xs,GPsamples[0,:],'-',color='0.25')
    if not constraints is None:
        plot.bovy_plot(mjd_g,g,'k.',zorder=5,ms=10,overplot=True)
    plot.bovy_text(r'$\mathrm{'+basefilename+r'}$',title=True)
    #pyplot.fill_between(xs,thismean-sc.sqrt(sc.diagonal(thiscovar)),thismean+sc.sqrt(sc.diagonal(thiscovar)),color='.75')
    for ii in range(nGP):
        pyplot.plot(xs,GPsamples[ii,:],'-',color=str(0.25+ii*.5/(nGP-1)))
    #pyplot.plot(xs,thismean,'k-',linewidth=2)
    if not constraints == []:
        pyplot.errorbar(6.15,-0.25,yerr=meanErr_g,color='k')
    pyplot.xlabel(r'$\mathrm{MJD-constant}\ [\mathrm{yr}]$')
    pyplot.ylabel(r'$m-\langle m\rangle\ [\mathrm{mag}]$')
    pyplot.xlim(-0.1,6.5)
    pyplot.ylim(-1.,1.)
    plot._add_ticks()
    plot.bovy_end_print(basefilename+'_full.ps')

    plot.bovy_print()
    pyplot.plot(xs,GPsamples[0,:],'-',color='0.25')
    if not constraints is None:
        plot.bovy_plot(mjd_g,g,'k.',zorder=5,ms=10,overplot=True)
    #plot.bovy_text(r'$\mathrm{SDSSJ203817.37+003029.8}$',title=True)
    #pyplot.fill_between(xs,thismean-sc.sqrt(sc.diagonal(thiscovar)),thismean+sc.sqrt(sc.diagonal(thiscovar)),color='.75')
    for ii in range(1,nGP):
        pyplot.plot(xs,GPsamples[ii,:],'-',color=str(0.25+ii*.5/(nGP-1)))
    #pyplot.plot(xs,thismean,'k-',linewidth=2)
    if not constraints == []:
        pyplot.errorbar(6.15,-0.25,yerr=meanErr_g,color='k')
    pyplot.xlabel(r'$\mathrm{MJD-constant}\ [\mathrm{yr}]$')
    pyplot.ylabel(r'$m-\langle m\rangle\ [\mathrm{mag}]$')
    pyplot.xlim(3,6.5)
    plot._add_ticks()
    plot.bovy_end_print(basefilename+'_zoom.ps')


    plot.bovy_print()
    pyplot.figure()
    for ii in range(nGP):
        pyplot.loglog(sc.arange(1.,len(GPsamples[ii,:])/2)*(xs[1]-xs[0]),
                      2.*sc.var(GPsamples[ii,:])-2.*sc.correlate(GPsamples[ii,:]-sc.mean(GPsamples[ii,:]),GPsamples[ii,:]-sc.mean(GPsamples[ii,:]),"same")[1:len(GPsamples[ii,:])/2][::-1]/len(GPsamples[ii,:]),
                      color=str(0.25+ii*.5/(nGP-1)))
    xline= [(xs[1]-xs[0]),xs[len(xs)/2]]
    pyplot.loglog(xline,nu.exp(-3.02045715)*nu.array(xline)**(0.5868293),'k--')
    pyplot.xlabel(r'$\Delta t\ [\mathrm{yr}]$')
    pyplot.ylabel(r'$\mathrm{structure\ function}$')
    #plot.bovy_text(r'$\mathrm{SDSSJ203817.37+003029.8}$',title=True)
    plot.bovy_end_print(basefilename+'_structfunc.ps')


if __name__ == '__main__':
    if len(sys.argv) > 2:
        if 'like' in sys.argv[1]:
            constraints=[]
        else:
            constraints= None
        showExample(band=sys.argv[2],basefilename=sys.argv[1],constraints=constraints)
    elif len(sys.argv) > 1:
        if 'like' in sys.argv[1]:
            constraints=[]
        else:
            constraints= None
        showExample(basefilename=sys.argv[1],constraints=constraints)
    else:
        showExample()
