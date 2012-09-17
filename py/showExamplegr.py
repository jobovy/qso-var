import sys
import numpy as nu
from numpy import array
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

def showExamplegr(filename='../data/SDSSJ203817.37+003029.8.fits',
                constraints=None,basefilename='SDSSJ203817.37+003029.8'):
    """
    NAME:
       showExamplegr
    PURPOSE:
       show an example of a power-law structure function GP covariance function
       fit to an SDSS quasar for g and r
    INPUT:
       filename - filename with the data
       constraints - if None, use all constraints, if [] use no constraints!
       basefilename - basefilename for plots
    OUTPUT:
       writes several plots
    HISTORY:
       2010-08-11 - Written - Bovy (NYU)
    """
    file= fu.table_fields(filename)
    mjd_g= nu.array(file.mjd_g)/365.25
    g= nu.array(file.g)
    err_g= nu.array(file.err_g)
    mjd_r= nu.array(file.mjd_r)/365.25
    r= nu.array(file.r)
    err_r= nu.array(file.err_r)

    mask= (mjd_g != 0)*(g < 20.6)*(g > 19.7)#Adjust for non-g
    g= g[mask]
    g-= nu.mean(g)
    err_g= err_g[mask]
    mjd_g= mjd_g[mask]
    mjd_g-= nu.amin(mjd_g)
    meanErr_g= nu.mean(err_g)

    r= r[mask]
    r-= nu.mean(r)
    err_r= err_r[mask]
    mjd_r= mjd_r[mask]
    mjd_r-= nu.amin(mjd_r)
    meanErr_r= nu.mean(err_r)
    
    meanErr_gr= nu.mean(nu.sqrt(err_r**2.+err_g**2.))

    nu.random.seed(4)
    nGP=5
    nx=201
    params_mean= ()
    from powerlawSFgr import covarFunc
    params= {'logGamma': array([-9.27821954]), 'logGammagr': array([-19.85172122]), 'gamma': array([ 0.45932629]), 'gammagr': array([ 0.29784021])}
    cf= covarFunc(**params)
    params_covar= (cf)
    ndata= len(g)
    if constraints is None:
        listx= [(t,'g') for t in mjd_g]
        listx.extend([(t,'r') for t in mjd_r])
        listy= [m for m in g]
        listy.extend([m for m in r])
        listy= nu.array(listy)
        noise= [m for m in err_g]
        noise.extend([m for m in err_r])
        noise= nu.array(noise)
        trainSet= trainingSet(listx=listx,listy=listy,noise=noise)
        constraints= trainSet
    else:
        constraints= nu.array([])

    useconstraints= constraints
    txs= nu.linspace(-0.1,6.5,nx)
    xs= [(txs[ii],'g') for ii in range(nx)]
    xs.extend([(txs[ii],'r') for ii in range(nx)])
    GPsamples= eval_gp(xs,mean,covar,(),params_covar,nGP=nGP,constraints=useconstraints,tiny_cholesky=.000001)
    thismean= calc_constrained_mean(xs,mean,params_mean,covar,params_covar,useconstraints)
    thiscovar= calc_constrained_covar(xs,covar,params_covar,useconstraints)
    plot.bovy_print()
    pyplot.plot(txs,GPsamples[0,:nx],'-',color='0.25')
    if isinstance(constraints,trainingSet):
        pyplot.plot(mjd_g,g,'k.',zorder=5,ms=10)
    plot.bovy_text(r'$\mathrm{'+basefilename+r'}$',title=True)
    #pyplot.fill_between(xs,thismean-sc.sqrt(sc.diagonal(thiscovar)),thismean+sc.sqrt(sc.diagonal(thiscovar)),color='.75')
    for ii in range(1,nGP):
        pyplot.plot(txs,GPsamples[ii,:nx],'-',color=str(0.25+ii*.5/(nGP-1)))
    pyplot.plot(txs,thismean[:nx],'k-',linewidth=2)
    if isinstance(constraints,trainingSet):
        pyplot.errorbar(6.15,-0.25,yerr=meanErr_g,color='k')
    pyplot.xlabel(r'$\mathrm{MJD-constant}\ [\mathrm{yr}]$')
    pyplot.ylabel(r'$g-\langle g\rangle\ [\mathrm{mag}]$')
    pyplot.xlim(-0.1,6.5)
    plot._add_ticks()
    plot.bovy_end_print(basefilename+'_full.ps')


    plot.bovy_print()
    pyplot.figure()
    pyplot.plot(txs,GPsamples[0,nx-1:-1],'-',color='0.25')
    if isinstance(constraints,trainingSet):
        pyplot.plot(mjd_r,r,'k.',zorder=5,ms=10)
    plot.bovy_text(r'$\mathrm{'+basefilename+r'}$',title=True)
    #pyplot.fill_between(xs,thismean-sc.sqrt(sc.diagonal(thiscovar)),thismean+sc.sqrt(sc.diagonal(thiscovar)),color='.75')
    for ii in range(1,nGP):
        pyplot.plot(txs,GPsamples[ii,nx-1:-1],'-',color=str(0.25+ii*.5/(nGP-1)))
    pyplot.plot(txs,thismean[nx-1:-1],'k-',linewidth=2)
    if isinstance(constraints,trainingSet):
        pyplot.errorbar(6.15,-0.25,yerr=meanErr_r,color='k')
    pyplot.xlabel(r'$\mathrm{MJD-constant}\ [\mathrm{yr}]$')
    pyplot.ylabel(r'$r-\langle r\rangle\ [\mathrm{mag}]$')
    pyplot.xlim(-0.1,6.5)
    plot._add_ticks()
    plot.bovy_end_print(basefilename+'_fullr.ps')


    plot.bovy_print()
    pyplot.figure()
    ii= 0
    colors= nu.array([GPsamples[ii,jj]-GPsamples[ii,jj+nx] for jj in range(nx)])
    colors= colors-nu.mean(colors)
    pyplot.plot(txs,colors,'-',color='0.25')
    if isinstance(constraints,trainingSet):
        plot.bovy_plot(mjd_g,g-r,'k.',zorder=5,ms=10,overplot=True)
    plot.bovy_text(r'$\mathrm{'+basefilename+r'}$',title=True)
    #pyplot.fill_between(xs,thismean-sc.sqrt(sc.diagonal(thiscovar)),thismean+sc.sqrt(sc.diagonal(thiscovar)),color='.75')
    for ii in range(1,nGP):
        colors= nu.array([GPsamples[ii,jj]-GPsamples[ii,jj+nx] for jj in range(nx)])
        colors= colors-nu.mean(colors)
        pyplot.plot(txs,colors,'-',color=str(0.25+ii*.5/(nGP-1)))
    plotthismean= nu.zeros(nx)
    for ii in range(nx):
        plotthismean[ii]= thismean[ii]-thismean[ii+nx]
    pyplot.plot(txs,plotthismean,'k-',linewidth=2)
    if isinstance(constraints,trainingSet):
        pyplot.errorbar(6.15,-0.18,yerr=meanErr_gr,color='k')
    pyplot.xlabel(r'$\mathrm{MJD-constant}\ [\mathrm{yr}]$')
    pyplot.ylabel(r'$g-r- \langle g - r \rangle\ [\mathrm{mag}]$')
    pyplot.xlim(-0.1,6.5)
    if isinstance(constraints,trainingSet):
        pyplot.ylim(-0.25,.25)
    else:
        pass #pyplot.ylim(-10.,10.)
    plot._add_ticks()
    plot.bovy_end_print(basefilename+'_color.ps')


    plot.bovy_print()
    pyplot.figure()
    ii= 2
    colors= nu.array([GPsamples[ii,jj]-GPsamples[ii,jj+nx] for jj in range(nx)])
    pyplot.plot(colors,GPsamples[ii,:nx],'-',color='0.25')
    if isinstance(constraints,trainingSet):
        pyplot.plot(g-r,g,'k.',zorder=5,ms=10)
    plot.bovy_text(r'$\mathrm{'+basefilename+r'}$',title=True)
    #pyplot.fill_between(xs,thismean-sc.sqrt(sc.diagonal(thiscovar)),thismean+sc.sqrt(sc.diagonal(thiscovar)),color='.75')
    for ii in range(1,nGP):
        colors= nu.array([GPsamples[ii,jj]-GPsamples[ii,jj+nx] for jj in range(nx)])
        pyplot.plot(colors,GPsamples[ii,0:nx],'-',color=str(0.25+ii*.5/(nGP-1)))
    colors= nu.array([thismean[jj]-thismean[jj+nx] for jj in range(nx)])
    pyplot.plot(colors,thismean[:nx],'k-',linewidth=2)
    if isinstance(constraints,trainingSet):
        pyplot.errorbar(.13,-0.3,yerr=meanErr_g,xerr=meanErr_gr,color='k')
    pyplot.xlabel(r'$g-r-\langle g-r\rangle\ [\mathrm{mag}]$')
    pyplot.ylabel(r'$g-\langle g\rangle\ [\mathrm{mag}]$')
    pyplot.xlim(-.2,.2)
    plot._add_ticks()
    plot.bovy_end_print(basefilename+'_ggr.ps')

    return

    plot.bovy_print()
    pyplot.plot(xs,GPsamples[0,:],'-',color='0.25')
    if not constraints == []:
        plot.bovy_plot(mjd_g,g,'k.',zorder=5,ms=10,overplot=True)
    #plot.bovy_text(r'$\mathrm{SDSSJ203817.37+003029.8}$',title=True)
    #pyplot.fill_between(xs,thismean-sc.sqrt(sc.diagonal(thiscovar)),thismean+sc.sqrt(sc.diagonal(thiscovar)),color='.75')
    for ii in range(1,nGP):
        pyplot.plot(xs,GPsamples[ii,:],'-',color=str(0.25+ii*.5/(nGP-1)))
    #pyplot.plot(xs,thismean,'k-',linewidth=2)
    if not constraints == []:
        pyplot.errorbar(6.15,-0.1,yerr=meanErr_g,color='k')
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
        showExamplegr(band=sys.argv[2],basefilename=sys.argv[1],constraints=constraints)
    elif len(sys.argv) > 1:
        if 'like' in sys.argv[1]:
            constraints=[]
        else:
            constraints= None
        showExamplegr(basefilename=sys.argv[1],constraints=constraints)
    else:
        showExamplegr()
