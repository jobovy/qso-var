import sys
import numpy as nu
import scipy as sc
import inspect
from scipy import signal, optimize
from matplotlib import pyplot
#from astrometry.util import pyfits_utils as fu
import pyfits
import galpy.util.bovy_plot as bovy_plot
try:
    from flexgp.trainingSet import trainingSet
    from flexgp.trainGP import trainGP, marginalLikelihood, pack_params
    from flexgp.eval_gp import *
    from flexgp.sampleGP import sampleGP
except ImportError:
    from gp.trainingSet import trainingSet
    from gp.trainGP import trainGP, marginalLikelihood, pack_params
    from gp.eval_gp import *
    from gp.sampleGP import sampleGP
_BANDDICT= {'u':0,'g':1,'r':2,'i':3,'z':4}
_ERASESTR= "                                                                                "
def mean_func(x,params):
    return params.evaluate(x)
def covar_func(x,y,params):
    """covar: covariance function"""
    return params.evaluate(x,y)
class VarQso():
    """Class that represents a quasar light-curve"""
    def __init__(self,*args,**kwargs):
        """
        NAME:
           __init__
        PURPOSE:
           initialize
        INPUT:
           filename - name of the file that holds the data for this object
           band - bands to load
           flux= if True, load fluxes (default=load luptitudes; only relevant 
                 for file-based initialization)
        OUTPUT:
        HISTORY:
           2010-12-21 - Written - Bovy (NYU)
           2011-07-01 - Added flux= keyword - Bovy (NYU)
        """
        if kwargs.has_key('band'):
            band= kwargs['band']
        else:
            band= 'ugriz'
        if kwargs.has_key('medianize'):
            medianize= kwargs['medianize']
        else:
            medianize= True
        if kwargs.has_key('flux') and kwargs['flux']:
            self.flux= True
        else:
            self.flux= False
        self.mjd= {}
        self.m= {}
        self.err_m= {}
        self.meanErr= {}
        if len(args) == 1: #File-based initialization
            file= _load_fits(args[0])
            #file= fu.table_fields(args[0])
            for b in band:
                self._init_oneband(file,b,medianize,self.flux)
        elif len(args) == 3: #(mjd,band),m,e_m
            for b in band:
                self._init_arrays(args,b,medianize,self.flux)
        try:
            minMJD= nu.array([nu.amin(self.mjd[b]) for b in band])
        except ValueError:
            minMJD=nu.array([0.])
        minMJD= nu.amin(minMJD)
        for b in band:
            self.mjd[b]= self.mjd[b]-minMJD

    def _init_arrays(self,args,band,medianize,flux):
        if ((not isinstance(args[0][0],(int,float)) and 
             len(args[0][0]) == 1) or \
                isinstance(args[0][0],(int,float))):
            self.mjd[band]= nu.array([x for x in args[0]])
            self.m[band]= nu.array([x for x in args[1]])
            self.err_m[band]= nu.array([x for x in args[2]])
        else:
            self.mjd[band]= nu.array([x[0] for x in args[0] if x[1] == band])
            self.m[band]= nu.array([args[1][ii] for ii in range(len(args[0]))\
                                        if args[0][ii][1] == band])
            self.err_m[band]= nu.array([args[2][ii] for ii in range(len(args[0]))\
                                            if args[0][ii][1] == band])
        if medianize:
            self.medianize(band)
        if flux:
            from sdsspy.util import lups2nmgy
            nmgy, ivar= lups2nmgy(self.m[band],err=self.err_m[band],
                                  band=_BANDDICT[band])
            self.m[band]= nmgy
            self.err_m[band]= 1./nu.sqrt(ivar)
            self.err_m[band]/= nu.mean(self.m[band]) #BOVY: not quite correct
            self.m[band]/= nu.mean(self.m[band])
        self.meanErr[band]= nu.sqrt(nu.mean(self.err_m[band]**2.))

    def _init_oneband(self,file,band,medianize,flux):
        if band == 'u':
            self.mjd[band]= nu.array(file.mjd_u)
            self.m[band]= nu.array(file.u)
            self.err_m[band]= nu.array(file.err_u)
        if band == 'g':
            self.mjd[band]= nu.array(file.mjd_g)
            self.m[band]= nu.array(file.g)
            self.err_m[band]= nu.array(file.err_g)
        if band == 'r':
            self.mjd[band]= nu.array(file.mjd_r)
            self.m[band]= nu.array(file.r)
            self.err_m[band]= nu.array(file.err_r) 
        if band == 'i':
            self.mjd[band]= nu.array(file.mjd_i)
            self.m[band]= nu.array(file.i)
            self.err_m[band]= nu.array(file.err_i)
        if band == 'z':
            self.mjd[band]= nu.array(file.mjd_z)
            self.m[band]= nu.array(file.z)
            self.err_m[band]= nu.array(file.err_z)
        mask= (self.mjd[band] != 0) #Gets rid of co-add runs
        self.mjd[band]= self.mjd[band][mask]
        self.m[band]= self.m[band][mask]
        self.err_m[band]= self.err_m[band][mask]
        self.mjd[band]= self.mjd[band]/365.25
        if medianize:
            self.medianize(band)
        if flux:
            from sdsspy.util import lups2nmgy
            nmgy, ivar= lups2nmgy(self.m[band],err=self.err_m[band],
                                  band=_BANDDICT[band])
            self.m[band]= nmgy
            self.err_m[band]= 1./nu.sqrt(ivar)
            self.err_m[band]/= nu.mean(self.m[band]) #BOVY: not quite correct
            self.m[band]/= nu.mean(self.m[band])
        self.meanErr[band]= nu.sqrt(nu.mean(self.err_m[band]**2.))

    def medianize(self,band):
        """
        NAME:
           medianize
        PURPOSE:
           medianize the light-curve to get rid of outliers
        INPUT:
           band - band to medianize
        OUTPUT:
           (none)
        HISTORY:
           2010-12-21 - 
        """
        med= signal.medfilt(self.m[band],kernel_size=7)
        mask= (nu.fabs(med-self.m[band]) <= 0.25)
        self.mjd[band]= self.mjd[band][mask]
        self.m[band]= self.m[band][mask]
        self.err_m[band]= self.err_m[band][mask]
    
    def nepochs(self,band):
        """
        NAME:
           nepochs
        PURPOSE:
           return the number of epochs
        INPUT:
           band - return the number of epochs in this band
        OUTPUT:
           minimum number of epochs
        HISTORY:
           2010-12-26 - Written - Bovy (NYU)
        """
        nepochs= [len(self.mjd[b]) for b in band]
        return nu.amin(nepochs)
        
    def determine_seasons(self,duration,band,minnepochs=10):
        """
        NAME:
           determine_seasons
        PURPOSE:
           determine the observing seasons
        INPUT:
           duration - duration of the season in yr
           band - do the determination for this band
           minnepochs - minimum number of epochs in a season
        OUTPUT:
           array of MJDs at the start of a season
        HISTORY:
           2012-10-11 - Written - Bovy (IAS)
        """
        #For each data point, see how many data-points are within duration of it
        nepochs= self.nepochs(band)
        ndata_in_duration= sc.zeros(nepochs,dtype='int')
        this_mjd= sc.sort(self.mjd[band])
        for ii in range(nepochs):
            ndata_in_duration[ii]= sc.sum(((self.mjd[band]-this_mjd[ii]) < duration)*((self.mjd[band]-this_mjd[ii]) >= 0.))
            #If there is a gap in the window, exclude it
            if ndata_in_duration[ii] > minnepochs \
                    and sc.any((sc.roll(this_mjd[ii:ii+ndata_in_duration[ii]-1],-1)-this_mjd[ii:ii+ndata_in_duration[ii]-1]) > 0.3):
                ndata_in_duration[ii]= 0.
        #Find the local maxima of this array, bigger than minnepochs
        thisn, ii= ndata_in_duration[0], 0
        max_indx= [0]
        for ii in range(1,nepochs):
            if ndata_in_duration[ii] > thisn:
                max_indx.append(ii)
            thisn= ndata_in_duration[ii]
        cand_mjd= this_mjd[max_indx]
        cand= ndata_in_duration[max_indx]
        return cand_mjd[(cand > minnepochs)]

    def determine_seasons_duration(self,duration,band,minnepochs=10):
        """
        NAME:
           determine_seasons_duration
        PURPOSE:
           determine the observing seasons' duration
        INPUT:
           duration - duration of the season in yr (upper limit)
           band - do the determination for this band
           minnepochs - minimum number of epochs in a season
        OUTPUT:
           array of MJDs at the start of a season
        HISTORY:
           2012-10-11 - Written - Bovy (IAS)
        """
        #For each data point, see how many data-points are within duration of it, and determine the duration
        nepochs= self.nepochs(band)
        ndata_in_duration= sc.zeros(nepochs,dtype='int')
        season_duration= sc.zeros(nepochs)
        this_mjd= sc.sort(self.mjd[band])
        for ii in range(nepochs):
            ndata_in_duration[ii]= sc.sum(((self.mjd[band]-this_mjd[ii]) < duration)*((self.mjd[band]-this_mjd[ii]) >= 0.))
            season_duration[ii]= sc.amax(self.mjd[band][((self.mjd[band]-this_mjd[ii]) < duration)*((self.mjd[band]-this_mjd[ii]) >= 0.)])-this_mjd[ii]
            if ndata_in_duration[ii] > minnepochs \
                    and sc.any((sc.roll(this_mjd[ii:ii+ndata_in_duration[ii]-1],-1)-this_mjd[ii:ii+ndata_in_duration[ii]-1]) > 0.3):
                ndata_in_duration[ii]= 0.
        #Find the local maxima of this array, bigger than minnepochs
        thisn, ii= ndata_in_duration[0], 0
        max_indx= [0]
        for ii in range(1,nepochs):
            if ndata_in_duration[ii] > thisn:
                max_indx.append(ii)
            thisn= ndata_in_duration[ii]
        cand_mjd= this_mjd[max_indx]
        cand_duration= season_duration[max_indx]
        cand= ndata_in_duration[max_indx]
        return cand_duration[(cand > minnepochs)]

    def skew(self,taus,band,minnepochs=15,duration=150./365.25,
             s2=False,s3=False,**kwargs):
        """
        NAME:
           skew
        PURPOSE:
           determine the skew of the quasar lightcurve
        INPUT:
           taus - lags to determine the skew at (in yr)
           band - do the determination for this band
           minnepochs= minimum number of epochs in a season
           s2= if True, return S2 rather than the skew
           s3= if True, return S3 rather than the skew
           +self.fit kwargs (if fit needs to be done)
        OUTPUT:
           skew(taus)
        HISTORY:
           2012-10-11 - Written - Bovy (IAS)
        """
        #First determine the seasons
        start_mjds= self.determine_seasons(duration,band,
                                           minnepochs=minnepochs)
        if len(start_mjds) == 0:
            raise RuntimeError("Object does not have a season with a sufficient number of epochs ...")
        #Then interpolate the seasonal lightcurve using Wiener filter
        hasfit= hasattr(self,'LCparams')
        if not hasfit:
            self.fit(band,**kwargs)
        nwindows= len(start_mjds)
        xs, pm, pv= [], [], []
        for ii in range(nwindows):
            thisxs= sc.arange(start_mjds[ii]-0.025,start_mjds[ii]+taus[-1]+0.025,taus[1]-taus[0])
            thispm, thispv= self.predict(thisxs,band)
            xs.append(thisxs)
            pm.append(thispm)
            pv.append(nu.diag(thispv))
#        return (xs,pm,pv)
        #Calculate skew
        if s2:
            return S2(xs,pm,pv,sc.arange(len(taus)))
        elif s3:
            return S3(xs,pm,pv,sc.arange(len(taus)))
        else:
            return skew(xs,pm,pv,sc.arange(len(taus)))
                    
    def plotSF(self,band='r',nGP=5,nx=201,plotMean=False,**kwargs):
        """
        NAME:
           plotSF
        PURPOSE:
           plot the structure function or color structure function
        INPUT:
           band - band to plot (e.g., 'r' for SF, 'gr' for CSF)
           nGP - number of GPs to plot
           nx - number of gridpoints for theoryLC
           plotMean - also plot the sf of the mean of the GP
        OUTPUT:
        HISTORY:
           2010-12-21 - Started - Bovy (NYU)
        """
        hasfit= hasattr(self,'LCparams')
        if hasfit and not hasattr(self,'lc'):
            self.lc= [TheoryLC(trainSet=self._build_trainset(self.fitband),
                               band=self.fitband,
                               params=self.LCparams,
                               type=self.LCtype,
                               mean=self.LCmean,
                               nx=nx,
                               extramean=[nu.mean(self.m[b]) for b in self.fitband])
                      for ii in range(nGP)]
        else:
            if hasfit: print "Warning: re-using self.lc"
        #Calculate data empirical SF
        sfxs, sf, esf, n= self._calc_sf(band,min=0.01,max=10)
        minsize, maxsize= 5, 15
        minn= nu.nanmin(n)
        maxn= nu.nanmax(n)
        for ii in range(len(n)): 
            if nu.isnan(n[ii]): n[ii]= (nu.nanmin(n)+nu.nanmax(n))/2.
        #Plot it
        pyplot.figure()
        for ii in range(len(sf)):
            if nu.isnan(sf[ii]): continue
            if sf[ii] <= 0.:
                n[ii]= (nu.nanmin(n)+nu.nanmax(n))/2. #Remove from consideratio
                continue
            pyplot.loglog(sfxs[ii],sf[ii],'k.',zorder=5,
                          ms=(minsize+(n[ii]-minn)/(maxn-minn)*(maxsize-minsize)))
        offset_x= nu.exp(nu.log(sfxs[nu.argmin(n)])+0.13)-sfxs[nu.argmin(n)]
        offset_y= nu.exp(nu.log(sf[nu.argmin(n)])+0.13)-sf[nu.argmin(n)]
        bovy_plot.bovy_text(sfxs[nu.argmin(n)]+offset_x,
                            sf[nu.argmin(n)]+offset_y,
                            str(int(minn)))
        offset_x= nu.exp(nu.log(sfxs[nu.argmax(n)])+0.13)-sfxs[nu.argmax(n)]
        offset_y= nu.exp(nu.log(sf[nu.argmax(n)])+0.13)-sf[nu.argmax(n)]
        bovy_plot.bovy_text(sfxs[nu.argmax(n)]+offset_x,
                            sf[nu.argmax(n)]+offset_y,
                            str(int(maxn)))
        if hasfit:
            for lc in self.lc:
                lc.plotSF(band,**kwargs)
        pyplot.errorbar(sfxs,sf,yerr=esf,ecolor='k',fmt=None)
        if hasfit and plotMean:
            pmean, pvar= self.predict(lc.txs,band)
            if len(band) == 1:
                pyplot.loglog(sc.arange(1.,len(pmean)/2)*
                              (lc.txs[1]-lc.txs[0]),
                              2.*sc.var(pmean)
                              -2.*sc.correlate(pmean
                                               -sc.mean(pmean),
                                               pmean
                                           -sc.mean(pmean),"same")\
                              [1:len(pmean)/2][::-1]/\
                                  len(pmean),'-',color='0.65')
            else:
                print "WARNING: multiple bands plotMean not implemented..."
                return None
        #for color-sf, plot data error separately
        if len(band) == 2:
            datacsf= 2.*self.meanErr[band[0]]**2.\
                +2.*self.meanErr[band[1]]**2.,
            pyplot.fill_between(nu.array([0.01,10.]),
                                nu.array([datacsf,datacsf]).reshape(2),
                                y2=nu.array([0.0001,0.0001]),
                                color='0.9',
                                facecolor='0.9')
        pyplot.xlabel(r'$\Delta t\ (\mathrm{lag})\ [\mathrm{yr}]$')
        if self.flux:
            pyplot.ylabel(r'$SF_{'+band+'}\ [\mathrm{dimensionless}]$')
        else:
            pyplot.ylabel(r'$SF_{'+band+'}\ [\mathrm{mag}^2]$')
        pyplot.xlim(0.01,10)
        pyplot.ylim(0.0001,1.)

    def _calc_sf(self,band,min=None,max=None):
        sfxs= []
        sf= []
        esf= []
        n= []
        allmjddiffs= []
        alldiffs= []
        if len(band) == 1:
            for ii in range(self.nepochs(band)):
                for jj in range(ii+1,self.nepochs(band)):
                    allmjddiffs.append(self.mjd[band][ii]-self.mjd[band][jj])
                    alldiffs.append(self.m[band][ii]-self.m[band][jj])
        elif len(band) == 2:
            #match mjds
            xs= []
            ys= []
            for ii in range(len(self.mjd[band[0]])):
                mjd= self.mjd[band[0]][ii]
                match= ((mjd-self.mjd[band[1]])**2. \
                            < (1./365.25/24.)**2)
                if True in match:
                    xs.append(mjd)
                    ys.append(self.m[band[0]][ii]-
                              self.m[band[1]][match])
            xs= nu.array(xs)
            ys= nu.array(ys)
            for ii in range(len(xs)):
                for jj in range(ii+1,len(xs)):
                    allmjddiffs.append(xs[ii]-xs[jj])
                    alldiffs.append(ys[ii]-ys[jj])
        allmjddiffs= nu.array(allmjddiffs)
        alldiffs= nu.array(alldiffs)
        allmjddiffs= nu.fabs(allmjddiffs)
        if min is None: mindiff= nu.amin(allmjddiffs)
        else: mindiff= min
        if max is None: maxdiff= nu.amax(allmjddiffs)
        else: maxdiff= max
        #Bin logarithmically
        nbins= 10
        bins= nu.exp(nu.linspace(nu.log(mindiff),nu.log(maxdiff),nbins+1))
        #bins= nu.append(bins,maxdiff)
        for ii in range(nbins):
            mask= (allmjddiffs >= bins[ii])*(allmjddiffs <= bins[ii+1])
            try:
                thismjddiffs= allmjddiffs[mask]
                thisdiffs= alldiffs[mask]
                sfxs.append(nu.mean(thismjddiffs))
                sf.append(nu.mean(thisdiffs**2.))
                esf.append(nu.sqrt(nu.var(thisdiffs**2.)/len(thisdiffs)))
                if len(thisdiffs) == 0: n.append(nu.nan)
                else: n.append(len(thisdiffs))
            except TypeError:
                sf.append(nu.nan)
                sfxs.append(nu.nan)
                esf.append(nu.nan)
                n.append(nu.nan)
        if len(band) == 1:
            return (nu.array(sfxs),nu.array(sf)-2.*self.meanErr[band]**2.,
                    nu.array(esf),n)
        elif len(band) == 2:
            return (nu.array(sfxs),
                    nu.array(sf),
                    nu.array(esf),n)

    def plot(self,band='r',nGP=5,nx=201,plotMeanVar=False,nsigma=1,
             **kwargs):
        """
        NAME:
           plot
        PURPOSE:
           plot the data
        INPUT:
           band - band to plot (e.g., 'r'; 'ug' plots u-g)
           nGP - number of GPs to plot
           nx - number of gridpoints for theoryLC
           plotMeanVar= if True, plot mean and variance of GP (as snake)
           nsigma= number of sigmas to plot in plotMeanVar (int)
        OUTPUT:
           (none)
        HISTORY:
           2010-12-21 - Started - Bovy (NYU)
        """
        hasfit= hasattr(self,'LCparams')
        if hasfit and not hasattr(self,'lc'):
            self.lc= [TheoryLC(trainSet=self._build_trainset(self.fitband),
                               band=self.fitband,
                               params=self.LCparams,
                               type=self.LCtype,
                               mean= self.LCmean,
                               nx=nx,
                               extramean=[nu.mean(self.m[b]) for b in self.fitband])
                      for ii in range(nGP)]
        else:
            if hasfit: print "Warning: re-using self.lc"
        if len(band) == 1:
            bovy_plot.bovy_plot(self.mjd[band],
                                self.m[band],
                                'k.',zorder=5,ms=10)
            pyplot.errorbar(nu.amax(self.mjd[band])+0.3\
                                -0.1*(nu.amax(self.mjd[band])\
                                          +0.6-nu.amin(self.mjd[band])),
                            nu.amax(self.m[band])+0.1\
                                -0.2*(nu.amax(self.m[band])\
                                          +0.2-nu.amin(self.m[band])),
                            yerr=self.meanErr[band],color='k',
                            elinewidth=3,capsize=5)
            if hasfit:
                for lc in self.lc:
                    lc.plot(band,**kwargs)
            if hasfit and plotMeanVar:
                pmean, pvar= self.predict(lc.txs,band)
                bovy_plot.bovy_plot(lc.txs,pmean,'k--',overplot=True,
                                    lw=2.)
                pyplot.fill_between(lc.txs,pmean-nu.sqrt(nu.diagonal(pvar)),
                                    pmean+nu.sqrt(nu.diagonal(pvar)),
                                    color='.75')
                colord, cc= (1.-0.75)/nsigma, 1.
                while nsigma > 1:
                    pyplot.fill_between(lc.txs,pmean-nsigma*\
                                            nu.sqrt(nu.diagonal(pvar)),
                                        pmean-(nsigma-1.)*\
                                            nu.sqrt(nu.diagonal(pvar)),
                                        color='%f' % (1.-colord*cc))
                    pyplot.fill_between(lc.txs,pmean+nsigma*\
                                            nu.sqrt(nu.diagonal(pvar)),
                                        pmean+(nsigma-1.)*\
                                            nu.sqrt(nu.diagonal(pvar)),
                                        color='%f' % (1.-colord*cc))
                    cc+= 1.
                    nsigma-= 1
            pyplot.xlabel(r'$\mathrm{MJD-constant}\ [\mathrm{yr}]$')
            if self.flux:
                pyplot.ylabel(r'$f_'+band+'/\langle f_{'+band+r'} \rangle$')
            else:
                pyplot.ylabel(r'$'+band+'_0\ [\mathrm{mag}]$')
            pyplot.xlim(nu.amin(self.mjd[band])-0.3,
                        nu.amax(self.mjd[band])+0.3)
            pyplot.ylim(nu.amax(self.m[band])+0.1,
                        nu.amin(self.m[band])-0.1)
            bovy_plot._add_ticks()
        elif len(band) == 2:
            #match mjds
            xs= []
            ys= []
            for ii in range(len(self.mjd[band[0]])):
                mjd= self.mjd[band[0]][ii]
                match= ((mjd-self.mjd[band[1]])**2. \
                            < (1./365.25/24.)**2)
                if True in match:
                    xs.append(mjd)
                    if self.flux:
                        ys.append(self.m[band[0]][ii]/
                                  self.m[band[1]][match])
                    else:
                        ys.append(self.m[band[0]][ii]-
                                  self.m[band[1]][match])
            xs= nu.array(xs)
            ys= nu.array(ys)
            bovy_plot.bovy_plot(xs,ys,
                                'k.',zorder=5,ms=10)
            if self.flux:
                print "Errorbar plotting not implemented yet"
            else:
                pyplot.errorbar(nu.amax(xs)+0.3\
                                    -0.1*(nu.amax(xs)\
                                              +0.6-nu.amin(xs)),
                                nu.amax(ys)+0.1\
                                    -0.2*(nu.amax(ys)\
                                              +0.2-nu.amin(ys)),
                                yerr=nu.sqrt(self.meanErr[band[0]]**2.+
                                             self.meanErr[band[1]]**2.),
                                color='k',
                                elinewidth=3,capsize=5)
            if hasfit:
                for lc in self.lc:
                    lc.plot(band,**kwargs)
            if hasfit and plotMeanVar:
                pmean, pvar= self.predict(lc.txs,band)
                bovy_plot.bovy_plot(lc.txs,pmean,'k--',overplot=True,
                                    lw=2.)
                pyplot.fill_between(lc.txs,pmean-nu.sqrt(nu.diagonal(pvar)),
                                    pmean+nu.sqrt(nu.diagonal(pvar)),
                                    color='.75')
                colord, cc= (1.-0.75)/nsigma, 1.
                while nsigma > 1:
                    pyplot.fill_between(lc.txs,pmean-nsigma*\
                                            nu.sqrt(nu.diagonal(pvar)),
                                        pmean-(nsigma-1.)*\
                                            nu.sqrt(nu.diagonal(pvar)),
                                        color='%f' % (1.-colord*cc))
                    pyplot.fill_between(lc.txs,pmean+nsigma*\
                                            nu.sqrt(nu.diagonal(pvar)),
                                        pmean+(nsigma-1.)*\
                                            nu.sqrt(nu.diagonal(pvar)),
                                        color='%f' % (1.-colord*cc))
                    cc+= 1.
                    nsigma-= 1
            pyplot.xlabel(r'$\mathrm{MJD-constant}\ [\mathrm{yr}]$')
            if self.flux:
                pyplot.ylabel(r'$f_'+band[0]+' / f_'+band[1]+'$')
            else:
                pyplot.ylabel(r'$('+band[0]+' - '+band[1]+')_0\ [\mathrm{mag}]$')
            pyplot.xlim(nu.amin(xs)-0.3,
                        nu.amax(xs)+0.3)
            pyplot.ylim(nu.amax(ys)+0.1,
                        nu.amin(ys)-0.1)
            bovy_plot._add_ticks()
        else:
            raise IOError("'band' must have length 1 or 2")

    def plotColorMag(self,band='r',color='gr',nGP=5,nx=201,
                     plotMean=False,**kwargs):
        """
        NAME:
           plotColorMag
        PURPOSE:
           plot the data in color-magnitude
        INPUT:
           band - band to plot (e.g., 'r')
           color - color to plot (e.g., 'gr')
           nGP - number of GPs to plot
           nx - number of gridpoints for theoryLC
           plotMean - also plot the mean from the GP
        OUTPUT:
           (none)
        HISTORY:
           2011-04-27 - Started - Bovy (NYU)
        """
        hasfit= hasattr(self,'LCparams')
        if hasfit and not hasattr(self,'lc'):
            self.lc= [TheoryLC(trainSet=self._build_trainset(self.fitband),
                               band=self.fitband,
                               params=self.LCparams,
                               type=self.LCtype,
                               mean=self.LCmean,
                               nx=nx,
                               extramean=[nu.mean(self.m[b]) for b in self.fitband])
                      for ii in range(nGP)]
        else:
            print "Warning: re-using self.lc"
        #match mjds
        xs= []
        ys= []
        zs= []
        for ii in range(len(self.mjd[band[0]])):
            mjd= self.mjd[band[0]][ii]
            matchy= ((mjd-self.mjd[color[0]])**2. \
                         < (1./365.25/24.)**2)
            matchz= ((mjd-self.mjd[color[1]])**2. \
                         < (1./365.25/24.)**2)
            if True in matchy and True in matchz:
                xs.append(self.m[band[0]][ii])
                ys.append(self.m[color[0]][matchy])
                zs.append(self.m[color[1]][matchz])
        xs= nu.array(xs)
        ys= nu.array(ys)
        zs= nu.array(zs)
        if self.flux:
            bovy_plot.bovy_plot(xs.flatten(),ys.flatten()/zs.flatten(),
                                'k.',zorder=5,ms=10)
        else:
            bovy_plot.bovy_plot(xs.flatten(),ys.flatten()-zs.flatten(),
                                'k.',zorder=5,ms=10)
        if self.flux:
            print "Errorbar plotting not implemented yet"
        else:
            pyplot.errorbar(nu.amax(xs)+0.3\
                                -0.1*(nu.amax(xs)\
                                          +0.6-nu.amin(xs)),
                            nu.amax(ys-zs)+0.1\
                                -0.2*(nu.amax(ys-zs)\
                                          +0.2-nu.amin(ys-zs)),
                            yerr=nu.sqrt(self.meanErr[color[0]]**2.+
                                         self.meanErr[color[1]]**2.),
                            xerr=self.meanErr[band[0]],
                            color='k',
                            elinewidth=3,capsize=5)
        if hasfit:
            for lc in self.lc:
                lc.plotColorMag(band,color,**kwargs)
            if hasfit and plotMean:
                exs= [(ti,band) for ti in lc.txs]
                pmean, pvar= self.predict(exs,band)
                pmeanc, pvarc= self.predict(lc.txs,color)
                bovy_plot.bovy_plot(pmean,pmeanc,'k--',overplot=True,
                                    lw=2.)
        pyplot.xlabel(r'$'+band[0]+'_0\ [\mathrm{mag}]$')
        if self.flux:
            pyplot.ylabel(r'$'+color[0]+' \ '+color[1]+'$')
        else:
            pyplot.ylabel(r'$('+color[0]+' - '+color[1]+')_0\ [\mathrm{mag}]$')
        pyplot.xlim(nu.amax(xs)+0.3,
                    nu.amin(xs)-0.3)
        if self.flux:
            pyplot.ylim(nu.amin(ys-zs)-0.1,
                        nu.amax(ys-zs)+0.1)
        else:
            pyplot.ylim(nu.amin(ys/zs)-0.1,
                        nu.amax(ys/zs)+0.1)
        bovy_plot._add_ticks()
        
    def mjd_overlap(self,band='gr'):
        """
        NAME:
           mjd_overlap
        PURPOSE:
           return the index with matching mjd's (dictionary)
        INPUT:
           band - bands to find the overlap between
        OUTPUT:
           dictionary of indices
        HISTORY:
           2011-06-02 - Written - Bovy (NYU)
        """
        if len(band) < 2: return None
        depoch= 71.7/60./60./24./365.25 #dt between consecutive bands (from http://www.sdss3.org/instruments/camera.php)
        #initialize:
        out= {}
        for b in band:
            out[b]= []
        #go through all of the mjds
        for ii in range(len(self.mjd[band[0]])):
            mjd= self.mjd[band[0]][ii]
            trialindices= [ii]
            for b in band[1:len(band)]:
                thismjds= nu.fabs(self.mjd[b]-mjd)
                indx= (thismjds < depoch*5.) #5 should be enough
                if True in indx:
                    trialindices.append(nu.arange(len(self.mjd[b]))[indx][0])
            if len(trialindices) == len(band):
                for jj in range(len(band)):
                    out[band[jj]].append(trialindices[jj])
        return out

    def resample(self,xs,band='ugriz',errors=True,noconstraints=False,
                 wedge=False,wedgerate=0.25*365.25,wedgetau=200./365.25,
                 addlagged=False,lag=None,laggedxs=None):
        """
        NAME:
           resample
        PURPOSE:
           resample a lightcurve
        INPUT:
           xs - xs to resample (includes band info if relevant); a subset
                of these are taken that are in band
           band - filters to sample
           noconstraints - don't constrain the resampled lightcurve using the data
           wedge= if True, sample from wedge process (default: False)
           wedgerate= rate of wedges
           wedgetau= tau of wedges
           ##Following only work for single band
           addlagged= False: if True, add a lagged version
           lag= lag in yr if addlagged
           laggedxs= xs for lagged curve (pre-lag, i.e., observer's frame)
        OUTPUT:
           another VarQso object
        HISTORY:
           2010-12-30 - Written - Bovy (NYU)
        """
        #First filter the xs
        if len(band) == 1 and not isinstance(xs[0],(int,float)) and \
                len(xs[0]) == 2:
            if not band in self.fitband: 
                print "Error: requested band is not the fitted band"
                return None
            if len(self.fitband) == 1:
                xs= [x[0] for x in xs if x[1] == band]
            else:
                xs= [x for x in xs if x[1] == band]
            if xs == []: 
                print "Error: no xs corresponding to the requested band"
                return None
        elif len(band) == 1 and ((not isinstance(xs[0],(int,float)) and 
                                 len(xs[0]) == 1) or \
                                     isinstance(xs[0],(int,float))):
            pass
        elif len(band) > 1 and ((not isinstance(xs[0],(int,float)) and 
                                 len(xs[0]) == 1) or \
                                    isinstance(xs[0],(int,float))):
            xs= [(x,b) for x in xs for b in band]
        elif len(band) > 1 and not isinstance(xs[0],(int,float)) and \
                len(xs[0]) == 2:
            xs= [x for x in xs if x[1] in band]
        if addlagged:
            xs= list(xs)
            xs.extend(list(laggedxs-lag))
            xs= nu.array(xs)
        #Now get ready for drawing
        cf= self.LC.cf
        mf= self.LC.mf
        tiny_cholesky=10.**-4.*covar_func(0.,0.,(cf))
        mean={}
        for b in band:
            mean[b]= nu.mean(self.m[b])
        #Now draw
        try:
            if wedge:
                #Perform wedge sampling
                start= nu.amin(xs)-wedgetau
                dts= nu.random.exponential(scale=1./wedgerate,size=10000)#to be sure
                times= nu.cumsum(dts)+start
                #trim times
                times= times[(times < nu.amax(xs))]
                dataSF= 2.*(covar_func(0.,0.,(cf))-covar_func(0.,0.025,(cf)))
                amp= nu.sqrt(dataSF/wedgerate/0.025) #V=SF at 
                if not self.flux:
                    amp*= -1.
                GPsample= nu.zeros(xs.shape)
                for ii in range(len(times)):
                    GPsample+= wedge_func(xs,times[ii],amp,tau=wedgetau)
                GPsample-= nu.mean(GPsample)
            elif noconstraints:
                GPsample= eval_gp(xs,mean_func,covar_func,(mf),(cf),nGP=1,
                                  constraints=None,
                                  tiny_cholesky=tiny_cholesky).reshape(len(xs))
            else:
                GPsample= eval_gp(xs,mean_func,covar_func,(mf),(cf),nGP=1,
                                  constraints=self._build_trainset(self.fitband),
                                  tiny_cholesky=tiny_cholesky).reshape(len(xs))
        except nu.linalg.linalg.LinAlgError:
            raise
            #ACTUALLY THIS NEVER HAPPENS AND THE CODE BLOCK BELOW IS UNTESTED
            #Fails are because MJDs are so close together that they are perfectly
            #correlated (e.g., SDSS ugriz observations)
            #Remove the offending mjds and fill them in later, typically observations
            #spaced less than a day apart are completely correlated
            print "Warning: constrained covariance matrix was not positive-definite"
            print "Trying to deal with it ..."
            mjds= [x[0] for x in xs]
            sortindx= sorted(range(len(mjds)), key = mjds.__getitem__)
            #print mjds
            #print sortindx
            mjds= [mjds[s] for s in sortindx]
            xs= [xs[s] for s in sortindx]
            groups= []
            ii= 0
            thisxs= xs
            while ii < len(mjds):
                cmp_mjds= nu.array([m-mjds[ii] for m in mjds])
                cmp_mjds[ii]= -1.
                indx= (nu.fabs(cmp_mjds) < 1./365.25)
                corr_mjds= nu.array(mjds)[indx]
                try:
                    if len(corr_mjds) > 0:
                        addn= len(corr_mjds)
                        for jj in range(addn):
                            thisxs.pop(jj+1)
                        groups.append((ii,addn))
                    else: addn= 0
                except TypeError:
                    addn= 0
                ii+= 1+addn
            #print thisxs
            GPsample= eval_gp(thisxs,mean_func,covar_func,(mf),(cf),nGP=1,
                              constraints=self._build_trainset(self.fitband),
                              tiny_cholesky=0.00001).reshape(len(xs))
            GPsample= list(GPsample)
            #Fill in the blanks
            for g in groups:
                for jj in range(g[1]):
                    GPsample.insert(g[0]+1,GPsample[g[0]])
            GPsample= nu.array(GPsample)
            #print GPsample
        #Add the mean
        if len(band) == 1 and ((not isinstance(xs[0],(int,float)) and 
                                len(xs[0]) == 1) or \
                                   isinstance(xs[0],(int,float))):
            GPsample+= mean[band]
        else:
            for b in band:
                for ii in range(len(xs)):
                    if xs[ii][1] == b:
                        GPsample[ii]+= mean[b]
        #Add errors
        if errors:
            errs= []
            for ii in range(len(xs)):
                if len(band) == 1:
                    errs.append(self.meanErr[band])
                    GPsample[ii]+= nu.random.randn()*self.meanErr[band]
                else:
                    errs.append(self.meanErr[xs[ii][1]])
                    GPsample[ii]+= nu.random.randn()*self.meanErr[xs[ii][1]]
        else:
            errs= nu.zeros(len(GPsample))
        #Setup VarQso object
        if addlagged:
            #Split xs
            xs= xs[0:len(xs)-len(laggedxs)]
            vA= VarQso(xs,GPsample[0:len(GPsample)-len(laggedxs)],
                       errs[0:len(GPsample)-len(laggedxs)],
                       band=band,medianize=False)
            vB= VarQso(laggedxs,GPsample[len(GPsample)-len(laggedxs):len(GPsample)],
                       errs[len(GPsample)-len(laggedxs):len(GPsample)],
                       band=band,medianize=False)
            return (vA,vB)
        else:
            return VarQso(xs,GPsample,errs,band=band,medianize=False)

    def fit(self,band,type='powerlawSF',loglike=False,mean='zero',
            init_params=None,fix=None):
        """
        NAME:
           fit
        PURPOSE:
           fit a structure function model to the data
        INPUT:
           band - band to fit
           type - type of structure function model to fit
           loglike - return loglike in the parameters dict
           mean - type of mean fitting ('zero','const',...)
           init_params - initial parameters for fit
           fix= None or list of parameters to hold fixed
        OUTPUT:
           best-fit parameters dictionary
        HISTORY:
           2010-12-21 - Written - Bovy (NYU)
        """
        trainSet= self._build_trainset(band)
        self.fitband= band
        self.LC= LCmodel(trainSet=trainSet,type=type,mean=mean,
                         init_params=init_params)
        self.LCparams= self.LC.fit(fix=fix)
        self.LCtype= type
        self.LCmean= mean
        if loglike:
            (params,packing)= pack_params(self.LC.cf,self.LC.mf)
            covarFuncName= inspect.getmodule(self.LC.cf).__name__
            thisCovarClass= __import__(covarFuncName)
            loglike= marginalLikelihood(params,trainSet,packing,
                                        thisCovarClass,None,False)
            self.LCparams['loglike']= -loglike
        return self.LCparams

    def sampleGP(self,band=None,type='powerlawSF',mean='zero',
                 nsamples=100,step=None,fix=None,metropolis=False,
                 markovpy=False,burnin=20):
        """
        NAME:
           sampleGP
        PURPOSE:
           sample the structure function model parameters
        INPUT:
           band - band to sample
           type= type of structure function model to sample
           mean= type of mean fitting ('zero','const',...)
           nsamples= number of samples desired
           fix= None or list of parameters to hold fixed
           metropolis= if True, use Metropolis sampling
           markovpy= if True, use markovpy sampling
           burnin= number of samples to discard as burn-in
        OUTPUT:
           (none; get the samples using self.get_sampleGP())
        HISTORY:
           2011-06-13 - Written - Bovy (NYU)
        """
        hasfit= hasattr(self,'LCparams')
        if not hasfit: #build trainingset and optimize before sampling
            if band is None:
                raise IOError("'band' must be set")
            trainSet= self._build_trainset(band)
            self.fitband= band
            self.LC= LCmodel(trainSet=trainSet,type=type,mean=mean,
                             init_params=self.LCparams)
            self.LCparams= self.LC.fit(fix=fix)
            self.LCtype= type
            self.LCmean= mean
        #Step defaults
        if step is None and not markovpy:
            if self.LCtype == 'powerlawSF':
                if len(self.fitband) > 1:
                    step= [0.3,0.3,0.4,0.4]
                else:
                    step= [0.3,0.4]
            elif self.LCtype == 'periodicDRW':
                step= [0.3,0.2,0.2] #l,p,a
            elif self.LCtype == 'KS11':
                step= [0.1,0.3,0.4]
            elif self.LCtype == 'DRW':
                step= [1.,1.]
            elif self.LCtype == 'reverbPLgr':
                step= 0.3
        elif step is None:
            step= 0.01
        out= self.LC.sampleGP(nsamples=nsamples+burnin,
                              step=step,fix=fix,
                              metropolis=metropolis,markovpy=markovpy)
        self.LCparamsSamples= out[burnin:len(out)]

    def get_sampleGP(self):
        """
        NAME:
           get_sampleGP
        PURPOSE:
           return the samples from sampling the GP hyperparameters
        INPUT:
           (none)
        OUTPUT:
           list of dictionaries of parameters
        HISTORY:
           2011-06-13 - Written - Bovy (NYU)
        """
        if not hasattr(self,'LCparamsSamples'):
            raise AttributeError("You must sample first using self.sampleGP")
        return self.LCparamsSamples

    def set_sampleGP(self,samples):
        """
        NAME:
           set_sampleGP
        PURPOSE:
           set the samples from sampling the GP hyperparameters
        INPUT:
           samples - samples dictionary
        OUTPUT:
           (none)
        HISTORY:
           2011-06-18 - Written - Bovy (NYU)
        """
        self.LCparamsSamples= samples

    def plot_sampleGP(self,d1='logA',d2='gamma',xrange=None,yrange=None,
                      xlabel=None,ylabel=None,bins=None,logx=False,logy=False,
                      whitemax=False):
        """
        NAME:
           plot_sampleGP
        PURPOSE:
           plot the output from sampleGP
        INPUT:
           d1= parameter to plot on the x-axis (e.g., d1='logA')
           d2= parameter to plot on the y-axis (e.g., d2='gamma')
           xrange=, yrange=
           xlabel=, ylabel=
           bins= number of bins to use (in one dim)
           logx= logy= if True, take the log of the x/y input
           whitemax= if True, use a white cross to indicate the maximum
        OUTPUT:
           plot to output device
        HISTORY:
           2011-06-13- Written - Boyv (NYU)
        """
        if not hasattr(self,'LCparamsSamples'):
            raise AttributeError("You must sample first using self.sampleGP")
        xs, ys= [], []        
        for sample in self.LCparamsSamples:
            if d1 == 'logA': xs.append(sample[d1][0]/2.)
            elif d1 == 's': xs.append(sample[d1][0]-1.)
            else: xs.append(sample[d1][0])
            if d2 == 'logA': ys.append(sample[d2][0]/2.)
            elif d2 == 's': ys.append(sample[d2][0]-1.)
            else: ys.append(sample[d2][0])
        xlabelnounits= self._plot_sampleGP_label(xlabel,d1,nounits=True)
        ylabelnounits= self._plot_sampleGP_label(ylabel,d2,nounits=True)
        xlabel= self._plot_sampleGP_label(xlabel,d1)
        ylabel= self._plot_sampleGP_label(ylabel,d2)
        xs= nu.array(xs)
        ys= nu.array(ys)
        if logx: xs= nu.log(xs)
        if logy: ys= nu.log(ys)
        bovy_plot.scatterplot(xs,ys,
                              'k,',onedhists=True,
                              bins= self._plot_sampleGP_bins(bins,len(xs)),
                              xlabel=xlabel,ylabel=ylabel,
                              xrange=self._plot_sampleGP_range(xrange,d1,logx),
                              yrange=self._plot_sampleGP_range(yrange,d2,logy))
        #Overplot the best fit
        if d1 == 'logA': xb= self.LCparams[d1][0]/2.
        elif d1 == 's': xb= self.LCparams[d1][0]-1.
        else: xb= self.LCparams[d1][0]
        if d2 == 'logA': yb= self.LCparams[d2][0]/2.
        elif d2 == 's': yb= self.LCparams[d2][0]-1.
        else: yb= self.LCparams[d2][0]
        if logx: xb= nu.log(xb)
        if logy: yb= nu.log(yb)
        if whitemax:
            bovy_plot.bovy_plot(xb,yb,'wx',overplot=True,ms=10.,mew=2.)
        else:
            bovy_plot.bovy_plot(xb,yb,'rx',overplot=True,ms=10.,mew=2.)
        #Add statistics about the distribution
        d1mean= nu.mean(xs)
        d1std= nu.std(xs)
        d2mean= nu.mean(ys)
        d2std= nu.std(ys)
        corr= nu.corrcoef(xs,ys)[0,1]
        bovy_plot.bovy_text(r'$\langle$'+xlabelnounits+
                            r'$\rangle= %5.2f \pm %5.2f$' % (d1mean,d1std)
                            +'\n'+
                            r'$\langle $'+ylabelnounits+
                            r'$\rangle= %5.2f \pm %5.2f$' % (d2mean,d2std)
                            +'\n'+
                            r'$\rho = %5.2f$' % corr,
                            top_left=True)
        return None

    def _plot_sampleGP_bins(self,bins,n):
        if not bins is None: return bins
        guess= int(nu.sqrt(n))
        return guess

    def _plot_sampleGP_label(self,label,d1,nounits=False):
        if not label is None: return label
        if d1 == 'logA': 
            if nounits:
                return r'$\log A$'
            else:
                return r'$\log A\ (\mathrm{amplitude\ at\ 1\ yr})$'
        elif d1 == 'logP': 
            if nounits:
                return r'$\log P$'
            else:
                return r'$\log P\ (\mathrm{period})$'
        elif d1 == 'gamma':
            if nounits:
                return r'$\gamma$'
            else:
                return r'$\gamma\ (\mathrm{power-law\ exponent})$'
        elif d1 == 'loga2': 
            if nounits:
                return r'$\log a^2$'
            else:
                return r'$\log a^2\ (\mathrm{variance})$'
        elif d1 == 'logl': 
            if nounits:
                return r'$\log \tau$'
            else:
                return r'$\log \tau\ (\mathrm{damping\ timescale})$'
        else: return d1

    def _plot_sampleGP_range(self,range,d1,logx):
        """Internal function to process the axis range for plot_sampleGP"""
        if not range is None: return range
        if d1 == 'logA': out= [-9.21/2.,0.]
        elif d1 == 'gamma':
            if logx:
                out= [-5.,1.]
            else:
                out= [0.,1.2]
        elif d1 == 's': out= [-1.1,.4]
        elif d1 == 'logl': out= [-10.,4]
        elif d1 == 'loga2': out= [-2.5,0.]
        else: return None
        return out
    
    def predict(self,t,band='g'):
        """
        NAME:
           predict
        PURPOSE:
           predict the value at a given time
        INPUT:
           t - time
           band= band to predict (default=g)
        OUTPUT:
           (prediction,variance)
        HISTORY:
           2011-06-03 - Written - Bovy (NYU)
        """
        #Predition is mean and variance of GP
        hasfit= hasattr(self,'LCparams')
        if not hasfit:
            print "Must have set the LC model ..."
            print "Returning ..."
            return
        if not isinstance(t,(list,nu.ndarray)): tlist= False
        else: tlist= True
        cf= self.LC.cf
        mf= self.LC.mf
        if len(band) > 2:
            print "Currently only one band is supported ..."
            print "Returning ..."
            return
        elif len(band) == 2:
            mean= nu.mean(self.m[band[0]])-nu.mean(self.m[band[1]])
            if tlist:
                nt= len(t)
                xs= [(ti,band[0]) for ti in t]
                xs.extend([(ti,band[1]) for ti in t])
            else:
                nt= 1
                xs= [(t,band[0]),(t,band[1])]
        else:
            mean= nu.mean(self.m[band])
            if tlist:
                xs= t
            else:
                xs= [t]
        pmean= calc_constrained_mean(xs,mf,(),cf,(),
                                     self._build_trainset(self.fitband))
        pvar= calc_constrained_covar(xs,cf,(),
                                     self._build_trainset(self.fitband))
        if len(band) == 2:
            pmean= [pmean[ii]-pmean[ii+nt] for ii in range(nt)]
            pvarout= nu.zeros((nt,nt))
            for ii in range(nt):
                for jj in range(nt):
                    pvarout[ii,jj]= pvar[ii,jj]+pvar[ii+nt,jj+nt]\
                        -pvar[ii+nt,jj]-pvar[ii,jj+nt]
            pvar= pvarout
        return (mean+pmean,pvar)

    def loglike(self,band,type=None,params=None,mean='zero'):
        """
        NAME:
           loglike
        PURPOSE:
           calculate the log likelihood of this lightcurve given a model 
           for the lightcurve SF/covariance
        INPUT:
           band
           type - covar/SF type
           params - parameters of the covariance/SF
        OUTPUT:
           log likelihood
        HISTORY:
           2011-01-09 - Written - Bovy (NYU)
        """
        if type is None:
            type= self.LCtype
        if params is None:
            params= self.LCparams
        if mean is None:
            mean= self.LCmean
        trainSet= self._build_trainset(band)
        LC= LCmodel(trainSet=trainSet,type=type,init_params=params,
                    mean=mean)
        (params,packing)= pack_params(LC.cf,LC.mf,None)
        covarFuncName= inspect.getmodule(LC.cf).__name__
        thisCovarClass= __import__(covarFuncName)
        return -marginalLikelihood(params,trainSet,packing,
                                   thisCovarClass,None,False)
        
    def cross_loglike(self,lag,vB,band,type=None,params=None,mean='zero'):
        """
        NAME:
           cross_loglike
        PURPOSE:
           calculate the log likelihood of this lightcurve and a lagged 
           second light curve given a model 
           for the lightcurve SF/covariance
        INPUT:
           lag - a lag
           vB - another VarQso object
           band
           type - covar/SF type
           params - parameters of the covariance/SF
        OUTPUT:
           log likelihood
        HISTORY:
           2011-01-09 - Written - Bovy (NYU)
        """
        if type is None:
            type= self.LCtype
        if params is None:
            params= self.LCparams
        if mean is None:
            mean= self.LCmean
        listx=list(self.mjd[band])
        listy=list(self.m[band]-nu.mean(self.m[band]))
        noise=list(self.err_m[band])
        listx.extend(list(vB.mjd[band]-lag))
        listy.extend(list(vB.m[band]-nu.mean(vB.m[band])))
        noise.extend(list(vB.err_m[band]))
        listx= nu.array(listx)
        listy= nu.array(listy)
        noise= nu.array(noise)
        trainSet= trainingSet(listx=listx,
                              listy=listy,
                              noise=noise)
        LC= LCmodel(trainSet=trainSet,type=type,init_params=params,
                    mean=mean)
        (params,packing)= pack_params(LC.cf,LC.mf,None)
        covarFuncName= inspect.getmodule(LC.cf).__name__
        thisCovarClass= __import__(covarFuncName)
        return -marginalLikelihood(params,trainSet,packing,
                                   thisCovarClass,None,False)
        
    def fit_lag(self,vB,band,type='powerlawSF',mean='zero',
                init_params=None,fix=None,
                init_lag=1.):
        """
        NAME:
           fit
        PURPOSE:
           fit a lag of a second quasar wrt this primary quasar; 
           structure function model is first fit to the data of this quasar
           alone, then a lag is fit
        INPUT:
           vB - second light curve (VarQso object)
           band - band to fit
           type - type of structure function model to fit
           loglike - return loglike in the parameters dict
           mean - type of mean fitting ('zero','const',...)
           init_params - initial parameters for fit
           fix= None or list of parameters to hold fixed
           init_lag= initial lag for the fit
        OUTPUT:
           lag in yr
        HISTORY:
           2010-12-21 - Written - Bovy (NYU)
        """
        hasfit= hasattr(self,'LCparams')
        if not hasfit: #build trainingset and optimize before sampling
            if band is None:
                raise IOError("'band' must be set")
            trainSet= self._build_trainset(band)
            self.fitband= band
            self.LC= LCmodel(trainSet=trainSet,type=type,mean=mean,
                             init_params=self.LCparams)
            self.LCparams= self.LC.fit(fix=fix)
            self.LCtype= type
            self.LCmean= mean
        lag= optimize.fmin_powell((lambda x: -self.cross_loglike(x,vB,band)),
                                  init_lag,args=())
        return lag

    def plot_like(self,band,type=None,params=None,mean='zero',
                  xrange=None,yrange=None,**kwargs):
        """
        NAME:
           plot_like
        PURPOSE:
           plot the  likelihood of this lightcurve given a model 
           for the lightcurve SF/covariance
        INPUT:
           band
           type - covar/SF type
           +bovy_plot.bovy_dens2d kwargs
        OUTPUT:
           log likelihood
        HISTORY:
           2011-01-09 - Written - Bovy (NYU)
        """
        if type is None:
            type= self.LCtype
        if mean is None:
            if not hasattr(self,'LCmean'):
                mean= 'zero'
            mean= self.LCmean
        if type == 'powerlawSF' and xrange is None:
            xrange= [-9.21,0.]
        if type == 'powerlawSF' and yrange is None:
            yrange= [0.,1.25]
        if hasattr(self,'binnedLike'):
            print "Warning: re-using previous binned likelihood"
        else:
            trainSet= self._build_trainset(band)
            LC= LCmodel(trainSet=trainSet,type=type,init_params=params,
                        mean=mean)
            (params,packing)= pack_params(LC.cf,LC.mf,None)
            covarFuncName= inspect.getmodule(LC.cf).__name__
            thisCovarClass= __import__(covarFuncName)
            #Now loop
            binnedLike= nu.zeros((51,51))
            xs= nu.linspace(xrange[0],xrange[1],binnedLike.shape[0])
            ys= nu.linspace(yrange[0],yrange[1],binnedLike.shape[1])
            for ii in range(binnedLike.shape[0]):
                for jj in range(binnedLike.shape[1]):
                    sys.stdout.write('\r'+"Working on %i / %i ...\r" %(ii*binnedLike.shape[1]+jj+1,binnedLike.shape[0]*binnedLike.shape[1]))
                    sys.stdout.flush()
                    if type == 'powerlawSF':
                        LC.cf._dict= {'logA':xs[ii],'gamma':ys[jj]}
                        (params,packing)= pack_params(LC.cf,LC.mf,None)
                        binnedLike[ii,jj]= -marginalLikelihood(params,trainSet,packing,
                                                               thisCovarClass,None,False)
            self.binnedLike= nu.exp(binnedLike)
            self.binnedLike/= nu.sum(self.binnedLike)*(xs[1]-xs[0])*(ys[1]-ys[0])
            sys.stdout.write('\r'+_ERASESTR+'\r')
            sys.stdout.flush()
        bovy_plot.bovy_dens2d(self.binnedLike.T,origin='lower',cmap='gist_yarg',
                              xrange=xrange,yrange=yrange,interpolation='nearest',
                              contours=True,levels=[0.68,0.95,0.99],cntrmass=True,**kwargs)

    def setLCmodel(self,params,band,type='powerlawSF',mean='zero'):
        """
        NAME:
           setLCmodel
        PURPOSE:
           set the lightcurve model for this data set
        INPUT:
           params - parameters of the model
           band - band the model is for
           type - type of lightcurve model
        OUTPUT:
           (none)
        HISTORY:
           2010-12-24 - Written - Bovy (NYU)
        """
        trainSet= self._build_trainset(band)
        self.LC= LCmodel(trainSet=trainSet,type=type,mean=mean,
                         init_params=params)
        self.LCparams= params
        self.LCtype= type
        self.LCmean= mean
        self.fitband= band

    def resetLCmodel(self):
        """
        NAME:
           resetLCmodel
        PURPOSE:
           reset self.lc, the theoretical LC model
        INPUT:
           (none)
        OUTPUT:
           (none)
        HISTORY:
           2010-12-24 - Written - Bovy (NYU)
        """
        self.__dict__.pop('lc')

    def _build_trainset(self,band):
        if len(band) == 1:
            trainSet= trainingSet(listx=self.mjd[band],
                                  listy=self.m[band]-nu.mean(self.m[band]),
                                  noise=self.err_m[band])
        else:
            listx= []
            listy= []
            noise= []
            for b in band:
                listx.extend([(t,b) for t in self.mjd[b]])
                listy.extend([m-nu.mean(self.m[b]) for m in self.m[b]])
                noise.extend([m for m in self.err_m[b]])
            trainSet= trainingSet(listx=listx,listy=listy,noise=noise)
        return trainSet

    def __getstate__(self):
        """Pickling routine"""
        outDict={}
        outDict['mjd']= self.mjd
        outDict['m']= self.m
        outDict['err_m']= self.err_m
        outDict['meanErr']= self.meanErr
        return outDict

    def __setstate__(self,state):
        self.mjd= state['mjd']
        self.m= state['m']
        self.err_m= state['err_m']
        self.meanErr= state['meanErr']
        return None

class LCmodel():
    """Class that represents a light-curve model"""
    def __init__(self,trainSet=None,type='powerlawSF',init_params=None,
                 mean='zero'):
        """
        NAME:
           __init__
        PURPOSE:
           initialize
        INPUT:
           trainSet - a GP training set
           type - type of covariance function to use
           init_params - parameters to initialize the covar func and mean func
           mean= type of mean fitting to do ('zero','const')
        OUTPUT:
           (none)
        HISTORY:
           2010-12-21 - Written - Bovy (NYU)
        """
        self.trainSet= trainSet
        try:
            nbands= len(trainSet.listx[0])
        except TypeError:
            nbands= 1
        self.nbands= nbands
        if self.nbands == 1:
            if type == 'powerlawSF':
                try:
                    from flexgp.powerlawSF import covarFunc 
                except ImportError:
                    from gp.powerlawSF import covarFunc 
                if init_params is None:
                    params= {'gamma': nu.array([ 0.49500723]), 
                             'logA': nu.array([-3.36044037])}
            elif type == 'brokenpowerlawSF':
                from brokenpowerlawSF import covarFunc 
                if init_params is None:
                    params= {'gamma1': nu.array([ 0.8]),
                             'gamma2': nu.array([ 0.4]),  
                             'breakt': nu.array([ 0.2]),
                             'logA': nu.array([-3.36044037])}
            elif type == 'scatter':
                from scatterCovariance import covarFunc 
                if init_params is None:
                    params= {'logA': nu.array([-3.36044037])}
            elif type == 'DRW':
                from DRW import covarFunc
                if init_params is None:
                    params= {'logl': nu.array([-1.37742591]), 
                             'loga2': nu.array([-3.47341754])}
            elif type == 'zero':
                from zeroCovariance import covarFunc
                if init_params is None:
                    params= {}
            elif type == 'periodicDRW':
                try:
                    from flexgp.periodicOU_ARD import covarFunc 
                except ImportError:
                    from gp.periodicOU_ARD import covarFunc  
                if init_params is None:
                    params= {'logl': nu.array([-1.37742591]), 
                             'logP':nu.array([nu.log(1./365.25)]),
                             'loga2': nu.array([-3.47341754])}
            if mean == 'zero':
                try:
                    from flexgp.zeroMean import meanFunc
                except ImportError:
                    from gp.zeroMean import meanFunc
                mean_params= {}
            elif mean == 'const':
                try:
                    from flexgp.constMean import meanFunc
                except ImportError:
                    from gp.constMean import meanFunc
                if init_params is None:
                    mean_params= {'m':0.}
        else: #multi-filter
            if type == 'powerlawSF':
                from powerlawSFmulti import covarFunc
                if init_params is None:
                    params= {'logA': nu.array([-3.59009776]), 
                             'logAgr': nu.array([-5.6487848]), 
                             'gamma': nu.array([ 0.45918053]), 
                             'gammagr': nu.array([ 0.21333858])}
            elif type == 'powerlawSFratios':
                from powerlawSFmultiratios import covarFunc
                if init_params is None:
                    params= {'logAr': nu.array([-6.59009776]), 
                             'logAri': nu.array([-10.6487848]), 
                             'gamma': nu.array([ 0.0305918053]), 
                             'gammagr': nu.array([ 0.0121333858])}
            elif type == 'zero':
                from zeroCovariance import covarFunc
                if init_params is None:
                    params= {}
            elif type == 'DRW':
                from OUgr import covarFunc
                if init_params is None:
                    params= {'logl': nu.array([-1.38968195]), 
                             'loglgr': nu.array([-2.46684501]), 
                             'logagr2': nu.array([-6.62320832]), 
                             'loga2': nu.array([-3.52099305])}
            elif type == 'KS11':
                from KS11 import covarFunc
                if init_params is None:
                    params= {'logA': nu.array([-3.59009776]), 
                             'gamma': nu.array([ 0.45918053]), 
                             's':nu.array([0.83])}
            elif type == 'KS11_nogamma':
                from KS11_nogamma import covarFunc
                if init_params is None:
                    params= {'logA': nu.array([-3.59009776]), 
                             's':nu.array([0.83])}
            elif type == 'KS11constscatter':
                from KS11constscatter import covarFunc
                if init_params is None:
                    params= {'logA': nu.array([-3.59009776]), 
                             'gamma': nu.array([ 0.45918053]), 
                             's':nu.array([0.83]),
                             'logsigma':nu.array([-3.])}
            elif type == 'KS11multiconstscatter':
                from KS11multiconstscatter import covarFunc
                if init_params is None:
                    params= {'logA': nu.array([-3.59009776]), 
                             'gamma': nu.array([ 0.45918053]), 
                             'sr':nu.array([0.83]),
                             'si':nu.array([1.]),
                             'sz':nu.array([1.]),
                             'logsigma':nu.array([-3.])}
            elif type == 'KS11multi':
                from KS11multi import covarFunc
                if init_params is None:
                    params= {'logA': nu.array([-3.59009776]), 
                             'gamma': nu.array([ 0.45918053]), 
                             'sr':nu.array([0.83]),
                             'si':nu.array([1.]),
                             'sz':nu.array([1.])}
            elif type == 'reverbPLgr':
                from reverbPLgr import covarFunc
                if init_params is None:
                    params= {'logA': nu.array([-5.59009776]), 
                             'gamma': nu.array([ 0.45918053]), 
                             'logtL':nu.array([0.]),
                             'logsigmatL':nu.array([-0.5]),
                             'beta':nu.array([0.3]),
                             'eps':nu.array([0.2]),
                             's':nu.array([0.83])}
            if mean == 'zero':
                try:
                    from flexgp.zeroMean import meanFunc
                except ImportError:
                    from gp.zeroMean import meanFunc
                mean_params= {}
            elif mean == 'const' \
                    and (type == 'KS11' or type == 'KS11constscatter' \
                             or type == 'KS11_nogamma' \
                             or type == 'reverbPLgr'):
                from KS11 import meanFunc
                if init_params is None:
                    mean_params= {'m':nu.array([0.]),
                                  'b':nu.array([0.]),
                                  's':nu.array([0.83])}
            elif mean == 'const' \
                    and (type == 'KS11multiconstscatter' \
                             or type == 'KS11multi'):
                from KS11multiconstscatter import meanFunc
                if init_params is None:
                    mean_params= {'m':nu.array([0.]),
                                  'br':nu.array([0.]),
                                  'bi':nu.array([0.]),
                                  'bz':nu.array([0.]),
                                  'sr':nu.array([0.83]),
                                  'si':nu.array([1.]),
                                  'sz':nu.array([1.])}
        if not init_params is None:
            params= init_params
            mean_params= init_params
        self.cf= covarFunc(**params)
        self.mf= meanFunc(**mean_params)

    def fit(self,fix=None):
        """
        NAME:
           fit
        PURPOSE:
           fit a LC model
        INPUT:
           fix= None or list of parameters to hold fixed
        OUTPUT:
           dictionary of best-fit parameters
        HISTORY:
           2010-12-21 - Written - Bovy (NYU)
        """
        (self.cf,self.mf)= trainGP(self.trainSet,self.cf,useDerivs=False,
                                     mean=self.mf,fix=fix)
        return dict(self.cf._dict,**self.mf._dict)

    def sampleGP(self,nsamples=100,step=None,fix=None,metropolis=False,
                 markovpy=False):
        """
        NAME:
           sampleGP
        PURPOSE:
           sample the GP light curve model parameters
        INPUT:
           fix= None or list of parameters to hold fixed
           metropolis= if True, use Metropolis sampling
           markovpy= if True, use markovpy sampling
        OUTPUT:
           list of dictionaries of best-fit parameters
        HISTORY:
           2011-06-13 - Written - Bovy (NYU)
        """
        listofcfmfs= sampleGP(self.trainSet,self.cf,nsamples=nsamples,
                              step=step, mean=self.mf,fix=fix,
                              metropolis=metropolis,markovpy=markovpy)
        return [dict(o[0]._dict,**o[1]._dict) for o in listofcfmfs]

class TheoryLC():
    """Theoretical light curve, equally spaced abcissae(!)"""
    def __init__(self,trainSet=None,type='powerlawSF',params=None,band='r',
                 nx=201,mean='zero',extramean=0.):
        """
        NAME:
           __init__
        PURPOSE:
           initialize
        INPUT:
           trainSet - trainSet object of constraints
           type - type of lightcurve model
           params - parameters of the model
           band - band the model is in
           nx - number of gridpoints in time
           mean - mean to add
           extramean - another mean to add
        OUTPUT:
           (none)
        HISTORY:
           2010-12-24 - Written - Bovy (NYU)
        """
        self.trainSet= trainSet
        if len(band) == 1:
            if type == 'powerlawSF':
                try:
                    from flexgp.powerlawSF import covarFunc 
                except ImportError:
                    from gp.powerlawSF import covarFunc 
            elif type == 'brokenpowerlawSF':
                from brokenpowerlawSF import covarFunc 
            elif type == 'DRW':
                from DRW import covarFunc
            elif type == 'scatter':
                from scatterCovariance import covarFunc 
            elif type == 'zero':
                from zeroCovariance import covarFunc
            elif type == 'periodicDRW':
                try:
                    from flexgp.periodicOU_ARD import covarFunc 
                except ImportError:
                    from gp.periodicOU_ARD import covarFunc 
            if mean == 'zero':
                try:
                    from flexgp.zeroMean import meanFunc
                except ImportError:
                    from gp.zeroMean import meanFunc
            elif mean == 'const':
                try:
                    from flexgp.constMean import meanFunc
                except ImportError:
                    from gp.constMean import meanFunc
        else: #multi-filter
            if type == 'powerlawSF':
                from powerlawSFmulti import covarFunc
            elif type == 'powerlawSFratios':
                from powerlawSFmultiratios import covarFunc
            elif type == 'DRW':
                from OUgr import covarFunc
            elif type == 'zero':
                from zeroCovariance import covarFunc
            elif type == 'KS11':
                from KS11 import covarFunc
            elif type == 'KS11_nogamma':
                from KS11_nogamma import covarFunc
            elif type == 'KS11constscatter':
                from KS11constscatter import covarFunc
            elif type == 'KS11multiconstscatter':
                from KS11multiconstscatter import covarFunc
            elif type == 'KS11multi':
                from KS11multi import covarFunc
            elif type == 'reverbPLgr':
                from reverbPLgr import covarFunc
            if mean == 'zero':
                try:
                    from flexgp.zeroMean import meanFunc
                except ImportError:
                    from gp.zeroMean import meanFunc
            elif mean == 'const' \
                    and (type == 'KS11' or type == 'KS11constscatter' \
                             or type == 'KS11_nogamma' \
                             or type == 'reverbPLgr'):
                from KS11 import meanFunc
            elif mean == 'const' \
                    and (type == 'KS11multiconstscatter' \
                             or type == 'KS11multi'):
                from KS11multiconstscatter import meanFunc
        self.params= params
        self.cf= covarFunc(**self.params)
        self.mf= meanFunc(**self.params)
        self.band= band
        self.extramean= extramean
        #Now draw
        self.nx= nx
        if len(self.band) == 1:
            self.txs= nu.linspace(nu.amin(trainSet.listx),
                                  nu.amax(trainSet.listx)+1.,self.nx)
            self.xs= self.txs
        else:
            mjds= [x[0] for x in trainSet.listx]
            self.txs= nu.linspace(nu.amin(mjds),nu.amax(mjds),self.nx)
            xs= []
            for b in self.band:
                xs.extend([(self.txs[ii],b) for ii in range(self.nx)])
            self.xs= xs
        self.GPsample= eval_gp(self.xs,mean_func,covar_func,(self.mf),
                               (self.cf),nGP=1,
                               constraints=self.trainSet,
                               tiny_cholesky=0.000001).reshape(len(self.xs))
        if len(self.band) == 1:
            self.GPsample+= self.extramean
        else:
            for ii in range(len(self.band)):
                self.GPsample[ii*nx:(ii+1)*nx]+= self.extramean[ii]

    def plot(self,band,*args,**kwargs):
        """
        NAME:
           plot
        PURPOSE:
           plot a sample from the GP
        INPUT:
           band - band to plot (e.g., 'r' or 'gr')
           bovy_plot args and kwargs
        OUTPUT:
           (none)
        HISTORY:
           2010-12-21 - Written - Bovy (NYU)
        """
        kwargs['overplot']=True
        if len(self.band) == 1:
            if band != self.band:
                print "Requested filter is not the fitted filter ..."
                print "Returning ..."
                return None
            bovy_plot.bovy_plot(self.txs,self.GPsample,*args,**kwargs)
        elif len(self.band) > 1:
            if len(band) == 1:#plot magnitude
                one= nu.array([self.GPsample[ii] for ii in \
                                   range(len(self.GPsample)) if 
                               self.xs[ii][1] == band])
                bovy_plot.bovy_plot(self.txs,one,*args,**kwargs)
            elif len(band) == 2:#plot color
                one= nu.array([self.GPsample[ii] for ii in \
                                   range(len(self.GPsample)) if 
                               self.xs[ii][1] == band[0]])
                two= nu.array([self.GPsample[ii] for ii in \
                                  range(len(self.GPsample)) if 
                              self.xs[ii][1] == band[1]])
                bovy_plot.bovy_plot(self.txs,one-two,*args,**kwargs)
            else:
                print "Plotting more than two filters is not supported ..."
                print "Returning ..."
                return None

    def plotColorMag(self,band,color,*args,**kwargs):
        """
        NAME:
           plot
        PURPOSE:
           plot a sample from the GP in color vs. magnitude space
        INPUT:
           band - band to plot (e.g., 'r')
           color - color to plot (e.g., 'gr')
           bovy_plot args and kwargs
        OUTPUT:
           (none)
        HISTORY:
           2011-04-27 - Written - Bovy (NYU)
        """
        kwargs['overplot']=True
        if len(self.band) == 1:
            print "Cannot plot a color if the fit only contained one band ..."
            print "Returning ..."
            return None
        elif len(self.band) > 1:
            if len(band) > 1:
                print "'band' must be a single filter ..."
                print "Returning ..."
                return None
            else: #plot magnitude
                mag= nu.array([self.GPsample[ii] for ii in \
                                   range(len(self.GPsample)) if 
                               self.xs[ii][1] == band])
                if len(color) != 2:
                    print "'color' must be two filters ..."
                    print "Returning ..."
                    return None
                else:
                    one= nu.array([self.GPsample[ii] for ii in \
                                       range(len(self.GPsample)) if 
                                   self.xs[ii][1] == color[0]])
                    two= nu.array([self.GPsample[ii] for ii in \
                                       range(len(self.GPsample)) if 
                                   self.xs[ii][1] == color[1]])
                    bovy_plot.bovy_plot(mag,one-two,*args,**kwargs)

    def plotSF(self,band,*args,**kwargs):
        """
        NAME:
           plotSF
        PURPOSE:
           plot the structure function or color structure function
        INPUT:
           band - band to plot (e.g., 'r' for SF, 'gr' for CSF)
        OUTPUT:
           (none)
        HISTORY:
           2010-12-21 - Started - Bovy (NYU)
        """
        if len(self.band) == 1:
            if band != self.band:
                print "Requested filter is not the fitted filter ..."
                print "Returning ..."
                return None
            pyplot.loglog(sc.arange(1.,len(self.GPsample)/2)*
                          (self.txs[1]-self.txs[0]),
                          2.*sc.var(self.GPsample)
                          -2.*sc.correlate(self.GPsample
                                           -sc.mean(self.GPsample),
                                           self.GPsample
                                           -sc.mean(self.GPsample),"same")\
                              [1:len(self.GPsample)/2][::-1]/\
                              len(self.GPsample),
                          *args,**kwargs)
            xs= sc.arange(1.,len(self.GPsample)/2)*\
                (self.txs[1]-self.txs[0])
            ys= []
            for x in xs:
                ys.append(covar_func(0.,0.,(self.cf))\
                              -covar_func(0.,x,(self.cf)))
            ys= 2.*nu.array(ys).reshape(len(xs))
            pyplot.loglog(xs,ys,'k--')
        elif len(self.band) > 1:
            if len(band) == 1:#plot magnitude
                one= nu.array([self.GPsample[ii] for ii in \
                                   range(len(self.GPsample)) if 
                               self.xs[ii][1] == band])
                xs= sc.arange(1.,len(one)/2)*\
                    (self.txs[1]-self.txs[0])
                ys= []
                for x in xs:
                    ys.append(covar_func((0.,band),(0.,band),(self.cf))\
                                  -covar_func((0.,band),(x,band),(self.cf)))
                ys= 2.*nu.array(ys).reshape(len(xs))
                pyplot.loglog(xs,ys,'k--')
            elif len(band) == 2:#plot color
                one= nu.array([self.GPsample[ii] for ii in \
                                   range(len(self.GPsample)) if 
                               self.xs[ii][1] == band[0]])
                two= nu.array([self.GPsample[ii] for ii in \
                                  range(len(self.GPsample)) if 
                              self.xs[ii][1] == band[1]])
                one= one-two
                xs= sc.arange(1.,len(one)/2)*\
                    (self.txs[1]-self.txs[0])
                ys= []
                for x in xs:
                    ys.append(-2.*covar_func((0.,band[0]),(x,band[0]),
                                             (self.cf))\
                                   -2.*covar_func((0.,band[1]),(x,band[1]),
                                                  (self.cf))
                               +4.*covar_func((0.,band[0]),(x,band[1]),
                                              (self.cf))
                               +2.*covar_func((0.,band[0]),(0,band[0]),
                                              (self.cf))\
                                   +2.*covar_func((0.,band[1]),(0,band[1]),
                                                  (self.cf))
                               -4.*covar_func((0.,band[0]),(0,band[1]),
                                              (self.cf)))
                ys= nu.array(ys).reshape(len(xs))
                pyplot.loglog(xs,ys,'k--')
            else:
                print "Plotting more than two filters is not supported ..."
                print "Returning ..."
                return None
            pyplot.loglog(sc.arange(1.,len(one)/2)*
                          (self.txs[1]-self.txs[0]),
                          2.*sc.var(one)
                          -2.*sc.correlate(one
                                           -sc.mean(one),
                                           one
                                           -sc.mean(one),"same")\
                              [1:len(one)/2][::-1]/\
                              len(one),
                          *args,**kwargs)

def S2(xs,pm,pv,taus,_retoutnorm=False):
    """
    NAME:
       S_2
    PURPOSE:
       calculate S_2 = <(h(t+taus)-h(t))^2>
    INPUT:
       xs - points at which function is evaluated
       pm - function
       pv - error in function squared
       taus - offsets (integers) to consider as lags
    OUTPUT:
       S_2
    HISTORY:
       2012-10-11 - Written - Bovy (IAS)
    """
    out= sc.zeros(len(taus))
    norm= sc.zeros(len(taus))
    if isinstance(pm,list): #Multiple windows
        for ii in range(len(pm)):
            thisout, thisnorm= S2(xs[ii],pm[ii],pv[ii],taus,_retoutnorm=True)
            out+= thisout
            norm+= thisnorm
        return out/norm
    expanded_pm= sc.resize(pm,2*len(pm)+1)
    expanded_pm[len(pm):2*len(pm)+1]= 0.
    expanded_pv= sc.resize(pv,2*len(pm)+1)
    expanded_pv[len(pm):2*len(pm)+1]= 10.**6.
    for ii in range(len(taus)):
        out[ii]= sc.sum((sc.roll(expanded_pm,taus[ii])-expanded_pm)**2./(expanded_pv+sc.roll(expanded_pv,taus[ii])))
        norm[ii]= sc.sum(1./(expanded_pv+sc.roll(expanded_pv,taus[ii])))
    if _retoutnorm:
        return (out,norm)
    else:
        return out
 
def S3(xs,pm,pv,taus,_retoutnorm=False):
    """
    NAME:
       S3
    PURPOSE:
       calculate S3 = <(h(t)-h(t+tau))^3>
    INPUT:
       xs - points at which function is evaluated
       pm - function
       pv - error in function squared
       taus - offsets (integers) to consider as lags
    OUTPUT:
       S_3
    HISTORY:
       2012-10-11 - Written - Bovy (IAS)
    """
    out= sc.zeros(len(taus))
    norm= sc.zeros(len(taus))
    if isinstance(pm,list): #Multiple windows
        for ii in range(len(pm)):
            thisout, thisnorm= S3(xs[ii],pm[ii],pv[ii],taus,_retoutnorm=True)
            out+= thisout
            norm+= thisnorm
        return out/norm
    expanded_pm= sc.resize(pm,2*len(pm)+1)
    expanded_pm[len(pm):2*len(pm)+1]= 0.
    expanded_pv= sc.resize(pv,2*len(pm)+1)
    expanded_pv[len(pm):2*len(pm)+1]= 10.**6.
    for ii in range(len(taus)):
        out[ii]= sc.sum((sc.roll(expanded_pm,taus[ii])-expanded_pm)**3./(expanded_pv+sc.roll(expanded_pv,taus[ii])))
        norm[ii]= sc.sum(1./(expanded_pv+sc.roll(expanded_pv,taus[ii])))
    if _retoutnorm:
        return (out, norm)
    else:
        return out/norm

def skew(xs,pm,pv,taus):
    """
    NAME:
       skew
    PURPOSE:
       calculate skew = S_3 / s_2**1.5
    INPUT:
       xs - points at which function is evaluated
       pm - function
       pv - error in function squared
       taus - offsets (integers) to consider as lags
    OUTPUT:
       skew
    HISTORY:
       2012-10-11 - Written - Bovy (IAS)
    """
    return S3(xs,pm,pv,taus)/S2(xs,pm,pv,taus)**1.5

def panstarrs_sampling(nseasons=1,
                       startmjd=(2455317.500000-2400000.5)): #May 2010
    """
    NAME:
       panstarrs_sampling
    PURPOSE:
       return a panstarrs-sampling
    INPUT:
       nseassons - number of seasons
       startmjd - 
    OUTPUT:
       list of (t,band)
    HISTORY:
       2010-12-30 - Written - Bovy (NYU)
    """
    out= [[(-90./365.25+ii+startmjd/365.25,'z'),
           (-30./365.25+ii+startmjd/365.25,'Y'),
           (0./365.25+ii+startmjd/365.25,'i'),
           (5./365.25+ii+startmjd/365.25,'r'),
           (10./365.25+ii+startmjd/365.25,'g'),
           (30./365.25+ii+startmjd/365.25,'i'),
           (35./365.25+ii+startmjd/365.25,'r'),
           (40./365.25+ii+startmjd/365.25,'g'),
           (60./365.25+ii+startmjd/365.25,'Y'),
           (90./365.25+ii+startmjd/365.25,'z')] 
          for ii in range(nseasons)]
    return [item for sublist in out for item in sublist]

def sdss_sampling(startmjd=(2451635.5-2400000.5)): #April 2000
    """
    NAME:
       sdss_sampling
    PURPOSE:
       return one random SDSS epoch
    INPUT:
       startmjd - start MJD
    OUTPUT:
       list
    HISTORY:
       2011-01-31 - Bovy - NYU
    """
    rnd= nu.random.uniform()
    start= (2451635.5-2400000.5)/365.25 #April 2000
    end= (2453923.500000-2400000.5)/365.25 #July 2007
    depoch= 71.7/60./60./24./365.25 #dt between consecutive bands (from http://www.sdss3.org/instruments/camera.php)
    epoch= startmjd/365.25+rnd*(end-start)
    return [(epoch,'r'), #Robert
            (epoch+depoch,'i'), #is
            (epoch+2.*depoch,'u'), #under
            (epoch+3.*depoch,'z'), #ze
            (epoch+4.*depoch,'g')] #Gunn

def _load_fits(file,ext=1):
    """Loads fits file's data and returns it as a numpy.recarray with 
    lowercase field names"""
    hdulist= pyfits.open(file)
    out= hdulist[ext].data
    hdulist.close()
    return _as_recarray(out)

def _as_recarray(recarray):
    """go from FITS_rec to recarray"""
    newdtype = nu.dtype(recarray.dtype.descr)
    newdtype.names= tuple([n.lower() for n in newdtype.names])
    newrecarray = nu.recarray(recarray.shape, dtype=newdtype)
    for field in recarray.dtype.fields:
        newrecarray[field.lower()] = recarray.field(field)
    return newrecarray

def wedge_func(t,tstart,amp,tau=200./365.25):
    tend= tstart+tau
    if isinstance(t,nu.ndarray):
        out= nu.empty(t.shape)
        out= -amp/tau*(t-tend)
        out[(t < tstart)]= 0.
        out[(t > tend)]= 0.
        return out
    else:
        if t < tstart or t > tend:
            return 0.
        else:
            return -amp/tau*(t-tend)
