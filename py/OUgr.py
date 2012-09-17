import scipy
from gp.covarianceClass import *
from gp.fast_cholesky import fast_cholesky_invert
class covarFunc (covariance):
    """
    covarFunc OUgr: Ornstein-Uhlenbeck covariance function for both in-band 
    and cross-band correlations
    """
    def __init__(self,**kwargs):
        """
        NAME:
           __init__
        PURPOSE:
           initialize an OU_ARD object
        OPTIONAL KEYWORD INPUTS:
           loga2= or a=
           logl= or l= double, or list/array of doubles
           dim= dimension
        OUTPUT:
        HISTORY:
           2010-06-27 - Written - Bovy (NYU)
        """
        self._dict= {}
        if kwargs.has_key('loga2'):
            self.loga2= kwargs['loga2']
        elif kwargs.has_key('a'):
            self.loga2= 2.*scipy.log(kwargs['a'])
        else:
            self.loga2= 0.
        self._dict['loga2']= self.loga2
        if kwargs.has_key('logagr2'):
            self.logagr2= kwargs['logagr2']
        elif kwargs.has_key('agr'):
            self.logagr2= 2.*scipy.log(kwargs['agr'])
        else:
            self.logagr2= 0.
        self._dict['logagr2']= self.logagr2
        if kwargs.has_key('logl'):
            self.logl= kwargs['logl']
        elif kwargs.has_key('l'):
            self.logl= scipy.log(kwargs['l'])
        else:
            self.logl= numpy.array([0.])
        self._dict['logl']= self.logl
        if kwargs.has_key('loglgr'):
            self.loglgr= kwargs['loglgr']
        elif kwargs.has_key('lgr'):
            self.loglgr= scipy.log(kwargs['lgr'])
        else:
            self.loglgr= numpy.array([0.])
        self._dict['loglgr']= self.loglgr
        #Define shortcuts
        self.a2= scipy.exp(self.loga2)
        self.l= scipy.exp(self.logl)
        self.agr2= scipy.exp(self.logagr2)
        self.lgr= scipy.exp(self.loglgr)

    def evaluate(self,x,xp):
        """
        NAME:
           evaluate
        PURPOSE:
           evaluate the Ornstein-Uhlenbeck covariance function
        INPUT:
           x - one point
           xp - another point
        OUTPUT:
           covariance
        HISTORY:
           2010-06-27 - Written - Bovy (NYU)
        """
        if not (x[1] == xp[1]):
            return (self.a2-self._sf(numpy.fabs(x[0]-xp[0])))\
                +1./2.*(self._csf(numpy.fabs(x[0]-xp[0]))-self.agr2)
        else:
            return (self.a2-self._sf(numpy.fabs(x[0]-xp[0])))
    
    def _sf(self,x):
        return self.a2*(1.-numpy.exp(-x/self.l))

    def _csf(self,x):
        return self.agr2*(1.-numpy.exp(-x/self.lgr))

    def _list_params(self):
        """
        NAME:
           list_params
        PURPOSE:
           list all of the hyper-parameters of this covariance function
        INPUT:
        OUTPUT:
           (list of hyper-parameters (['a','l']),
           dimensionality of the hyper-parameter
        HISTORY:
           2010-02-15 - Written - Bovy (NYU)
        """
        return self._dict.keys()

