import scipy
from gp.covarianceClass import *
from gp.meanClass import *
from gp.fast_cholesky import fast_cholesky_invert
_MAXH= 10.
class covarFunc (covariance):
    """
    covarFunc KS11constscatter: covariance function with a power-law
                  structure function for band g, and linear relation between 
                  g and r (eg), and constant scatter in this relation


    """
    def __init__(self,**kwargs):
        """
        NAME:
           __init__
        PURPOSE:
           initialize a KS11 object
        OPTIONAL KEYWORD INPUTS:
           logA= or A= double, or list/array of doubles
           gamma=
           s=
           B= 
           sigma= or logsigma= 
        OUTPUT:
        HISTORY:
           2011-04-27 - Written - Bovy (NYU)
        """
        self._dict= {}
        if kwargs.has_key('logA'):
            self.logA= kwargs['logA']
        elif kwargs.has_key('A'):
            self.logA= scipy.log(kwargs['A'])
        else:
            self.logA= 0.
        self._dict['logA']= self.logA
        if kwargs.has_key('logsigma'):
            self.logsigma= kwargs['logsigma']
        elif kwargs.has_key('sigma'):
            self.logsigma= scipy.log(kwargs['sigma'])
        else:
            self.logsigma= 0.
        self._dict['logsigma']= self.logsigma
        if kwargs.has_key('gamma'):
            self.gamma= kwargs['gamma']
        else:
            self.gamma= 0.
        self._dict['gamma']= self.gamma
        if kwargs.has_key('s'):
            self.s= kwargs['s']
        else:
            self.s= 0.
        self._dict['s']= self.s
        #Define shortcuts
        self.A= scipy.exp(self.logA)
        self.sigma= scipy.exp(self.logsigma)

    def evaluate(self,x,xp):
        """
        NAME:
           evaluate
        PURPOSE:
           evaluate the power-law SF covariance function
        INPUT:
           x - one point
           xp - another point
        OUTPUT:
           covariance
        HISTORY:
           2010-08-11 - Written - Bovy (NYU)
        """
        if self.gamma > 2. or self.gamma < 0.: return -9999.99
        extraterm= 0.
        if not (x[1] == xp[1]):
            prefactor= self.s
        elif x[1] == xp[1] and x[1] != 'g':
            prefactor= self.s**2.
            if numpy.fabs(x[0]-xp[0]) < 1./365./24.:
                extraterm= self.sigma
        else:
            prefactor= 1.
        return prefactor*0.5*(self._sf(_MAXH)-self._sf(numpy.fabs(x[0]-xp[0])))+extraterm
    
    def _sf(self,x):
        return self.A*x**self.gamma

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

