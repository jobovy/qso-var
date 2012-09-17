import scipy
from gp.covarianceClass import *
from gp.meanClass import *
class covarFunc (covariance):
    """
    covarFunc KS11_nogamma: covariance function with a scatter
                  structure function for band g, and linear relation between 
                  g and r (eg)


    """
    def __init__(self,**kwargs):
        """
        NAME:
           __init__
        PURPOSE:
           initialize a KS11 object
        OPTIONAL KEYWORD INPUTS:
           logA= or A= double, or list/array of doubles
           s=
           B= 
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
        if kwargs.has_key('s'):
            self.s= kwargs['s']
        else:
            self.s= 0.
        self._dict['s']= self.s
        #Define shortcuts
        self.A= scipy.exp(self.logA)

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
        if not (x[1] == xp[1]):
            prefactor= self.s
        elif x[1] == xp[1] and x[1] != 'g':
            prefactor= self.s**2.
        else:
            prefactor= 1.
        x= numpy.fabs(x[0]-xp[0])
        if x < 1v.pl./365./24.:
            return self.A/2.*prefactor
        else:
            return 0.

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

