import scipy
from gp.covarianceClass import *
from gp.fast_cholesky import fast_cholesky_invert
_MAXH= 10.
class covarFunc (covariance):
    """
    covarFunc powerlawSFmulti: covariance function with a power-law
                          structure function for bands g and r (eg)

    powerlawSF(x,x')= A-2 (Gamma|x-x'|)^gamma/2 - Bgr

    """
    def __init__(self,**kwargs):
        """
        NAME:
           __init__
        PURPOSE:
           initialize a powerlawSFgr object
        OPTIONAL KEYWORD INPUTS:
           logGamma= or Gamma= double, or list/array of doubles
           gamma=
           gammagr=
        OUTPUT:
        HISTORY:
           2010-08-11 - Written - Bovy (NYU)
        """
        self._dict= {}
        if kwargs.has_key('logA'):
            self.logA= kwargs['logA']
        elif kwargs.has_key('A'):
            self.logA= scipy.log(kwargs['A'])
        else:
            self.logA= 0.
        self._dict['logA']= self.logA
        if kwargs.has_key('logAgr'):
            self.logAgr= kwargs['logAgr']
        elif kwargs.has_key('Agr'):
            self.logAgr= scipy.log(kwargs['Agr'])
        else:
            self.logAgr= 0.
        self._dict['logAgr']= self.logAgr
        if kwargs.has_key('gamma'):
            self.gamma= kwargs['gamma']
        else:
            self.gamma= 0.
        self._dict['gamma']= self.gamma
        if kwargs.has_key('gammagr'):
            self.gammagr= kwargs['gammagr']
        else:
            self.gammagr= 0.
        self._dict['gammagr']= self.gammagr
        #Define shortcuts
        self.A= scipy.exp(self.logA)
        self.Agr= scipy.exp(self.logAgr)

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
        if self.gammagr > 2. or self.gammagr < 0.: return -9999.99
        if not (x[1] == xp[1]):
            return 0.5*(self._sf(_MAXH)-self._sf(numpy.fabs(x[0]-xp[0])))\
                +1./4.*(self._csf(numpy.fabs(x[0]-xp[0]))-self._csf(_MAXH))
        else:
            return 0.5*(self._sf(_MAXH)-self._sf(numpy.fabs(x[0]-xp[0])))
    
    def _sf(self,x):
        return self.A*x**self.gamma

    def _csf(self,x):
        return self.Agr*x**self.gammagr

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

    def isDomainFinite(self):
        """
        NAME:
           isDomainFinite
        PURPOSE:
           return dictionary that says whether the hyperparameters' domains are finite
        INPUT:
        OUTPUT:
           boolean list
        HISTORY:
           2011-06-13 - Written - Bovy (NYU)
        """
        out= covariance.isDomainFinite(self)
        out['gamma']= [True,True]
        out['gammagr']= [True,True]
        return out

    def paramsDomain(self):
        """
        NAME:
           paramsDomain
        PURPOSE:
           return dictionary that has each hyperparameter's domain 
           (irrelevant for hyperparameters with infinite domains)
        INPUT:
        OUTPUT:
           dictionary of lists
        HISTORY:
           2011-06-13 - Written - Bovy (NYU)
        """
        out= covariance.paramsDomain(self)
        out['gamma']= [0.,2.]
        out['gammagr']= [0.,2.]
        return out

    def create_method(self):
        """
        NAME:
           create_method
        PURPOSE:
           return dictionary that has each hyperparameter's create_method
           for slice sampling
        INPUT:
        OUTPUT:
           dictionary of methods
        HISTORY:
           2011-06-13 - Written - Bovy (NYU)
        """
        out= covariance.create_method(self)
        out['gamma']= 'whole'
        out['gammagr']= 'whole'
        return out
