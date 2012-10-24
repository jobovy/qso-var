import scipy
from flexgp.covarianceClass import *
_MAXH= 10.
class covarFunc (covariance):
    """
    covarFunc twopowerlawSF: covariance function with a power-law
                          structure function

    powerlawSF(x,x')= s^2-A/2 (alpha |x-x'|^gamma_1 + (1-alpha) |x-x'|^gamma_2)

    """
    def __init__(self,**kwargs):
        """
        NAME:
           __init__
        PURPOSE:
           initialize a powerlawSF object
        OPTIONAL KEYWORD INPUTS:
           logA= or A=
           gamma=
        OUTPUT:
        HISTORY:
           2010-06-21 - Written - Bovy (NYU)
           2010-12-07 - sill = SF(\infty) - Bovy
        """
        self._dict= {}
        if kwargs.has_key('logA'):
            self.logA= kwargs['logA']
        elif kwargs.has_key('A'):
            self.logA= scipy.log(kwargs['A'])
        else:
            self.logA= 0.
        self._dict['logA']= self.logA
        if kwargs.has_key('gamma1'):
            self.gamma1= kwargs['gamma1']
        else:
            self.gamma1= 0.
        self._dict['gamma1']= self.gamma1
        if kwargs.has_key('gamma2'):
            self.gamma2= kwargs['gamma2']
        else:
            self.gamma2= 0.
        self._dict['gamma2']= self.gamma2
        if kwargs.has_key('alpha'):
            self.alpha= kwargs['alpha']
        else:
            self.alpha= 0.5
        self._dict['alpha']= self.alpha
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
           2010-06-21 - Written - Bovy (NYU)
        """
        if self.gamma1 > 2. or self.gamma1 < 0.: return -9999.99
        if self.gamma2 > 2. or self.gamma2 < 0.: return -9999.99
        if self.alpha > 1. or self.alpha < 0.: return -9999.99
        if not isinstance(x,numpy.ndarray):
            x= numpy.array(x)
        if not isinstance(xp,numpy.ndarray):
            xp= numpy.array(xp)
        return 0.5*(self._sf(_MAXH)-self._sf(x-xp))

    def _sf(self,x):
        if numpy.fabs(x) > _MAXH: return self._sf(_MAXH)
        return self.A*(self.alpha*numpy.fabs(x)**self.gamma1
                       +(1.-self.alpha)*numpy.fabs(x)**self.gamma2)
    
    def deriv(self,x,xp,key=None,covarValue=None):
        """
        NAME:
           deriv
        PURPOSE:
           derivative of the covariance function wrt a hyper-parameter or
           a pseudo-input if key == None
        INPUT:
           x - one point
           xp - another point
           key - key corresponding to the desired hyper-parameter in _dict or 'pseudo' if key == None
        OPTIONAL INPUT:
           covarValue - value of the covariance for these two points
        OUTPUT:
           derivative
        HISTORY:
           2010-07-05 - Written - Bovy (NYU)
        """
        raise NotImplementedError("'deriv' not implemented for twopowerlawSF")
        if not covarValue == None:
            covarValue= self.evaluate(x,xp)
        if key == 'logA':
            return -self.A/2.*numpy.fabs(x-xp)**self.gamma
        elif key == 'gamma':
            if x == xp:
                return 0.
            else:
                return -self.A/2.*numpy.fabs(x-xp)**self.gamma*numpy.log(numpy.fabs(x-xp))
        
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
        out['gamma1']= [True,True]
        out['gamma2']= [True,True]
        out['alpha']= [True,True]
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
        out['gamma1']= [0.,2.]
        out['gamma2']= [0.,2.]
        out['alpha']= [0.,1.]
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
        out['gamma1']= 'whole'
        out['gamma2']= 'whole'
        out['alpha']= 'whole'
        return out
