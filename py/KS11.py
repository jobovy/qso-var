import scipy
from gp.covarianceClass import *
from gp.meanClass import *
from gp.fast_cholesky import fast_cholesky_invert
_MAXH= 10.
class meanFunc (mean):
    """
    meanFunc KS11: mean function of the KS11 linear relation between g and r 
    model
    """
    def __init__(self,**kwargs):
        """
        NAME:
           __init__
        PURPOSE:
           initialize a KS11 mean object
        OPTIONAL KEYWORD INPUTS:
           m=
           s=
           B= 
        OUTPUT:
        HISTORY:
           2011-05-19 - Written - Bovy (NYU)
        """
        self._dict= {}
        if kwargs.has_key('b'):
            self.b= kwargs['b']
        else:
            self.b= 0.
        self._dict['b']= self.b
        if kwargs.has_key('s'):
            self.s= kwargs['s']
        else:
            self.s= 0.
        self._dict['s']= self.s
        if kwargs.has_key('m'):
            self.m= kwargs['m']
        else:
            self.m= 0.
        self._dict['m']= self.m

    def evaluate(self,x):
        """
        NAME:
           evaluate
        PURPOSE:
           evaluate the mean
        INPUT:
           x - one point
        OUTPUT:
           mean
        HISTORY:
           2011-05-19 - Written - Bovy (NYU)
        """
        if x[1] != 'g':
            return self.s*self.m+self.b
        else:
            return self.m
    
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

        
class covarFunc (covariance):
    """
    covarFunc KS11: covariance function with a power-law
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
           gamma=
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
        if not (x[1] == xp[1]):
            prefactor= self.s
        elif x[1] == xp[1] and x[1] != 'g':
            prefactor= self.s**2.
        else:
            prefactor= 1.
        return prefactor*0.5*(self._sf(_MAXH)-self._sf(numpy.fabs(x[0]-xp[0])))
    
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
        out['gamma']= [True,True] #All but gamma are infinite
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
        out['gamma']= [0.,2.] #All but gamma are infinite
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
        out['gamma']= 'whole' #All but gamma are default, stepping out=quick
        return out
