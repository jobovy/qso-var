import scipy
from gp.covarianceClass import *
from gp.fast_cholesky import fast_cholesky_invert
_MAXH= 10.
class covarFunc (covariance):
    """
    covarFunc powerlawSFgr: covariance function with a power-law
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
        if kwargs.has_key('logGamma'):
            self.logGamma= kwargs['logGamma']
        elif kwargs.has_key('Gamma'):
            self.logGamma= scipy.log(kwargs['Gamma'])
        else:
            self.logGamma= 0.
        self._dict['logGamma']= self.logGamma
        if kwargs.has_key('logGammagr'):
            self.logGammagr= kwargs['logGammagr']
        elif kwargs.has_key('Gammagr'):
            self.logGammagr= scipy.log(kwargs['Gammagr'])
        else:
            self.logGammagr= 0.
        self._dict['logGammagr']= self.logGammagr
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
        self.Gamma= scipy.exp(self.logGamma)
        self.Gammagr= scipy.exp(self.logGammagr)

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
        #if (x[1] == 'g' and xp[1] == 'r') or (x[1] == 'r' and xp[1] == 'g'):
        if not (x[1] == xp[1]):
            return 0.5*(self._sf(_MAXH)-self._sf(numpy.fabs(x[0]-xp[0])))\
                +1./4.*(self._csf(numpy.fabs(x[0]-xp[0]))-self._csf(_MAXH))
        else:
            return 0.5*(self._sf(_MAXH)-self._sf(numpy.fabs(x[0]-xp[0])))
    
    def _sf(self,x):
        if x > _MAXH: return self._sf(_MAXH)
        return (self.Gamma*x)**self.gamma

    def _csf(self,x):
        if x > _MAXH: return self._csf(_MAXH)
        return (self.Gammagr*x)**self.gammagr

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
        raise AttributeError #Not implemented
        if not covarValue == None:
            covarValue= self.evaluate(x,xp)
        if key == 'logA':
            return -self.A/2.*numpy.fabs(x-xp)**self.gamma
        elif key == 'logs2':
            return self.s2
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

