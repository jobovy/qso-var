import math
import scipy
from gp.covarianceClass import *
from gp.meanClass import *
from gp.fast_cholesky import fast_cholesky_invert
_DEBUG= True
_MAXH= 10.
class covarFunc (covariance):
    """
    covarFunc reverbPLgr: covariance function with a power-law
                  structure function for band g, linear relation between 
                  g and r (eg), and a lag in g


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
           beta=
           eps=
           logtL= or tL=
           logsigmat= or sigmat=
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
        if kwargs.has_key('logtL'):
            self.logtL= kwargs['logtL']
        elif kwargs.has_key('tL'):
            self.logtL= scipy.log(kwargs['tL'])
        else:
            self.logtL= 0.
        self._dict['logtL']= self.logtL
        if kwargs.has_key('logsigmatL'):
            self.logsigmatL= kwargs['logsigmatL']
        elif kwargs.has_key('sigmatL'):
            self.logsigmatL= scipy.log(kwargs['sigmatL'])
        else:
            self.logsigmatL= 0.
        self._dict['logsigmatL']= self.logsigmatL
        if kwargs.has_key('beta'):
            self.beta= kwargs['beta']
        else:
            self.beta= 0.1
        self._dict['beta']= self.beta
        if kwargs.has_key('eps'):
            self.eps= kwargs['eps']
        else:
            self.eps= 0.1
        self._dict['eps']= self.eps
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
        self.tL= scipy.exp(self.logtL)
        self.sigmatL= scipy.exp(self.logsigmatL)
        if _DEBUG: print self._dict

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
           2011-06-20 - Written - Bovy (NYU)
        """
        if self.gamma > 2. or self.gamma < 0.: return -9999.99
        if self.beta > 2. or self.beta < 0.: return -9999.99
        if self.eps < 0.: return -9999.99
        if self.logsigmatL < math.log(1./365./24.) \
                or self.logsigmatL > math.log(10.): return -9999.99
        if self.logtL < math.log(1./365./24.) \
                or self.logtL > math.log(10.): return -9999.99
        if self.logsigmatL > self.logtL: return -9999.99
        if not (x[1] == xp[1]) and x[1] == 'g':
            return self._evaluate_gr(x[0]-xp[0])
        if not (x[1] == xp[1]) and xp[1] == 'g': #not used since symmetric
            return self._evaluate_rg(x[0]-xp[0])
        elif x[1] == xp[1] and x[1] != 'g':
            return self._evaluate_rr(math.fabs(x[0]-xp[0]))
        else:
            return self._evaluate_gg(x[0]-xp[0])

    def _evaluate_rr(self,dt):
        return 0.5*self.s**2.*(self._sf(_MAXH)-self._sf(dt))

    def _evaluate_gr(self,dt):
        return ((1.-self.beta)*0.5*self.s*(self._sf(_MAXH)
                                           -self._sf(math.fabs(dt)))
                +self.beta*0.5*self.s*self.eps*self._one_transfer_int(dt,plus=True))

    def _evaluate_rg(self,dt):
        return ((1.-self.beta)*0.5*self.s*(self._sf(_MAXH)
                                           -self._sf(math.fabs(dt)))
                +self.beta*0.5*self.s*self.eps*self._one_transfer_int(dt,plus=False))

    def _evaluate_gg(self,dt):
        return ((1.-self.beta)**2.*0.5*(self._sf(_MAXH)-self._sf(math.fabs(dt)))
                +0.5*self.beta**2.*self.eps**2.*self._two_transfer_int(dt)
                +(1.-self.beta)*self.beta*0.5*self.eps*(self._one_transfer_int(dt,plus=True)+self._one_transfer_int(dt,plus=False)))


    def _two_transfer_int(self,dt):
        return self._sf(_MAXH)-1./8./self.sigmatL**2.*self.A/(self.gamma+1.)\
            /(self.gamma+2.)*(math.fabs(dt+2.*self.sigmatL)**(self.gamma+2.)
                              -2.*math.fabs(dt)**(self.gamma+2.)
                              +math.fabs(dt-2.*self.sigmatL)**(self.gamma+2.))

    def _one_transfer_int(self,dt,plus=True):
        if plus:
            dtplus= dt+self.tL+self.sigmatL
            dtminus= dt+self.tL-self.sigmatL
        else:
            dtplus= dt-self.tL+self.sigmatL
            dtminus= dt-self.tL-self.sigmatL
        if dtminus <= 0. and dtplus >= 0.:
            out= math.fabs(dtplus)**(self.gamma+1.)\
                +math.fabs(dtminus)**(self.gamma+1.)
        elif dtminus > 0.:
            out= math.fabs(dtplus)**(self.gamma+1.)\
                -math.fabs(dtminus)**(self.gamma+1.)
        else:
                out= math.fabs(dtminus)**(self.gamma+1.)\
                    -math.fabs(dtplus)**(self.gamma+1.)
        return self._sf(_MAXH)-1./4./self.sigmatL*self.A/(self.gamma+1.)*out
            
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
        out['gamma']= [True,True]
        out['beta']= [True,True]
        out['eps']= [True,False]
        out['logsigmatL']= [True,True]
        out['logtL']= [True,True]
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
        out['beta']= [0.,1.]
        out['eps']= [0.,0.]
        out['logsigmatL']= [math.log(1./365./24.),math.log(10.)]
        out['logtL']= [math.log(1./365./24.),math.log(10.)]
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
