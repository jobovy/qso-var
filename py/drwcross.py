import scipy
from gp.covarianceClass import *
from gp.fast_cholesky import fast_cholesky_invert
class covarFunc (covariance):
    """
    covarFunc drw-cross: covariance function with a DRW in-band
                          structure function for bands g and r (eg)

    DRW = S exp(-|t-t'|/tau

    """
    def __init__(self,**kwargs):
        """
        NAME:
           __init__
        PURPOSE:
           initialize a drw-cross object
        OPTIONAL KEYWORD INPUTS:
           logS= or S=
           logtau= or tau=
           B= or logB=
        OUTPUT:
        HISTORY:
           2010-08-11 - Written - Bovy (NYU)
        """
        self._dict= {}
        if kwargs.has_key('logS'):
            self.logS= kwargs['logS']
        elif kwargs.has_key('S'):
            self.logS= scipy.log(kwargs['S'])
        else:
            self.logS= 0.
        self._dict['logS']= self.logS
        if kwargs.has_key('logtau'):
            self.logtau= kwargs['logtau']
        elif kwargs.has_key('tau'):
            self.logtau= scipy.log(kwargs['tau'])
        else:
            self.logtau= 0.
        self._dict['logtau']= self.logtau
        if kwargs.has_key('logB'):
            self.logB= kwargs['logB']
        elif kwargs.has_key('B'):
            self.logB= scipy.log(kwargs['B'])
        else:
            self.logB= 0.
        self._dict['logB']= self.logB
        if kwargs.has_key('logC'):
            self.logC= kwargs['logC']
        elif kwargs.has_key('C'):
            self.logC= scipy.log(kwargs['C'])
        else:
            self.logC= 0.
        self._dict['logC']= self.logC
        if kwargs.has_key('gammagr'):
            self.gammagr= kwargs['gammagr']
        else:
            self.gammagr= 0.
        self._dict['gammagr']= self.gammagr
        if kwargs.has_key('logGammagr'):
            self.logGammagr= kwargs['logGammagr']
        elif kwargs.has_key('Gammagr'):
            self.logGammagr= scipy.log(kwargs['Gammagr'])
        else:
            self.logGammagr= 0.
        self._dict['logGammagr']= self.logGammagr
        #Define shortcuts
        self.S= scipy.exp(self.logS)
        self.B= scipy.exp(self.logB)
        self.C= scipy.exp(self.logC)
        self.Gammagr= scipy.exp(self.logGammagr)
        self.tau= scipy.exp(self.logtau)

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
        if (x[1] == 'g' and xp[1] == 'r') or (x[1] == 'r' and xp[1] == 'g'):
            thisx= x
            thisxp= (xp[0],x[1])
            return self.evaluate(thisx,thisxp)\
                -self.B+self.C*numpy.exp(-(self.Gammagr*numpy.fabs(x[0]-xp[0]))**(self.gammagr))
#                -self.B+(self.Gammagr*numpy.fabs(x[0]-xp[0]))**(self.gammagr)
#                -self.B
        else:
            return self.S*numpy.exp(-numpy.fabs(x[0]-xp[0])/self.tau)
    
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

