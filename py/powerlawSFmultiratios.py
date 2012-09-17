import scipy
from gp.covarianceClass import *
from gp.fast_cholesky import fast_cholesky_invert
_MAXH= 10.
_AGR= 0.21
_AIR= -0.07
_AZR= -0.26
_ACIZRI= 0.21
_ACGRRI= 0.14
_ACGIRI= 0.27
_ACGZRI= 0.49
_ACRZRI= 0.24
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
        if kwargs.has_key('logAr'):
            self.logAr= kwargs['logAr']
        elif kwargs.has_key('Ar'):
            self.logAr= scipy.log(kwargs['Ar'])
        else:
            self.logAr= 0.
        self._dict['logAr']= self.logAr
        if kwargs.has_key('logAri'):
            self.logAri= kwargs['logAri']
        elif kwargs.has_key('Ari'):
            self.logAri= scipy.log(kwargs['Ari'])
        else:
            self.logAri= 0.
        self._dict['logAri']= self.logAri
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
        self.Ar= scipy.exp(self.logAr)
        self.Ari= scipy.exp(self.logAri)

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
        #Get the relevant coefficiencts
        if x[1] == 'g':
            ASF1= self.logAr+2.*_AGR
        elif x[1] == 'r':
            ASF1= self.logAr
        elif x[1] == 'i':
            ASF1= self.logAr+2.*_AIR
        elif x[1] == 'z':
            ASF1= self.logAr+2.*_AZR
        if xp[1] == 'g':
            ASF2= self.logAr+2.*_AGR
        elif xp[1] == 'r':
            ASF2= self.logAr
        elif xp[1] == 'i':
            ASF2= self.logAr+2.*_AIR
        elif xp[1] == 'z':
            ASF2= self.logAr+2.*_AZR
        if not (x[1] == xp[1]):
            if (x[1] == 'g' and xp[1] == 'r') \
                    or (x[1] == 'r' and xp[1] == 'g'):
                AC= self.logAri+2.*_ACGRRI
            elif (x[1] == 'i' and xp[1] == 'z') \
                    or (x[1] == 'z' and xp[1] == 'i'):
                AC= self.logAri+2.*_ACIZRI
            elif (x[1] == 'i' and xp[1] == 'g') \
                    or (x[1] == 'g' and xp[1] == 'i'):
                AC= self.logAri+2.*_ACGIRI
            elif (x[1] == 'g' and xp[1] == 'z') \
                    or (x[1] == 'z' and xp[1] == 'g'):
                AC= self.logAri+2.*_ACGZRI
            elif (x[1] == 'r' and xp[1] == 'z') \
                    or (x[1] == 'z' and xp[1] == 'r'):
                AC= self.logAri+2.*_ACRZRI
            elif (x[1] == 'i' and xp[1] == 'r') \
                    or (x[1] == 'r' and xp[1] == 'i'):
                AC= self.logAri
            AC= numpy.exp(AC)
        else:
            AC= 0.
        ASF1= numpy.exp(ASF1)
        ASF2= numpy.exp(ASF2)
        return 0.25*(self._sf(_MAXH,ASF1)-self._sf(numpy.fabs(x[0]-xp[0]),ASF1)+
                     self._sf(_MAXH,ASF2)-self._sf(numpy.fabs(x[0]-xp[0]),ASF2))\
                     +1./4.*(self._csf(numpy.fabs(x[0]-xp[0]),AC)-self._csf(_MAXH,AC))
    
    def _sf(self,x,A):
        if x < 1./365./24.:
            return 0.
        return A*x**self.gamma

    def _csf(self,x,A):
        if x < 1./365./24.:
            return 0.
        return A*x**self.gammagr

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

