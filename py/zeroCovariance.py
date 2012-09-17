import scipy
from gp.covarianceClass import *
from gp.meanClass import *
from gp.fast_cholesky import fast_cholesky_invert
_MAXH= 10.
class covarFunc (covariance):
    """
    covarFunc zeroCovariance: zero covariance function
    """
    def __init__(self,**kwargs):
        """
        NAME:
           __init__
        PURPOSE:
           initialize
        OPTIONAL KEYWORD INPUTS:
        OUTPUT:
        HISTORY:
           2011-06-10 - Written - Bovy (NYU)
        """
        self._dict= {}

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

