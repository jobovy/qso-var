import scipy
from gp.covarianceClass import *
from gp.meanClass import *
from gp.fast_cholesky import fast_cholesky_invert
_MAXH= 10.
class meanFunc (mean):
    """
    meanFunc KS11multi:
    mean function of the KS11 linear relation between g and r,i,z
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
           sr, si, sz=
           br, bi, bz= 
        OUTPUT:
        HISTORY:
           2011-05-19 - Written - Bovy (NYU)
        """
        self._dict= {}
        if kwargs.has_key('br'):
            self.br= kwargs['br']
        else:
            self.br= 0.
        self._dict['br']= self.br
        if kwargs.has_key('bi'):
            self.bi= kwargs['bi']
        else:
            self.bi= 0.
        self._dict['bi']= self.bi
        if kwargs.has_key('bz'):
            self.bz= kwargs['bz']
        else:
            self.bz= 0.
        self._dict['bz']= self.bz
        if kwargs.has_key('sr'):
            self.sr= kwargs['sr']
        else:
            self.sr= 0.
        self._dict['sr']= self.sr
        if kwargs.has_key('si'):
            self.si= kwargs['si']
        else:
            self.si= 0.
        self._dict['si']= self.si
        if kwargs.has_key('sz'):
            self.sz= kwargs['sz']
        else:
            self.sz= 0.
        self._dict['sz']= self.sz
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
        if x[1] == 'r':
            return self.sr*self.m+self.br
        elif x[1] == 'i':
            return self.si*self.m+self.bi
        if x[1] == 'z':
            return self.sz*self.m+self.bz
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
    covarFunc KS11multi: covariance function with a power-law
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
           sr, si, sz=
           Br, Bi, Bz= 
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
        if kwargs.has_key('sr'):
            self.sr= kwargs['sr']
        else:
            self.sr= 0.
        self._dict['sr']= self.sr
        if kwargs.has_key('si'):
            self.si= kwargs['si']
        else:
            self.si= 0.
        self._dict['si']= self.si
        if kwargs.has_key('sz'):
            self.sz= kwargs['sz']
        else:
            self.sz= 0.
        self._dict['sz']= self.sz
        #Define shortcuts
        self.A= scipy.exp(self.logA)
        self.ss= {'g':{'g':1.,'r':self.sr,'i':self.si,'z':self.sz},
                  'r':{'g':self.sr,'r':self.sr**2.,'i':self.si*self.sr,
                       'z':self.sz*self.sr},
                  'i':{'g':self.si,'r':self.si*self.sr,'i':self.si**2.,
                       'z':self.sz*self.si},
                  'z':{'g':self.sz,'r':self.sz*self.sr,'i':self.si*self.sz,
                       'z':self.sz**2.}}

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
        prefactor= self.ss[x[1]][xp[1]]
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

