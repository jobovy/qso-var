import os, os.path
import sys
import cPickle as pickle
import numpy
from galpy.util import save_pickles
from skewQSO import get_options, _ERASESTR
#options
parser= get_options()
options,args= parser.parse_args()
#Parse
if not options.split is None:
    ras= numpy.arange(options.split)+1
else:
    ras= [0,1,2,3,20,21,22,23,-2]
#cmd
savefilebase= args[0]
#savefiles= [args[0]+'_%i' % ra for ra in ras]
spl= args[0].split('.')
newname= ''
for jj in range(len(spl)-1):
    newname+= spl[jj]
    if not jj == len(spl)-2: newname+= '.'
savefiles= [newname+'_%i.' % ra + spl[-1] for ra in ras]
ii= 0
if os.path.exists(savefiles[ii]):
    savefile= open(savefiles[ii],'rb')
    skews= pickle.load(savefile)
    gaussskews= pickle.load(savefile)
    type= pickle.load(savefile)
    band= pickle.load(savefile)
    mean= pickle.load(savefile)
    taus= pickle.load(savefile)      
    savefile.close()
else:
    raise IOError("Savefile %s does not exist ..." % savefiles[ii])
outskews= [skews]
outgaussskews= [gaussskews]
for ii in range(1,len(savefiles)):
    if os.path.exists(savefiles[ii]):
        savefile= open(savefiles[ii],'rb')
        skews= pickle.load(savefile)
        gaussskews= pickle.load(savefile)
        savefile.close()
    else:
        raise IOError("Savefile %s does not exist ..." % savefiles[ii])
    outskews.append(skews)
    outgaussskews.append(gaussskews)
#Now load into out dict
out= dict(outskews[0])
outgauss= dict(outgaussskews[0])
for ii in range(1,len(savefiles)):
    out.update(outskews[ii])
    outgauss.update(outgaussskews[ii])
#Save
save_pickles(args[0],out,outgauss,type,band,mean,taus)
