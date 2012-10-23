import os, os.path
import sys
import subprocess
import pypar
import numpy
from skewQSO import get_options
myid= pypar.rank()
node= pypar.get_processor_name()
#options
parser= get_options()
options,args= parser.parse_args()
#Parse cpu
if not options.split is None:
    ras= numpy.arange(options.split)
else:
    ras= [0,1,2,3,20,21,22,23,-2]
ra= ras[myid]
#cmd
cmd= [os.getenv('PYTHON'),
      'skewQSO.py',
      args[0]+'_%i' % ra,
      '-b',options.band,
      '--type='+options.type,
      '--mean='+options.mean,
      '--fitsfile='+options.fitsfile,
      '-n %i' % options.nsamples,
      '--saveevery=%i' % options.saveevery,
      '--rah=%i' % ra,
      '--dtau=%f'% options.dtau,
      '--taumax=%f' % options.taumax,
      '--wedgerate=%f' % options.wedgerate]
if options.wedge: cmd.append('--wedge')
if options.star: cmd.append('--star')
if options.rrlyrae: cmd.append('--rrlyrae')
if not options.split is None:
    cmd.append('--split=%i' % options.split)
#Now run
tryCall= True
while tryCall:
    try:
        subprocess.check_call(cmd)
        tryCall= False
    except (OSError,subprocess.CalledProcessError):
        print >> sys.stdout, node+": subprocess.check_call failed, re-trying ..."
        sys.stdout.flush()
print "Finished %f" % ra
