import subprocess
import pypar
from skewQSO import get_options
myid= pypar.rank()
node= pypar.get_processor_name()
ras= [0,1,2,3,20,21,22,23,-2]
ra= ras[myid]
#cmd
cmd= [os.getenv('PYTHON'),
      'skewQSO.py',
      args[0]+'_%i' % ra,
      '-b '+options.band,
      '--type='+options.type,
      '--mean='+options.mean,
      '--fitsfile='+options.fitsfile,
      '-n %i' % options.nsamples,
      '--saveevery=%i' % options.saveevery,
      '--rah=%i' % ra,
      '--dtau=%f'+options.dtau,
      '--taumax=%f' % options.taumax,
      '--wedgerate=%f' % options.wedgerate]
if options.wedge: cmd.append('--wedge')
#Now run
tryCall= True
while tryCall:
    try:
        subprocess.check_call(cmd)
        tryCall= False
    except (OSError,subprocess.CalledProcessError):
        if options.mpi:
            print >> sys.stdout, node+": subprocess.check_call failed, re-trying ..."
            sys.stdout.flush()
        else:
            print >> sys.stdout, "subprocess.check_call failed, re-trying ..."
    print "Finished %f" % ra
