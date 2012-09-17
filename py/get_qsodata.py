import re
import sys
import os.path
import subprocess
import numpy as nu
import astrometry.util.pyfits_utils as pyfits_utils
import pyfits
from optparse import OptionParser
_PYTHON='/usr/local/epd/bin/python'
_CASJOBS='/global/data/scr/jb2777/astrometry/src/astrometry/util/casjobs.py'
def get_qsodata():
    """
    NAME:
       get_qsodata
    PURPOSE:
       get the QSO sample from the CAS
    INPUT:
       parser - from optParser
    OUTPUT:
    HISTORY:
       2010-12-21 - Written - Bovy (NYU)
    """
    caspwdfile= open('caspwd','r')
    casusr= caspwdfile.readline().rstrip()
    caspwd= caspwdfile.readline().rstrip()
    caspwdfile.close()
    #Read qso file
    s82file= '../data/S82qsos.fits'
    s82qsos= pyfits_utils.table_fields(s82file)[::-1]
    nqso= len(s82qsos.ra)
    done= []
    for qso in s82qsos:
        print "Working on QSO "+qso.oname.strip().replace(' ', '')
        if qso.oname.strip().replace(' ', '') in done:
            print qso.oname.strip().replace(' ', '')+" already done!"
            return
        done.append(qso.oname.strip().replace(' ', ''))
        tmpsavefilename= '../data/s82qsos/'+qso.oname.strip().replace(' ', '')+'.fit'
        if os.path.exists(tmpsavefilename):
            print "file "+tmpsavefilename+" exists"
            print "Delete file "+tmpsavefilename+" before running this to update the sample from the CAS"
        else:
            dbname= prepare_sql(qso)
            subprocess.call([_PYTHON,_CASJOBS,casusr,caspwd,'querywait',
                             '@tmp.sql'])
            subprocess.call([_PYTHON,_CASJOBS,casusr,caspwd,
                             'outputdownloaddelete',dbname,tmpsavefilename])

def prepare_sql(qso):
    output_f = open('tmp1.sql', 'w')
    subprocess.call(["sed",'s/QSORA/'+str(qso.ra).strip()+'/g',"s82qso.sql"],stdout=output_f)
    output_f.close()
    output_f = open('tmp2.sql', 'w')
    subprocess.call(["sed",'s/QSODEC/'+str(qso.dec).strip()+'/g','tmp1.sql'],stdout=output_f)
    output_f.close()
    output_f = open('tmp.sql', 'w')
    dbname= re.split(r'\.|\+|\-|\[|\]',str(qso.oname).strip().replace(' ',''))
    d= ''
    for db in dbname:
        d+= db
    dbname=d
    if d[0] in '0123456789':
        d= d[1:-1]
    subprocess.call(["sed",'s/DBNAME/'+d+'/g','tmp2.sql'],stdout=output_f)
    output_f.close()
    return d
       
if __name__ == '__main__':
    get_qsodata()
