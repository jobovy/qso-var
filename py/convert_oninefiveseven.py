import sys
import os, os.path
import csv
import numpy
import fitsio
def convert_oninefiveseven(type='vanderriest'):
    #Read the file
    if type == 'vanderriest':
        file= '../data/0957.dat'
    elif type == 'kundic_g':
        file= '../data/L0957_g.DAT'
    elif type == 'kundic_r':
        file= '../data/L0957_r.DAT'
    dialect= csv.excel
    dialect.skipinitialspace=True
    reader= csv.reader(open(file,'r'),delimiter=' ',
                       dialect=dialect)
    jd, amag, amagerr, bmag, bmagerr= [], [], [], [], []
    for row in reader:
        if row[0][0] == '#':
            continue
        jd.append(float(row[0]))
        amag.append(float(row[1]))
        amagerr.append(float(row[2]))
        bmag.append(float(row[3]))
        bmagerr.append(float(row[4]))
    #A
    ndata= len(jd)
    if type == 'kundic_g':
        out= numpy.recarray((ndata,),
                            dtype=[('mjd_g', 'f8'),
                                   ('g', 'f8'),
                                   ('err_g', 'f8')])
        out.mjd_g= 2440000-2400000.5+numpy.array(jd) #not sure why I'm bothering
        out.g= amag
        out.err_g= amagerr
    else:
        out= numpy.recarray((ndata,),
                            dtype=[('mjd_r', 'f8'),
                                   ('r', 'f8'),
                                   ('err_r', 'f8')])
        out.mjd_r= 2440000-2400000.5+numpy.array(jd) #not sure why I'm bothering
        out.r= amag
        out.err_r= amagerr
    if type == 'vanderriest':
        fitsio.write('../data/0957-A.fits',out,clobber=True)
    elif type == 'kundic_g':
        fitsio.write('../data/L0957-A_g.fits',out,clobber=True)
    elif type == 'kundic_r':
        fitsio.write('../data/L0957-A_r.fits',out,clobber=True)
    #B
    if type == 'kundic_g':
        out= numpy.recarray((ndata,),
                            dtype=[('mjd_g', 'f8'),
                                   ('g', 'f8'),
                                   ('err_g', 'f8')])
        out.mjd_g= 2440000-2400000.5+numpy.array(jd) #not sure why I'm bothering
        out.g= bmag
        out.err_g= bmagerr
    else:
        out= numpy.recarray((ndata,),
                            dtype=[('mjd_r', 'f8'),
                                   ('r', 'f8'),
                                   ('err_r', 'f8')])
        out.mjd_r= 2440000-2400000.5+numpy.array(jd) #not sure why I'm bothering
        out.r= bmag
        out.err_r= bmagerr
    if type == 'vanderriest':
        fitsio.write('../data/0957-B.fits',out,clobber=True)
    elif type == 'kundic_g':
        fitsio.write('../data/L0957-B_g.fits',out,clobber=True)
    elif type == 'kundic_r':
        fitsio.write('../data/L0957-B_r.fits',out,clobber=True)
    return None

if __name__ == '__main__':
    if len(sys.argv) > 1:
        convert_oninefiveseven(type=sys.argv[1])
    else:
        convert_oninefiveseven()
