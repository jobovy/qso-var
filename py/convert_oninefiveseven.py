import sys
import os, os.path
import csv
import numpy
import fitsio
def convert_oninefiveseven():
    #Read the file
    file= '../data/0957.dat'
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
    out= numpy.recarray((ndata,),
                        dtype=[('mjd_r', 'f8'),
                               ('r', 'f8'),
                               ('err_r', 'f8')])
    out.mjd_r= 2440000-2400000.5+numpy.array(jd) #not sure why I'm bothering
    out.r= amag
    out.err_r= amagerr
    fitsio.write('../data/0957-A.fits',out,clobber=True)
    #B
    out= numpy.recarray((ndata,),
                        dtype=[('mjd_r', 'f8'),
                               ('r', 'f8'),
                               ('err_r', 'f8')])
    out.mjd_r= 2440000-2400000.5+numpy.array(jd) #not sure why I'm bothering
    out.r= bmag
    out.err_r= bmagerr
    fitsio.write('../data/0957-B.fits',out,clobber=True)
    return None

if __name__ == '__main__':
    convert_oninefiveseven()
