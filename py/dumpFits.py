import os, os.path
import numpy
import cPickle as pickle
from optparse import OptionParser
import pyfits
from astrometry.util import pyfits_utils as fu
from plotFits import open_qsos
def dumpFits(parser):
    (options,args)= parser.parse_args()
    if len(args) == 0:
        parser.print_help()
        return
    if options.outfilename is None:
        print "-o filename options needs to be set ..."
        print "Returning ..."
        return None
    #load fits
    if os.path.exists(args[0]):
        savefile= open(args[0],'rb')
        params= pickle.load(savefile)
        type= pickle.load(savefile)
        band= pickle.load(savefile)
        mean= pickle.load(savefile)
        savefile.close()
    else:
        raise IOError("input file does not exist ...")
    #dump fits, first the arrays
    cols= []
    #Get basic quasar stuff
    qsos= open_qsos(options.inputfilename)
    qsoDict= {}
    ii=0
    for qso in qsos:
        qsoDict[qso.oname.strip().replace(' ', '')+'.fit']= ii
        ii+= 1
    paramsKeys= params[params.keys()[0]].keys()
    qsoKeys= params.keys()
    for key in paramsKeys:
        thiscol= []
        for qso in qsoKeys:
            try:
                thiscol.append(params[qso][key][0])
            except IndexError:
                thiscol.append(params[qso][key])
        cols.append(pyfits.Column(name=key,format='D',
                                  array=numpy.array(thiscol)))
    if not options.star and not options.nuvx and not options.uvx \
            and not options.rrlyrae and not options.nuvxall:
        #RA, Dec, name, z, key
        thiscol= []
        for qso in qsoKeys:
            thiscol.append(qsos[qsoDict[qso]].ra)
        cols.append(pyfits.Column(name='ra',format='D',
                                  array=numpy.array(thiscol)))
        thiscol= []
        for qso in qsoKeys:
            thiscol.append(qsos[qsoDict[qso]].dec)
        cols.append(pyfits.Column(name='dec',format='D',
                                  array=numpy.array(thiscol)))
        thiscol= []
        for qso in qsoKeys:
            thiscol.append(qsos[qsoDict[qso]].oname)
        cols.append(pyfits.Column(name='oname',format='30A',
                                  array=numpy.array(thiscol)))
        thiscol= []
        for qso in qsoKeys:
            thiscol.append(qsos[qsoDict[qso]].z)
        cols.append(pyfits.Column(name='z',format='D',
                                  array=numpy.array(thiscol)))
    else:
        if options.star:
            dir= '../data/star/'
        elif options.nuvx:
            dir= '../data/nuvx/'
        elif options.uvx:
            dir= '../data/uvx/'
        elif options.nuvxall:
            dir= '../data/nuvx_all/'
        elif options.rrlyrae:
            dir= '../data/rrlyrae/'
        #RA, Dec, key
        racol= []
        deccol= []
        if options.rrlyrae: idcol= []
        for qso in qsoKeys:
            #open relevant file
            file= fu.table_fields(os.path.join(dir,qso))
            indx= (file.mjd_g == 0.)
            if True in indx:
                deccol.append(file.dec[indx][0])
                racol.append(file.ra[indx][0])
            else:
                deccol.append(file.dec[0])
                racol.append(file.ra[0])
            if options.rrlyrae: idcol.append(file.id[0])
        cols.append(pyfits.Column(name='dec',format='D',
                                  array=numpy.array(deccol)))
        cols.append(pyfits.Column(name='ra',format='D',
                                  array=numpy.array(racol)))
        if options.rrlyrae:
            cols.append(pyfits.Column(name='ID',format='K',
                                      array=numpy.array(idcol)))
    thiscol= []
    for qso in qsoKeys:
        thiscol.append(qso)
    cols.append(pyfits.Column(name='key',format='30A',
                              array=numpy.array(thiscol)))
    columns= pyfits.ColDefs(cols)
    tbhdu= pyfits.new_table(columns)
    tbhdu.writeto(options.outfilename)
    #create header
    outfile= pyfits.open(options.outfilename)
    hdr= outfile[1].header
    hdr.update('type',type,'type of covariance fit')
    hdr.update('band',type,'band(s) fit')
    hdr.update('mean',type,'type of mean fit')
    return None

def get_options():
    usage = "usage: %prog [options] <filename>\n\nfilename= name of the file that contains the fits"
    parser = OptionParser(usage=usage)
    parser.add_option("-o",dest='outfilename', default=None,
                      help="Name of the fits file that the fits will be dumped to")
    parser.add_option("-i",dest='inputfilename', default='../data/S82qsos.fits',
                      help="Name of the file that holds data on all of the objects (RA,Dec)")
    parser.add_option("--star",action="store_true", dest="star",
                      default=False,
                      help="Fit stars")
    parser.add_option("--nuvx",action="store_true", dest="nuvx",
                      default=False,
                      help="Fit nUVX sample")
    parser.add_option("--uvx",action="store_true", dest="uvx",
                      default=False,
                      help="Fit UVX sample")
    parser.add_option("--nuvxall",action="store_true", dest="nuvxall",
                      default=False,
                      help="Fit nUVX_all sample")
    parser.add_option("--rrlyrae",action="store_true", dest="rrlyrae",
                      default=False,
                      help="Fit RR Lyrae sample")
    return parser

if __name__ == '__main__':
    dumpFits(get_options())
