import os, os.path
import cPickle as pickle
import numpy as nu
from optparse import OptionParser
from varqso import VarQso, panstarrs_sampling, sdss_sampling
from fitQSO import QSOfilenames
def resampleObjs(parser):
    (options,args)= parser.parse_args()
    if len(args) == 0:
        parser.print_help()
        return
    if os.path.exists(args[0]):
        print filename+" exists"
        print "Remove this file before running ..."
        print "Returning ..."
        return None
    if options.fitsfile is None:
        print "-f or --fitsfile must be set to the file containing the fits ..."
        print "Returning ..."
        return None
    #Load location of the data
    if options.sample == 'nuvx':
        dir= '../data/nuvx/'
    elif options.sample == 'qso':
        dir= '../data/s82qsos/'
    objs= QSOfilenames(dir=dir)
    #Load the fits
    if os.path.exists(options.fitsfile):
        fitsfile= open(options.fitsfile,'rb')
        params= pickle.load(fitsfile)
        type= pickle.load(fitsfile)
        band= pickle.load(fitsfile)
        fitsfile.close()
    else:
        print options.fitsfile+" does not exist ..."
        print "Returning ..."
        return None
    #Load sampling
    if options.sampling == 'PS1':
        sampling= panstarrs_sampling(3,startmjd=2.*365.25)
    elif options.sampling == 'SDSS-PS1':
        sampling= sdss_sampling(startmjd=-2.*365.25)
        sampling.extend(panstarrs_sampling(1,startmjd=2.*365.25))
    elif options.sampling == 'SDSS': pass
    else:
        print "Input to --sampling not understood ..."
        print "Returning ..."
        return None
    #Re-sample each source
    out= []
    savecount, count= 0, 0
    for obj in objs:
        key= os.path.basename(obj)
        print "Working on "+str(count)+" ("+str(savecount)+"): "+key
        savecount+= 1
        v= VarQso(obj)
        #Find fit
        try:
            thisfit= params[key]
        except KeyError:
            print "Fit not found, skipping this object ..."
            nepochs= v.nepochs(band)
            if nepochs < 20:
                print "Because #epochs < 20 ..."
            continue
        #Set LC model
        if options.nocolorvar:
            thisfit['gammagr']= 0.
            if 'logAgr' in thisfit.keys():
                thisfit['logAgr']= -7.
            else:
                thisfit['logAri']= -7.
            if options.sampling == 'SDSS':
                sampling= v.mjd['g']
                sampling= [(s,'g') for s in sampling]
            v.setLCmodel(thisfit,band,type)
            indx= v.mjd_overlap(band=band)
            refband= 'r'
            try:
                o= v.resample([(mjd,refband) for mjd in v.mjd[refband][indx[refband]]],
                              band=refband,errors=False)
            except nu.linalg.LinAlgError:
                print thisfit
                continue
            xs= []
            ys= []
            errs= []
            for b in band:
                #Add errors to the underlying
                xs.extend([(mjd,b) for mjd in v.mjd[b][indx[b]]])
                for ii in range(len(v.mjd[b][indx[b]])):
                    ys.append(o.m[refband][ii]+
                              nu.random.randn()*v.err_m[b][indx[b]][ii])
                    errs.append(v.err_m[b][indx[b]][ii])
            out.append([key,VarQso(xs,ys,errs,band=band,medianize=False)])
        else:
            v.setLCmodel(thisfit,band,type)
            #Resample
            out.append([key,v.resample(sampling,band=band)])
        count+= 1
    #Save
    outfile= open(args[0],'wb')
    pickle.dump(out,outfile)
    pickle.dump(band,outfile)
    outfile.close()
    return None

def get_options():
    usage = "usage: %prog [options] <savefilename>\n\nsavefilename= name of the file that the resampled lightcurves will be saved to"
    parser = OptionParser(usage=usage)
    parser.add_option("-f","--fitsfile",dest='fitsfile',default=None,
                      help="File that contains the fits")
    parser.add_option("--sample",dest='sample',default='nuvx',
                      help="sample to re-sample ('nuvx','qso')")
    parser.add_option("--sampling",dest='sampling',default='PS1',
                      help="kind of resampling to perform ('SDSS','PS1','SDSS-PS1')")
    parser.add_option("--nocolorvar",action="store_true", dest="nocolorvar",
                      default=False,
                      help="Remove color-variability before resampling, assume you are doing g, otherwise first band")
    return parser

if __name__ == '__main__':
    resampleObjs(get_options())
