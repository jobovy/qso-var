import os, os.path
import numpy
import cPickle as pickle
from optparse import OptionParser
from varqso import VarQso
_DEBUG=True
def QSOfilenames(dir='../data/s82qsos/'):
    return [os.path.join(dir,qso) for qso in os.listdir(dir)]
def fitQSO(parser):
    (options,args)= parser.parse_args()
    if len(args) == 0:
        parser.print_help()
        return
    if len(args) == 2:
        othersavefilename= args[1]
        othersavefile= open(othersavefilename,'rb')
        otherparams= pickle.load(othersavefile)
        othersavefile.close()
    else:
        otherparams= {}
    savefilename= args[0]
    if os.path.exists(savefilename):
        savefile= open(savefilename,'rb')
        params= pickle.load(savefile)
        type= pickle.load(savefile)
        band= pickle.load(savefile)
        mean= pickle.load(savefile)
        savefile.close()
        if params.has_key('.fit'): params.pop('.fit')
        for key in params.keys():
            if (params[key].has_key('gamma') and \
                    (params[key]['gamma'] < 0. or params[key]['gamma'] > 2.)) \
                    or (params[key].has_key('gammagr') \
                            and (params[key]['gammagr'] < 0. \
                                     or params[key]['gammagr'] > 2.)):
                print "Popping bad gamma ..."
                params.pop(key)
    else:
        params= {}
        type= options.type
        mean= options.mean
        band= options.band
    if options.star:
        dir= '../data/star/'
    elif options.nuvx:
        dir= '../data/nuvx/'
    elif options.nuvxall:
        dir= '../data/nuvx_all/'
    elif options.uvx:
        dir= '../data/uvx/'
    elif options.rrlyrae:
        dir= '../data/rrlyrae/'
    else:
        dir= '../data/s82qsos/'
    if options.resampled:
        if os.path.exists(options.infile):
            samplefile= open(options.infile,'rb')
            qsos= pickle.load(samplefile)
            samplefile.close()
        else:
            print "'--resampled' is set, but -i filename does not exist ..."
            print "Returning ..."
            return None
    else:
        qsos= QSOfilenames(dir=dir)
    savecount= 0
    count= len(params)
    for qso in qsos:
        if options.resampled:
            key= qso[0]
        else:
            key= os.path.basename(qso)
        if params.has_key(key) or otherparams.has_key(key):
            continue
        try:
            if int(key[5:7]) != options.rah and options.rah != -1:
                continue
        except ValueError:
            if options.rah == -2 or options.rah == -1:
                pass
            else:
                print "Skipping ValueError "+key
                continue
        print "Working on "+str(count)+": "+key
        if options.resampled:
            v= qso[1]
        else:
            v= VarQso(qso,flux=options.fitflux)
        if v.nepochs(band) < 20:
            print "This object does not have enough epochs ..."
            continue
        params[key]= v.fit(band=band,type=type,loglike=True,mean=mean)
        if _DEBUG:
            print params[key]
        if params[key]['loglike'] == -numpy.finfo(numpy.dtype(numpy.float64)).max:
            print "Popping bad fit ..."
            params.pop(key)
        if savecount == options.saveevery:
            print "Saving ..."
            save_pickles(params,type,band,mean,savefilename)
            savecount= 0
        savecount+= 1
        count+= 1
    save_pickles(params,type,band,mean,savefilename)
    print "All done"

def save_pickles(params,type,band,mean,savefilename):
    saving= True
    interrupted= False
    tmp_savefilename= savefilename+'.tmp'
    while saving:
        try:
            savefile= open(tmp_savefilename,'wb')
            pickle.dump(params,savefile)
            pickle.dump(type,savefile)
            pickle.dump(band,savefile)
            pickle.dump(mean,savefile)
            savefile.close()
            os.rename(tmp_savefilename,savefilename)
            saving= False
            if interrupted:
                raise KeyboardInterrupt
        except KeyboardInterrupt:
            if not saving:
                raise
            print "KeyboardInterrupt ignored while saving pickle ..."
            interrupted= True

def get_options():
    usage = "usage: %prog [options] <savefilename>\n\nsavefilename= name of the file that the fits will be saved to"
    parser = OptionParser(usage=usage)
    parser.add_option("-b","--band",dest='band',default='r',
                      help="band(s) to fit")
    parser.add_option("-t","--type",dest='type',default='powerlawSF',
                      help="Type of model to fit (powerlawSF, powerlawSFratios, or DRW)")
    parser.add_option("--mean",dest='mean',default='zero',
                      help="Type of mean to fit (zero, const)")
    parser.add_option("--saveevery",dest='saveevery',type='int',
                      default=100,
                      help="Save every --saveevery iterations")
    parser.add_option("--rah",dest='rah',type='int',
                      default=-1,
                      help="RA hour to consider (-1: all, -2: ValueError)")
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
    parser.add_option("--resampled",action="store_true", dest="resampled",
                      default=False,
                      help="Objects are 'resampled': stored in sav-file")
    parser.add_option("-i","--infile",dest='infile',
                      default=None,
                      help="Input file if --resampled")
    parser.add_option("--fitflux",action="store_true", dest="fitflux",
                      default=False,
                      help="Fit fluxes rather than magnitudes")
    return parser

if __name__ == '__main__':
    fitQSO(get_options())
