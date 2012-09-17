import os, os.path
import numpy
import cPickle as pickle
from scipy.cluster.vq import vq, kmeans
from optparse import OptionParser
def kmeansFits(parser):
    (options,args)= parser.parse_args()
    if len(args) == 0:
        parser.print_help()
        return
    if options.outfilename is None:
        print "-o filename options needs to be set ..."
        print "Returning ..."
        return None
    numpy.random.seed(seed=options.seed)
    #Restore fits
    savefilename= args[0]
    if os.path.exists(savefilename):
        savefile= open(savefilename,'rb')
        params= pickle.load(savefile)
        type= pickle.load(savefile)
        band= pickle.load(savefile)
        savefile.close()
    else:
        print "Input file does not exist ..."
        print "Returning ..."
        return
    #Prepare params for K-means
    print "Preparing data ..."
    if type == 'powerlawSF':
        if len(band) > 1:
            nparams= 4
        else:
            nparams= 2
    elif type == 'DRW':
        if len(band) == 1:
            nparams= 2
        else:
            print "DRW for multi-band fits not implemented yet ..."
            print "Returning ..."
            return
    elif type == 'KS11':
        nparams= 3
    elif type == 'scatter':
        nparams= 1
    ndata= len(params)
    kIn= numpy.zeros((ndata,nparams))
    if type == 'powerlawSF':
        #Stack as A,g,Ac,gc
        kIn[:,0]= numpy.array([p['logA'] for p in params.values()]).reshape(ndata)
        kIn[:,1]= numpy.array([p['gamma'] for p in params.values()]).reshape(ndata)
        if len(band) > 1:
            kIn[:,2]= numpy.array([p['logAgr'] for p in params.values()]).reshape(ndata)
            kIn[:,3]= numpy.array([p['gammagr'] for p in params.values()]).reshape(ndata)
    elif type == 'DRW':
        print "type == 'DRW' not implemented yet ..."
        print "Returning ..."
        return
    elif type == 'KS11':
        #Stack as A,g,s
        kIn[:,0]= numpy.array([p['logA'] for p in params.values()]).reshape(ndata)
        kIn[:,1]= numpy.array([p['gamma'] for p in params.values()]).reshape(ndata)
        kIn[:,2]= numpy.array([p['s'] for p in params.values()]).reshape(ndata)
    #Whiten, i.e., give unit variance
    print "Whitening data ..."
    whitenFactors= numpy.zeros(nparams)
    for ii in range(nparams):
        whitenFactors[ii]= numpy.std(kIn[:,ii])
        kIn[:,ii]/= whitenFactors[ii]
    #Ready to run K-means
    print "Running K-means ..."
    book, dist= kmeans(kIn,options.k)
    assign, dist= vq(kIn,book)
    #De-whiten the codebook
    for ii in range(nparams):
        book[:,ii]*= whitenFactors[ii]
    #Prepare for saving
    print "Preparing output for saving ..."
    outparams= []
    weights= []
    for kk in range(options.k):
        if type == 'powerlawSF':
            if len(band) > 1:
                outparams.append({'logA':book[kk,0],
                                  'gamma':book[kk,1],
                                  'logAgr':book[kk,2],
                                  'gammagr':book[kk,3]})
            else:
                outparams.append({'logA':book[kk,0],
                                  'gamma':book[kk,1]})
        elif type == 'DRW':
            print "DRW not implemented yet ..."
            print "Returning ..."
            return
        if type == 'KS11':
                outparams.append({'logA':book[kk,0],
                                  'gamma':book[kk,1],
                                  's':book[kk,2]})
        thisassign= assign[(assign == kk)]
        weights.append(len(thisassign))
    #Save
    print "Saving ..."
    if os.path.exists(options.outfilename):
        print options.outfilename+" exists ..."
        print "*Not* overwriting ..."
        print "Remove file before running ..."
        return
    if options.savefits:
        import pyfits
        cols= []
        if type == 'powerlawSF':
            colA= []
            colg= []
            for kk in range(options.k):
                colA.append(outparams[kk]['logA'])
                colg.append(outparams[kk]['gamma'])
            colA= numpy.array(colA)
            colg= numpy.array(colg)
            colw= numpy.log(numpy.array(weights))
            cols.append(pyfits.Column(name='logA',format='E',
                                      array=colA))
            cols.append(pyfits.Column(name='gamma',format='E',
                                      array=colg))
        elif type == 'KS11':
            colA= []
            colg= []
            cols= []
            for kk in range(options.k):
                colA.append(outparams[kk]['logA'])
                colg.append(outparams[kk]['gamma'])
                colg.append(outparams[kk]['s'])
            colA= numpy.array(colA)
            colg= numpy.array(colg)
            cols= numpy.array(colg)
            cols.append(pyfits.Column(name='logA',format='E',
                                      array=colA))
            cols.append(pyfits.Column(name='gamma',format='E',
                                      array=colg))
            cols.append(pyfits.Column(name='s',format='E',
                                      array=cols))           
        colw= numpy.log(numpy.array(weights))
        cols.append(pyfits.Column(name='logweight',format='E',
                                  array=colw))
        columns= pyfits.ColDefs(cols)
        tbhdu= pyfits.new_table(columns)
        tbhdu.writeto(options.outfilename)
    else:
        outfile= open(options.outfilename,'wb')
        pickle.dump(outparams,outfile)
        pickle.dump(weights,outfile)
        outfile.close()
    return

def get_options():
    usage = "usage: %prog [options] <filename>\n\nfilename= name of the file that contains the fits"
    parser = OptionParser(usage=usage)
    parser.add_option("-o",dest='outfilename', default=None,
                      help="Name of the file that the K-means solution will be saved to")
    parser.add_option("-k",dest='k', default=10,type='int',
                      help="'K' in K-means")
    parser.add_option("--seed",dest='seed', default=1,type='int',
                      help="seed for random number generator")
    parser.add_option("--savefits",action="store_true", 
                      default=False, dest="savefits",
                      help="Save as a fits file")
    return parser

if __name__ == '__main__':
    kmeansFits(get_options())
