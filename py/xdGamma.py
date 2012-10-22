#Do
#export PYTHONPATH=~/Repos/extreme-deconvolution-ngerrors/py:$PYTHONPATH 
import sys
import os, os.path
import numpy
import cPickle as pickle
from scipy.cluster.vq import vq, kmeans
from optparse import OptionParser
from extreme_deconvolution import extreme_deconvolution
from skewQSO import _ERASESTR
def xdGamma(parser):
    (options,args)= parser.parse_args()
    if len(args) == 0:
        parser.print_help()
        return
    if options.outfilename is None:
        print "-o filename options needs to be set ..."
        print "Returning ..."
        return None
    if os.path.exists(options.outfilename):
        print options.outfilename+" exists ..."
        print "*Not* overwriting ..."
        print "Remove file before running ..."
        return
    numpy.random.seed(seed=options.seed)
    #Restore samples
    savefilename= args[0]
    print "Reading data ..."
    if os.path.exists(savefilename):
        savefile= open(savefilename,'rb')
        samples= pickle.load(savefile)
        type= pickle.load(savefile)
        band= pickle.load(savefile)
        mean= pickle.load(savefile)
        savefile.close()
    else:
        print "Input file does not exist ..."
        print "Returning ..."
        return
    #Prepare samples for XD
    print "Preparing data ..."
    if type == 'powerlawSF':
        if len(band) > 1:
            print "multi-band not implemented yet ..."
            print "Returning ..."
            return
        else:
            nparams= 1
    elif type == 'DRW':
        print "DRW not implemented yet ..."
        print "Returning ..."
        return
    elif type == 'KS11':
        nparams= 1
    elif type == 'scatter':
        print "scatter not implemented yet ..."
        print "Returning ..."
        return
    ndata= len(samples)
    ydata= numpy.zeros((ndata,nparams))
    ycovar= numpy.zeros((ndata,nparams,nparams))
    ngamp= numpy.zeros((ndata,options.g))
    ngmean= numpy.zeros((ndata,options.g,nparams))
    ngcovar= numpy.zeros((ndata,options.g,nparams,nparams))
    for ii, key in enumerate(samples.keys()):
        sys.stdout.write('\r'+_ERASESTR+'\r')
        sys.stdout.flush()
        sys.stdout.write('\rWorking on preparing %i / %i\r' % (ii+1,ndata))
        sys.stdout.flush()
        if type == 'powerlawSF':
           #Stack as A,g,Ac,gc
            loggammas= []
            for sample in samples[key]:
                loggammas.append(numpy.log(sample['gamma'][0]))
            loggammas= numpy.array(loggammas)
            ydata[ii,0]= numpy.mean(loggammas)
            ycovar[ii,:,:]= numpy.var(loggammas)
            #Fit with g Gaussians
            thisydata= numpy.reshape(loggammas,(len(loggammas),nparams))
            thisycovar= numpy.zeros((len(loggammas),nparams))+numpy.var(loggammas)*10.**-4. #regularize
            thisxamp= numpy.ones(options.g)/options.g
            thisxcovar= numpy.ones((options.g,nparams,nparams))*numpy.var(loggammas)
            thisxmean= numpy.ones((options.g,nparams))*numpy.mean(loggammas)+numpy.std(loggammas)*numpy.random.normal(size=(options.g,nparams))/4.
            #print numpy.mean(loggammas), numpy.std(loggammas)
            extreme_deconvolution(thisydata,thisycovar,thisxamp,thisxmean,thisxcovar)
            ngamp[ii,:]= thisxamp
            ngmean[ii,:,:]= thisxmean
            ngcovar[ii,:,:,:]= thisxcovar
            if len(band) > 1:
                print "Multi-band not supported currently"
                print "Returning ..."
                return
        elif type == 'DRW':
                print "DRW not supported currently"
                print "Returning ..."
                return
        elif type == 'KS11':
            print "type == 'KS11' not implemented yet ..."
            print "Returning ..."
            return
    sys.stdout.write('\r'+_ERASESTR+'\r')
    sys.stdout.flush()
    #Outlier rejection
    #if type == 'powerlawSF':
    #    indx= (ydata[:,0] > -7.21)
    #    ydata= ydata[indx,:]
    #    ycovar= ycovar[indx,:,:]
    #Initial parameters for XD
    print "Running XD ..."
    xamp= numpy.ones(options.k)/float(options.k)
    xmean= numpy.zeros((options.k,nparams))
    for kk in range(options.k):
        xmean[kk,:]= numpy.mean(ydata,axis=0)\
            +numpy.random.normal()*numpy.std(ydata,axis=0)
    xcovar= numpy.zeros((options.k,nparams,nparams))
    for kk in range(options.k):
        xcovar[kk,:,:]= numpy.cov(ydata.T)
    extreme_deconvolution(ydata,ycovar,xamp,xmean,xcovar)
    if True:
        print xamp
        print xmean
        print xcovar
    #Prepare for saving
    print "Preparing output for saving ..."
    #Save
    print "Saving ..."
    if os.path.exists(options.outfilename):
        print options.outfilename+" exists ..."
        print "*Not* overwriting ..."
        print "Remove file before running ..."
        return
    if options.savefits:
        raise NotImplementedError("Fits saving not implemented yet")
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
        pickle.dump(xamp,outfile)
        pickle.dump(xmean,outfile)
        pickle.dump(xcovar,outfile)
        outfile.close()
    return

def get_options():
    usage = "usage: %prog [options] <filename>\n\nfilename= name of the file that contains the fits"
    parser = OptionParser(usage=usage)
    parser.add_option("-o",dest='outfilename', default=None,
                      help="Name of the file that the K-means solution will be saved to")
    parser.add_option("-k",dest='k', default=10,type='int',
                      help="'K' in K-means")
    parser.add_option("-g",dest='g', default=2,type='int',
                      help="Number of Gaussians to fit the samples with")
    parser.add_option("--seed",dest='seed', default=1,type='int',
                      help="seed for random number generator")
    parser.add_option("--savefits",action="store_true", 
                      default=False, dest="savefits",
                      help="Save as a fits file")
    return parser

if __name__ == '__main__':
    xdGamma(get_options())
