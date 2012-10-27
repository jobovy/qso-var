import os, os.path
import cPickle as pickle
from optparse import OptionParser
import numpy
from galpy.util import bovy_plot
from matplotlib import pyplot
from varqso import _load_fits
def plotSkewvs(parser):
    (options,args)= parser.parse_args()
    if len(args) == 0:
        parser.print_help()
        return
    savefilename= args[0]
    if os.path.exists(savefilename):
        savefile= open(savefilename,'rb')
        skews= pickle.load(savefile)
        gaussskews= pickle.load(savefile)
        type= pickle.load(savefile)
        band= pickle.load(savefile)
        mean= pickle.load(savefile)
        taus= pickle.load(savefile)      
        savefile.close()
    else:
        parser.print_help()
        return
    #Load additional catalogs
    if not options.fitsfile is None:
        fitsfile= open(options.fitsfile,'rb')
        params= pickle.load(fitsfile)
        fitsfile.close()
    #Match
    s82qsos= _load_fits('../data/S82qsos.fits')
    shenqsos= _load_fits('../data/dr7_bh_May09_2011-woname.fits')
    s82qsoDict= {}
    shenqsoDict= {}
    shenindexing= numpy.arange(len(shenqsos),dtype='int')
    ii=0
    shenmatches= 0
    for qso in s82qsos:
        s82qsoDict[qso.oname.strip().replace(' ', '')+'.fit']= ii
        ii+= 1       
        continue
        #Also match to Shen
        shenindx= (numpy.fabs(shenqsos.redshift-qso.z) < 10.**-4.)\
            *(numpy.fabs(shenqsos.ra-qso.ra) < 10.**-4.)
        if numpy.sum(shenindx) == 0: 
            shenqsoDict[qso.oname.strip().replace(' ', '')+'.fit']= -1
        elif numpy.sum(shenindx) > 1: 
            print "two solutions:"
            print qso.oname.strip(), shenqsos[shenindx].redshift, shenqsos[shenindx].ra
        else:
            shenqsoDict[qso.oname.strip().replace(' ', '')+'.fit']= shenindexing[shenindx]
            shenmatches+=1
        ii+= 1       
#    print "Found %i matches in Shen out of %i" % (shenmatches,len(shenqsos))
    ii=0
    for qso in shenqsos:
        if not qso.oname == '':
            shenqsoDict[qso.oname.strip().replace(' ', '')+'.fit']= ii
        else:
            shenqsoDict[qso.oname.strip().replace(' ', '')+'.fit']= -1
        ii+= 1       
    #Accumulate
    if options.type == 'z':
        lookin = 's82'
        quant= 'z'
        xlabel=r'$\mathrm{redshift}$'
        xrange=[0.,3.2]
        nbins= 16
        quants= numpy.linspace(0.3,3.,nbins)
    elif options.type == 'ledd':
        lookin = 'shen'
        quant= 'logedd_ratio'
        xlabel=r'$\log L / L_\mathrm{edd}$'
        nbins= 16
        quants= numpy.linspace(-2.,0.,nbins)
        xrange=[-2.,0.]
    elif options.type.lower() == 'mi':
        lookin = 'shen'
        quant= 'mi_z2'
        xlabel=r'$M_i(z=2)$'
        nbins= 16
        quants= numpy.linspace(-30.,-22.,nbins)
        xrange=[-22.,-29.]
    elif options.type.lower() == 'loga':
        lookin = 'fits'
        quant= 'logA'
        xlabel=r'$\log A\ (\mathrm{variability\ at\ 1\ yr}$'
        nbins= 16
        quants= numpy.linspace(-7.,-1.,nbins)#2x
        xrange=[-4.,0.]
    elif options.type.lower() == 'gamma':
        lookin = 'fits'
        quant= 'gamma'
        xlabel=r'$\gamma\ (\mathrm{power-law\ index}$'
        nbins= 16
        quants= numpy.linspace(0.,1.,nbins)#2x
        xrange=[0.,.5]
    tauindx= 5
    keys= skews.keys()
    allskews= numpy.zeros(nbins-1)
    allgaussskews= numpy.zeros((2,nbins-1))
    #Statistic
    statistic= numpy.median
    q= 0.95
    for jj in range(nbins-1):
        thisskews= []
        thisgaussskews= []
        for ii, key in enumerate(keys):
            if lookin == 's82':
                if s82qsos[s82qsoDict[key]][quant] < quants[jj] \
                        or s82qsos[s82qsoDict[key]][quant] >= quants[jj+1]:
                    continue
            elif lookin == 'fits' and quant == 'logA':
                thisA= params[key][quant]+params[key]['gamma']*numpy.log(1.+s82qsos[s82qsoDict[key]]['z'])
                if thisA < quants[jj] \
                        or thisA >= quants[jj+1]:
                    continue
            elif looking == 'fits':
                if params[key][quant] < quants[jj] \
                        or params[key][quant] >= quants[jj+1]:
                    continue
            else:
                if shenqsoDict[key] == -1:
                    continue
                if shenqsos[shenqsoDict[key]][quant] < quants[jj] \
                        or shenqsos[shenqsoDict[key]][quant] >= quants[jj+1]:
                    continue
            thisskews.append(-skews[key][tauindx]) #go to regular definition
            thisgaussskews.append(-gaussskews[key][:,tauindx]) #go to regular definition
        if len(thisskews) < 2:
            allskews[jj]= numpy.nan
            allgaussskews[:,jj]= numpy.nan
            continue
        thisskews= numpy.array(thisskews)
        thisgaussskews= numpy.array(thisgaussskews)
        #Fix nans
        indx= numpy.isnan(thisskews)
        thisskews= thisskews[True-indx]
        thisgaussskews= thisgaussskews[True-indx,:]
        #Save
        allskews[jj]= statistic(thisskews)
        mediangaussskew= statistic(thisgaussskews,axis=0)
        allgaussskews[:,jj]= quantile(mediangaussskew,q=q)
    quants+= (quants[1]-quants[0])/2.
    quants= quants[0:-1]
    if options.type == 'logA' or options.type == 'gamma':
        quants/= 2.
    #Plot
    bovy_plot.bovy_print(fig_width=7.)
    bovy_plot.bovy_plot(quants,allskews,'ko',
                        xlabel=xlabel,
                        ylabel=r'$\mathrm{skew}(\tau = %i\ \mathrm{days})$' % (int(taus[tauindx]*365.25)),
                        zorder=5,
                        yrange=[-.6,.6],
                        xrange=xrange)
    pyplot.fill_between(quants,allgaussskews[0,:],allgaussskews[1,:],
                        color='0.75',zorder=0)
    bovy_plot.bovy_end_print(options.plotfilename)
    return None

def quantile(t,q=0.68):
    """calculate the quantile of a distribution for this problem,
    t= ngauss,nt, quantile compute for each t"""
    out= numpy.zeros(2)
    sortedt= sorted(t)
    low= sortedt[int(numpy.floor((0.5-q/2.)*len(t)))]
    high= sortedt[int(numpy.floor((0.5+q/2.)*len(t)))]
    out[0]= low
    out[1]= high
    return out

def get_options():
    usage = "usage: %prog [options] <savefilename>\n\nsavefilename= name of the file that holds the skews"
    parser = OptionParser(usage=usage)
    parser.add_option("-o",dest='plotfilename',default=None,
                      help="Name for plotfile")
    parser.add_option("-i","--indx",dest='indx',default=None,type='float',
                      help="Just plot this index")
    parser.add_option("-b","--band",dest='band',default='r',
                      help="band(s) to sample")
    parser.add_option("-t","--type",dest='type',default='z',
                      help="Type of plot to make")
    parser.add_option("--mean",dest='mean',default='zero',
                      help="Type of mean to sample (zero, const)")
    parser.add_option("-f","--fitsfile",dest='fitsfile',
                      default=None,
                      help="File that holds the best-fits")
    return parser
    
if __name__ == '__main__':
    plotSkewvs(get_options())
