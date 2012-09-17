import os, os.path
import cPickle as pickle
from optparse import OptionParser
import numpy as nu
from varqso import _load_fits
import galpy.util.bovy_plot as bovy_plot
def plotFits(parser):
    (options,args)= parser.parse_args()
    if len(args) == 0:
        parser.print_help()
        return
    params= []
    for filename in args:
        if os.path.exists(filename):
            savefile= open(filename,'rb')  
            if options.kmeansOut:
                thisparams= pickle.load(savefile)
                thisDict= {}
                for kk in range(len(thisparams)):
                    thisDict[str(kk)]= thisparams[kk]
                params.append(thisDict)
            else:
                params.append(pickle.load(savefile))
            savefile.close()
        else:
            print filename+" does not exist ..."
            print "Returning ..."
            return
    if options.plottype == 'Ag':
        ys= nu.array([p['gamma'] for p in params[0].values()]).reshape(len(params[0]))
        if options.type == 'powerlawSFratios':
            xs= nu.array([p['logAr'] for p in params[0].values()]).reshape(len(params[0]))/2.
        else:
            xs= nu.array([p['logA'] for p in params[0].values()]).reshape(len(params[0]))/2.
        xrange=[-9.21/2.,0.]
        yrange=[0.,1.25]
        xlabel= r'$\log A_{'+options.band+r'}\ \mathrm{(amplitude\ at\ 1\ yr)}$'
        ylabel= r'$\gamma_{'+options.band+r'}\ \mathrm{(power\ law\ exponent)}$'
    elif options.plottype == 'As':
        ys= nu.array([p['s']-1. for p in params[0].values()]).reshape(len(params[0]))
        if options.type == 'powerlawSFratios':
            xs= nu.array([p['logAr'] for p in params[0].values()]).reshape(len(params[0]))/2.
        else:
            xs= nu.array([p['logA'] for p in params[0].values()]).reshape(len(params[0]))/2.
        xrange=[-9.21/2.,0.]
        yrange=[-0.9,0.4]
        xlabel= r'$\log A_{'+options.band+r'}\ \mathrm{(amplitude\ at\ 1\ yr)}$'
        ylabel= r'$s_{'+options.band+r'}$'
    if options.plottype == 'gs':
        ys= nu.array([p['s']-1. for p in params[0].values()]).reshape(len(params[0]))
        xs= nu.array([p['gamma'] for p in params[0].values()]).reshape(len(params[0]))/2.
        xrange=[0.,1.2]
        yrange=[-0.9,0.4]
        xlabel= r'$\gamma_{'+options.band+r'}\ \mathrm{(power-law\ index)}$'
        ylabel= r'$s_{'+options.band+r'}$'
    elif options.plottype == 'AA':
        ys= []
        xs= []
        for key in params[0].keys():
            try:
                ys.append(params[1][key]['logA']/2.)
                xs.append(params[0][key]['logA']/2.)
            except KeyError:
                continue
        xs= nu.array(xs).reshape(len(xs))
        ys= nu.array(ys).reshape(len(xs))
        xrange=[-9.21/2.,0.]
        yrange=xrange
        xlabel= r'$\log A_'+options.band[0]+r'\ \mathrm{(amplitude\ at\ 1\ yr)}$'
        ylabel= r'$\log A_'+options.band[1]+r'\ \mathrm{(amplitude\ at\ 1\ yr)}$'
        xlabelnounits= r'$\log A_'+options.band[0]+r'$'
        ylabelnounits= r'$\log A_'+options.band[1]+r'$'
    elif options.plottype == 'AAz':
        qsos= open_qsos()
        qsoDict= {}
        ii=0
        for qso in qsos:
            qsoDict[qso.oname.strip().replace(' ', '')+'.fit']= ii
            ii+= 1
        ys= []
        xs= []
        for key in params[0].keys():            
            try:
                xtmp= qsos[qsoDict[key]].z
            except KeyError:
                continue
            try:
                ytmp=params[1][key]['logA']/2.
            except KeyError:
                continue
            ys.append(params[0][key]['logA']/2-ytmp)
            xs.append(xtmp)
        xs= nu.array(xs).reshape(len(xs))
        ys= nu.array(ys).reshape(len(xs))
        yrange=[-1.,1.]
        xrange=[0.,5.5]
        ylabel= r'$\log A_r/A_g\ \mathrm{(amplitudes\ at\ 1\ yr)}$'
        xlabel= r'$z\ (\mathrm{redshift})$'
    elif options.plottype == 'gg':
        ys= []
        xs= []
        for key in params[0].keys():
            try:
                ys.append(params[1][key]['gamma'])
                xs.append(params[0][key]['gamma'])
            except KeyError:
                continue
        xs= nu.array(xs).reshape(len(xs))
        ys= nu.array(ys).reshape(len(xs))
        xrange=[0.,1.15]
        yrange=xrange
        xlabel= r'$\gamma_'+options.band[0]+r'\ \mathrm{(power\ law\ exponent)}$'
        ylabel= r'$\gamma_'+options.band[1]+r'\ \mathrm{(power\ law\ exponent)}$'
        xlabelnounits= r'$\gamma_'+options.band[0]+r'$'
        ylabelnounits= r'$\gamma_'+options.band[1]+r'$'
    elif options.plottype == 'Acgc':
        ys= nu.array([p['gammagr'] for p in params[0].values()]).reshape(len(params[0]))
        if options.type == 'powerlawSFratios':
            xs= nu.array([p['logAr'] for p in params[0].values()]).reshape(len(params[0]))/2.
        else:
            xs= nu.array([p['logAgr'] for p in params[0].values()]).reshape(len(params[0]))/2.
        xrange=[-9.21/2.,0.]
        yrange=[0.,1.25]
        xlabel= r'$\log A^c_{'+options.band+r'}\ \mathrm{(amplitude\ at\ 1\ yr)}$'
        ylabel= r'$\gamma^c_{'+options.band+r'}\ \mathrm{(power\ law\ exponent)}$'
    elif options.plottype == 'Acg':
        ys= nu.array([p['gamma'] for p in params[0].values()]).reshape(len(params[0]))
        if options.type == 'powerlawSFratios':
            xs= nu.array([p['logAr'] for p in params[0].values()]).reshape(len(params[0]))/2.
        else:
            xs= nu.array([p['logAgr'] for p in params[0].values()]).reshape(len(params[0]))/2.
        xrange=[-9.21/2.,0.]
        yrange=[0.,1.25]
        xlabel= r'$\log A^c_{'+options.band+r'}\ \mathrm{(amplitude\ at\ 1\ yr)}$'
        ylabel= r'$\gamma_{'+options.band+r'}\ \mathrm{(power\ law\ exponent)}$'
    elif options.plottype == 'AcAc':
        ys= []
        xs= []
        for key in params[0].keys():
            try:
                ys.append(params[1][key]['logAgr']/2.)
                xs.append(params[0][key]['logAgr']/2.)
            except KeyError:
                continue
        xs= nu.array(xs).reshape(len(xs))
        ys= nu.array(ys).reshape(len(xs))
        xrange=[-9.21/2.,0.]
        yrange=xrange
        xlabel= r'$\log A^c_'+options.band[0]+r'\ \mathrm{(amplitude\ at\ 1\ yr)}$'
        ylabel= r'$\log A^c_'+options.band[1]+r'\ \mathrm{(amplitude\ at\ 1\ yr)}$'
        xlabelnounits= r'$\log A^c_'+options.band[0]+r'$'
        ylabelnounits= r'$\log A^c_'+options.band[1]+r'$'
    elif options.plottype == 'gcgc':
        ys= []
        xs= []
        for key in params[0].keys():
            try:
                ys.append(params[1][key]['gammagr'])
                xs.append(params[0][key]['gammagr'])
            except KeyError:
                continue
        xs= nu.array(xs).reshape(len(xs))
        ys= nu.array(ys).reshape(len(xs))
        xrange=[0.,1.15]
        yrange=xrange
        xlabel= r'$\gamma^c_'+options.band[0]+r'\ \mathrm{(power\ law\ exponent)}$'
        ylabel= r'$\gamma^c_'+options.band[1]+r'\ \mathrm{(power\ law\ exponent)}$'
        xlabelnounits= r'$\gamma^c_'+options.band[0]+r'$'
        ylabelnounits= r'$\gamma^c_'+options.band[1]+r'$'
    elif options.plottype == 'AAc':
        if options.type == 'powerlawSFratios':
            xs= nu.array([p['logAr'] for p in params[0].values()]).reshape(len(params[0]))/2.
            ys= nu.array([p['logAri'] for p in params[0].values()]).reshape(len(params[0]))/2.
        else:
            xs= nu.array([p['logA'] for p in params[0].values()]).reshape(len(params[0]))/2.
            ys= nu.array([p['logAgr'] for p in params[0].values()]).reshape(len(params[0]))/2.
        xrange=[-9.21/2.,0.]
        yrange=[-9.21/2.,0.]
        xlabel= r'$\log A_{'+options.band+r'}\ \mathrm{(amplitude\ at\ 1\ yr)}$'
        ylabel= r'$\log A^c_{'+options.band+r'}\ \mathrm{(amplitude\ at\ 1\ yr)}$'
    elif options.plottype == 'ggc':
        if len(params) == 2:
            xs, ys= [], []
            for key in params[0].keys():
                try:
                    ys.append(params[1][key]['gammagr'])
                    xs.append(params[0][key]['gamma'])
                except KeyError:
                    continue
            xs= nu.array(xs).reshape(len(xs))
            ys= nu.array(ys).reshape(len(ys))           
        else:
            xs= nu.array([p['gamma'] for p in params[0].values()]).reshape(len(params[0]))
            ys= nu.array([p['gammagr'] for p in params[0].values()]).reshape(len(params[0]))/2.
        xrange=[0.,1.25]
        yrange=xrange
        xlabel= r'$\gamma_{'+options.band+r'}\ \mathrm{(power\ law\ exponent)}$'
        ylabel= r'$\gamma^c_{'+options.band+r'}\ \mathrm{(power\ law\ exponent)}$'
    elif options.plottype == 'Ai':
        qsos= open_qsos()
        qsoDict= {}
        ii=0
        for qso in qsos:
            qsoDict[qso.oname.strip().replace(' ', '')+'.fit']= ii
            ii+= 1
        ys= []
        xs= []
        for key in params[0].keys():
            try:
                xs.append(qsos[qsoDict[key]].mags[3])
            except KeyError:
                continue
            if options.type == 'powerlawSFratios':
                ys.append(params[0][key]['logAr'])
            else:
                ys.append(params[0][key]['logA'])
        ys= nu.array(ys)/2.
        xs= nu.array(xs)
        yrange=[-9.21/2.,0.]
        xrange=[15.0,21.3]
        ylabel= r'$\log A_{'+options.band+r'}\ \mathrm{(amplitude\ at\ 1\ yr)}$'
        xlabel= r'$i_0\ [\mathrm{mag}]$'
    elif options.plottype == 'gi':
        qsos= open_qsos()
        qsoDict= {}
        ii=0
        for qso in qsos:
            qsoDict[qso.oname.strip().replace(' ', '')+'.fit']= ii
            ii+= 1
        ys= []
        xs= []
        for key in params[0].keys():
            ys.append(params[0][key]['gamma'])
            xs.append(qsos[qsoDict[key]].mags[3])
        ys= nu.array(ys)
        xs= nu.array(xs)
        yrange=[0.,1.25]
        xrange=[15.0,21.3]
        ylabel= r'$\gamma_{'+options.band+r'}\ \mathrm{(power\ law\ exponent)}$'
        xlabel= r'$i_0\ [\mathrm{mag}]$'
    elif options.plottype == 'Az':
        qsos= open_qsos()
        qsoDict= {}
        ii=0
        for qso in qsos:
            qsoDict[qso.oname.strip().replace(' ', '')+'.fit']= ii
            ii+= 1
        ys= []
        xs= []
        for key in params[0].keys():            
            try:
                xs.append(qsos[qsoDict[key]].z)
            except KeyError:
                continue
            if options.type == 'powerlawSFratios':
                ys.append(params[0][key]['logAr'])
            else:
                ys.append(params[0][key]['logA'])
        ys= nu.array(ys)/2.
        xs= nu.array(xs)
        yrange=[-9.21/2.,0.]
        xrange=[0.,5.5]
        ylabel= r'$\log A_{'+options.band+r'}\ \mathrm{(amplitude\ at\ 1\ yr)}$'
        xlabel= r'$z\ (\mathrm{redshift})$'
    elif options.plottype == 'sz': #for KS11
        qsos= open_qsos()
        qsoDict= {}
        ii=0
        for qso in qsos:
            qsoDict[qso.oname.strip().replace(' ', '')+'.fit']= ii
            ii+= 1
        ys= []
        xs= []
        for key in params[0].keys():            
            try:
                xs.append(qsos[qsoDict[key]].z)
            except KeyError:
                continue
            if len(params) == 2:
                try:
                    ys.append(params[0][key]['s']/params[1][key]['s'])
                except KeyError:
                    xs.pop()
                    continue
            else:
                ys.append(params[0][key]['s'])
        ys= nu.array(ys)-1.
        xs= nu.array(xs)
        yrange=[-1.1,0.4]
        xrange=[0.,5.5]
        ylabel= r'$s_{'+options.band+r'}$'
        xlabel= r'$z\ (\mathrm{redshift})$'
    elif options.plottype == 'scatterz': #for KS11
        qsos= open_qsos()
        qsoDict= {}
        ii=0
        for qso in qsos:
            qsoDict[qso.oname.strip().replace(' ', '')+'.fit']= ii
            ii+= 1
        ys= []
        xs= []
        for key in params[0].keys():            
            try:
                xs.append(qsos[qsoDict[key]].z)
            except KeyError:
                continue
            ys.append(params[0][key]['logsigma'])
        ys= nu.array(ys)-1.
        xs= nu.array(xs)
        yrange=[-10,0.]
        xrange=[0.,5.5]
        ylabel= r'$\sigma^2_{'+options.band+r'}$'
        xlabel= r'$z\ (\mathrm{redshift})$'
    elif options.plottype == 'ss':
        ys, xs= [], []
        if options.band2 is None:
            for key in params[0].keys():
                try:
                    ys.append(params[1][key]['s'][0]-1.)
                except KeyError:
                    continue
                xs.append(params[0][key]['s'][0]-1.)
        else:
            for key in params[0].keys():
                try:
                    ys.append(params[1][key]['s'+options.band2][0]-1.)
                except KeyError:
                    continue
                xs.append(params[0][key]['s'+options.band][0]-1.)
        xs= nu.array(xs)
        ys= nu.array(ys)
        xrange=[-0.9,0.4]
        yrange=[-0.9,0.4]
        if options.band2 is None:
            xlabel= r'$s$'
            ylabel= r'$s$'
        else:
            xlabel= r'$s_{'+options.band+r'}$'
            ylabel= r'$s_{'+options.band2+r'}$'
        ylabelnounits= ylabel
        xlabelnounits= xlabel
    elif options.plottype == 'gz':
        qsos= open_qsos()
        qsoDict= {}
        ii=0
        for qso in qsos:
            qsoDict[qso.oname.strip().replace(' ', '')+'.fit']= ii
            ii+= 1
        ys= []
        xs= []
        for key in params[0].keys():
            ys.append(params[0][key]['gamma'])
            xs.append(qsos[qsoDict[key]].z)
        ys= nu.array(ys)
        xs= nu.array(xs)
        yrange=[0.,1.25]
        xrange=[0.,5.5]
        ylabel= r'$\gamma_{'+options.band+r'}\ \mathrm{(power\ law\ exponent)}$'
        xlabel= r'$z\ (\mathrm{redshift})$'
    elif options.plottype == 'Aci':
        qsos= open_qsos()
        qsoDict= {}
        ii=0
        for qso in qsos:
            qsoDict[qso.oname.strip().replace(' ', '')+'.fit']= ii
            ii+= 1
        ys= []
        xs= []
        for key in params[0].keys():
            ys.append(params[0][key]['logAgr'])
            xs.append(qsos[qsoDict[key]].mags[3])
        ys= nu.array(ys)/2.
        xs= nu.array(xs)
        yrange=[-9.21/2.,0.]
        xrange=[15.0,21.3]
        ylabel= r'$\log A^c_{'+options.band+r'}\ \mathrm{(amplitude\ at\ 1\ yr)}$'
        xlabel= r'$i_0\ [\mathrm{mag}]$'
    elif options.plottype == 'gci':
        qsos= open_qsos()
        qsoDict= {}
        ii=0
        for qso in qsos:
            qsoDict[qso.oname.strip().replace(' ', '')+'.fit']= ii
            ii+= 1
        ys= []
        xs= []
        for key in params[0].keys():
            ys.append(params[0][key]['gammagr'])
            xs.append(qsos[qsoDict[key]].mags[3])
        ys= nu.array(ys)
        xs= nu.array(xs)
        yrange=[0.,1.25]
        xrange=[15.0,21.3]
        ylabel= r'$\gamma^c_{'+options.band+r'}\ \mathrm{(power\ law\ exponent)}$'
        xlabel= r'$i_0\ [\mathrm{mag}]$'
    elif options.plottype == 'loglike2':
        ys= []
        xs= []
        for key in params[0].keys():
            if not params[1].has_key(key): continue
            ys.append(params[1][key]['loglike'])
            xs.append(params[0][key]['loglike'])
        ys= nu.array(ys)
        xs= nu.array(xs)
        yrange=[0.,200.]
        xrange=[0.,200.]
        ylabel= r'$\log L_{\mathrm{DRW}}\ \mathrm{in}\ '+options.band[0]+r'$'
        xlabel= r'$\log L_{\mathrm{PL}}\ \mathrm{in}\ '+options.band[0]+r'$'
    elif options.plottype == 'loglike3z':
        qsos= open_qsos()
        qsoDict= {}
        ii=0
        for qso in qsos:
            qsoDict[qso.oname.strip().replace(' ', '')+'.fit']= ii
            ii+= 1
        ys= []
        xs= []
        for key in params[0].keys():
            if not params[1].has_key(key): continue
            ys.append(params[1][key]['loglike']
                      -params[0][key]['loglike'])
            xs.append(qsos[qsoDict[key]].z)
        ys= nu.array(ys)
        xs= nu.array(xs)
        yrange=[-15.,15.]
        xrange=[0.,5.5]
        ylabel= r'$\log L_{\mathrm{DRW}}\ \mathrm{in}\ '+options.band[0]+r' - \log L_{\mathrm{PL}}\ \mathrm{in}\ '+options.band[0]+r'$'
        xlabel= r'$z\ (\mathrm{redshift})$'
    elif options.plottype == 'Acz':
        qsos= open_qsos()
        qsoDict= {}
        ii=0
        for qso in qsos:
            qsoDict[qso.oname.strip().replace(' ', '')+'.fit']= ii
            ii+= 1
        ys= []
        xs= []
        for key in params[0].keys():
            ys.append(params[0][key]['logAgr'])
            xs.append(qsos[qsoDict[key]].z)
        ys= nu.array(ys)/2.
        xs= nu.array(xs)
        yrange=[-9.21/2.,0.]
        xrange=[0.,5.5]
        ylabel= r'$\log A^c_{'+options.band+r'}\ \mathrm{(amplitude\ at\ 1\ yr)}$'
        xlabel= r'$z\ (\mathrm{redshift})$'
    elif options.plottype == 'gcz':
        qsos= open_qsos()
        qsoDict= {}
        ii=0
        for qso in qsos:
            qsoDict[qso.oname.strip().replace(' ', '')+'.fit']= ii
            ii+= 1
        ys= []
        xs= []
        for key in params[0].keys():
            ys.append(params[0][key]['gammagr'])
            xs.append(qsos[qsoDict[key]].z)
        ys= nu.array(ys)
        xs= nu.array(xs)
        yrange=[0.,1.25]
        xrange=[0.,5.5]
        ylabel= r'$\gamma^c_{'+options.band+r'}\ \mathrm{(power\ law\ exponent)}$'
        xlabel= r'$z\ (\mathrm{redshift})$'
    elif options.plottype == 'AsAc': #files: g, Ac, KS11
        xs= []
        ys= []
        if not options.imax is None:
            qsos= open_qsos()
            qsos= qsos[(qsos.mags[:,3] < options.imax)]
            qsoDict= {}
            ii=0
            for qso in qsos:
                qsoDict[qso.oname.strip().replace(' ', '')+'.fit']= ii
                ii+= 1               
        for key in params[0].keys():
            if not key in qsoDict.keys(): continue
            try:
                if options.type == 'powerlawSFratios':
                    Ac= params[1][key]['logAri']/2.+0.14
                    A= params[0][key]['logA']/2.
                else:
                    Ac= params[1][key]['logAgr']/2.
                    A= params[0][key]['logA']/2.                   
                s= params[2][key]['s']
                ys.append(Ac)
                xs.append(A+nu.log(1.-s))
            except KeyError:
                continue
        xs= nu.array(xs).reshape(len(xs))
        ys= nu.array(ys).reshape(len(ys))
        xrange=[-9.21/2.,0.]
        yrange=[-9.21/2.,0.]
        xlabel= r'$\log [ (1-s ) A_{'+options.band+r'}]\ \mathrm{(amplitude\ at\ 1\ yr)}$'
        ylabel= r'$\log A^c_{'+options.band+r'}\ \mathrm{(amplitude\ at\ 1\ yr)}$'
    elif options.plottype == 'st': #DRW
        ys= nu.array([p['logl'] for p in params[0].values()]).reshape(len(params[0]))
        ys/= nu.log(10.)
        ys+= nu.log10(365.)
        xs= nu.array([p['loga2'] for p in params[0].values()]).reshape(len(params[0]))
        xs/= nu.log(10.)
        xs+= nu.log10(2.)-ys
        #print nu.amin(ys), nu.amax(ys)
        xrange=[-6.,0.]
        yrange=[-.5,3.5]
        xlabel= r'$\log_{10} \sigma^2_{'+options.band+r'}\ \mathrm{(short-timescale\ variance)}$'
        ylabel= r'$\log_{10} \tau_{'+options.band+r'}\ \mathrm{(damping\ timescale\ in\ days)}$'
    elif options.plottype == 'a2l': #DRW
        ys= nu.array([p['logl'] for p in params[0].values()]).reshape(len(params[0]))
        xs= nu.array([p['loga2'] for p in params[0].values()]).reshape(len(params[0]))
        xrange=[-2.5,0.]
        yrange=[-10.,4.]
        xlabel= r'$\log a^2_{'+options.band+r'}\ \mathrm{(variance)}$'
        ylabel= r'$\log \tau_{'+options.band+r'}\ \mathrm{(damping\ timescale\ in\ yr)}$'
    elif options.plottype == 'Al': #DRW
        ys= nu.array([p['logl'] for p in params[0].values()]).reshape(len(params[0]))
        xs= nu.array([p['loga2'] for p in params[0].values()]).reshape(len(params[0]))
        xs= (nu.log(2.)+xs+nu.log(1.-nu.exp(-1./nu.exp(ys))))/2.
        xrange=[-9.21/2.,0.]
        yrange=[-10.,4.]
        xlabel= r'$\log A_{'+options.band+r'}\ \mathrm{(amplitude\ at\ 1\ yr)}$'
        ylabel= r'$\log \tau_{'+options.band+r'}\ \mathrm{(damping\ timescale\ in\ yr)}$'
    bovy_plot.bovy_print()
    bovy_plot.scatterplot(xs,ys,'k,',onedhists=True,
                          yrange=yrange,
                          xrange=xrange,bins=30,
                          xlabel=xlabel,
                          ylabel=ylabel)
    if options.plottype == 'AA' or options.plottype == 'gg' \
           or options.plottype == 'loglike2' \
           or options.plottype == 'AsAc' \
           or options.plottype == 'ss' \
           or options.plottype == 'AcAc':
        bovy_plot.bovy_plot(nu.array(xrange),nu.array(xrange),'0.5',
                            overplot=True)
    #also fit if desired
    if options.linearfit:
        indx= (xs > xrange[0])*(xs < xrange[1])\
            *(ys > yrange[0])*(ys < yrange[1])
        b= _fit_linear(xs[indx],ys[indx],removeoutliers=options.removeoutliers)
        bovy_plot.bovy_plot(nu.array(xrange),nu.array(xrange)+b,'0.25',
                            overplot=True)
        bovy_plot.bovy_text(ylabelnounits+r'$ = $'+
                            xlabelnounits+r'$ %4.2f$' % b,
                            bottom_right=True,color='0.25')
    bovy_plot.bovy_end_print(options.plotfilename)
    return None

def _fit_linear(xs,ys,removeoutliers=False):
    """linear fit to xs and ys"""
    ndata= len(xs)
    if removeoutliers:
        sigma= nu.std(ys-xs)
        indx= (nu.fabs(ys-xs) < 3*sigma)
        return nu.sum(ys[indx]-xs[indx])/ndata      
    else:
        return nu.sum(ys-xs)/ndata

def open_qsos(s82file='../data/S82qsos.fits'):
    return _load_fits(s82file)

def get_options():
    usage = "usage: %prog [options] <savefilename>\n\nsavefilename= name of the file that holds the fits"
    parser = OptionParser(usage=usage)
    parser.add_option("-b","--band",dest='band',default='r',
                      help="band(s) to fit")
    parser.add_option("-t","--type",dest='type',default='powerlawSF',
                      help="Type of model to fit (powerlawSF or DRW)")
    parser.add_option("-o",dest='plotfilename',default=None,
                      help="Name for the file that will hold the plot")
    parser.add_option("--plottype",dest='plottype',default='Ag',
                      help="Type of plot")
    parser.add_option("--imax",dest='imax',default=None,type='float',
                      help="Maximum i-band magnitude")
    parser.add_option("--kmeansOut",action="store_true", dest="kmeansOut",
                      default=False,
                      help="Input is actually the output of the K-means fitting routine")
    parser.add_option("--linearfit",action="store_true", dest="linearfit",
                      default=False,
                      help="Also perform a linear fit and plot it")
    parser.add_option("--removeoutliers",action="store_true",
                      dest="removeoutliers",
                      default=False,
                      help="Remove outliers when doing the linear fit")
    parser.add_option("--band2",dest='band2',default=None,
                      help="2nd band(s) to fit")

    return parser

if __name__ == '__main__':
    plotFits(get_options())
