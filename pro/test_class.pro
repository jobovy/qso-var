FUNCTION TEST_CLASS, class, truth, prior=prior, verbose=verbose, rrlyraeprior=rrlyraeprior, threshold=threshold
IF ~keyword_set(prior) THEN prior=0.05
IF ~keyword_set(rrlyraeprior) THEN rrlyraeprior= prior
IF ~keyword_set(threshold) THEN threshold= 0.5
truth= truth[where(truth.i GE 17.8 and (truth.ra GE 315. or truth.ra LE 60.))]
indx= lonarr(n_elements(class.key))-1L
;;prepare output
newstr= {qsologlike:0D, starloglike:0D,rrlyraeloglike:0D,logpqso:0D,key:''}
newstr= replicate(newstr,n_elements(truth.z))
out= struct_combine(truth,newstr)
for ii=0L, n_elements(class.key)-1 do begin
    logpqso= class[ii].qsologlike+alog(prior)-logsum([class[ii].qsologlike+alog(prior),class[ii].starloglike+alog(1.-prior-rrlyraeprior),class[ii].rrlyraeloglike+alog(rrlyraeprior)])
                                ;if class[ii].qsologlike GE
                                ;(class[ii].starloglike+lnprior) then
                                ;begin
    if logpqso GT alog(threshold) then begin
        if keyword_set(verbose) then print, pqso, where(strmatch(truth.oname,class[ii].key))
        indx[ii]= where(strmatch(truth.oname,class[ii].key))
    endif
    matchindx= where(strmatch(truth.oname,class[ii].key),cnt)
    if cnt ne 0 then out[matchindx].qsologlike= class[ii].qsologlike
    if cnt ne 0 then out[matchindx].starloglike= class[ii].starloglike
    if cnt ne 0 then out[matchindx].rrlyraeloglike= class[ii].rrlyraeloglike
    if cnt ne 0 then out[matchindx].key= class[ii].key
    if cnt ne 0 then out[matchindx].logpqso= logpqso
endfor
sample= truth[indx[where(indx GE 0.)]]
print, strtrim(string(n_elements(where(sample.z GE 0.1))),2)+" / "+$
  strtrim(string(n_elements(sample.z)),2)
print, double(n_elements(where(sample.z GE 0.1)))/n_elements(sample.z), $
  double(n_elements(where(sample.z GE 0.1)))/n_elements(where(truth.z GE 0.1))
;print, n_elements(where(truth.z GE 0.1))
RETURN, out
END
