PRO STAR_LC, n=n, seed=seed
IF ~keyword_set(seed) THEN seed= 1L
IF ~keyword_set(n) THEN n= 10000L
;;Read data
in= read_standardstar()
lc= mrdfits('../data/star_lc.fit',1)
;;subsample
indx= where(in.dec LT 1.25 AND in.dec GT -1.25 AND (in.ra LT 59. OR in.ra GT 310) AND $
            (in.mmed[1]+in.extinction[1]-in.mmed[2]-in.extinction[2]) GT 0.2 AND $
            (in.mmed[1]+in.extinction[1]-in.mmed[2]-in.extinction[2]) LT 0.48 AND $
            (in.mmed[1]+in.extinction[1]) GT 14. AND $
            (in.mmed[1]+in.extinction[1]) LT 20.2)
in= in[indx]
x= lindgen(n_elements(in.ra))
y= randomu(seed,n_elements(in.ra))
z= x[sort(y)]
z= z[0:n-1]
in= in[z]
;;Match for each datapoint
outDir= '../data/star/'
for ii=0L, n_elements(in.ra)-1 do begin
    print, "working on "+strtrim(string(ii),2)
    spherematch, in[ii].ra, in[ii].dec, lc.ra, lc.dec, 0.5/3600.,$
      iindx, lindx, maxmatch=0
    if lindx[0] eq -1 then continue
    print, n_elements(lindx)
    out= lc[lindx]
    mwrfits,out,outDir+strcompress(hogg_iau_name(in[ii].ra,in[ii].dec),/remove_all)+'.fit',/create
endfor
end
