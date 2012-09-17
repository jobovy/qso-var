PRO NUVX_LC
;;Read data
in= mrdfits('../data/nUVX.fit',1)
lc= mrdfits('../data/nuvx_lc.fit',1)
;;Match for each datapoint
outDir= '../data/nuvx/'
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
