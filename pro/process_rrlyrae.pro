PRO PROCESS_RRLYRAE
datadir='../data/table1_sesar2010/'
outdir='../data/rrlyrae/'
readfmt= 'F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F'
outStruct= {rrlyrae, ra:0D, dec:0D, $
            mjd_u:0D, u:0D, err_u:0D, $
            mjd_g:0D, g:0D, err_g:0D, $
            mjd_r:0D, r:0D, err_r:0D, $
            mjd_i:0D, i:0D, err_i:0D, $
            mjd_z:0D, z:0D, err_z:0D,ID:0L}
;;Get filenames
files= file_search(datadir+'*.dat',/test_regular)
nfiles= n_elements(files)
FOR ii=0L, nfiles-1 DO BEGIN
    ;;Read the file
    print, format = '("Working on ",i7," of ",i7,a1,$)', $
          ii+1,nfiles,string(13B)
    readcol, files[ii], ra, dec, mjd_u, u, err_u, $
      mjd_g, g, err_g, $
      mjd_r, r, err_r, $
      mjd_i, i, err_i, $
      mjd_z, z, err_z, format=readfmt, /silent
    badra= where(ra LT 0., cnt)
    if cnt ne 0 then ra[badra]+= 360.
    thisout= replicate(outStruct,n_elements(u))
    basefile= file_basename(files[ii])
    result= strsplit(basefile,".",/extract)
    print, result
    thisout.id= result[0]
    thisout.ra= ra
    thisout.dec= dec
    thisout.mjd_u= mjd_u
    thisout.u= u
    thisout.err_u= err_u
    thisout.mjd_g= mjd_g
    thisout.g= g
    thisout.err_g= err_g
    thisout.mjd_r= mjd_r
    thisout.r= r
    thisout.err_r= err_r
    thisout.mjd_i= mjd_i
    thisout.i= i
    thisout.err_i= err_i
    thisout.mjd_z= mjd_z
    thisout.z= z
    thisout.err_z= err_z
    mwrfits, thisout, outdir+strcompress(hogg_iau_name(ra[0],dec[0]),/remove_all)+'.fit',/create
ENDFOR
END
