FUNCTION READ_STANDARDSTAR, filename=filename
IF ~keyword_set(filename) THEN filename='../data/stripe82calibStars_v2.6.dat'

FMT= 'A,F,F,F,F,I,F,I,F,F,F,F,F,I,F,F,F,F,F,I,F,F,F,F,F,I,F,F,F,F,F,I,F,F,F,F,F'
READCOL,filename,F=FMT,$
  zevel,ra,dec,rarms,decrms,ntot,ar,$
  nobs_u, mmed_u, mmu_u, msig_u, mrms_u, chi2_u, $
  nobs_g, mmed_g, mmu_g, msig_g, mrms_g, chi2_g, $
  nobs_r, mmed_r, mmu_r, msig_r, mrms_r, chi2_r, $
  nobs_i, mmed_i, mmu_i, msig_i, mrms_i, chi2_i, $
  nobs_z, mmed_z, mmu_z, msig_z, mrms_z, chi2_z

outStruct= {ra:0D, dec:0D,rarms:0D,decrms:0D,ntot:0L,ar:0D, $
            nobs:lonarr(5), mmed:dblarr(5),mmu:dblarr(5),msig:dblarr(5),$
            mrms:dblarr(5),chi2:dblarr(5),extinction:dblarr(5)}
out= replicate(outStruct,n_elements(ra))
out.ra= ra
out.dec= dec
out.rarms= rarms
out.decrms= decrms
out.ntot= ntot
out.ar= ar
out.extinction[2]= out.ar
out.extinction[0]= 1.873*out.ar
out.extinction[1]= 1.377*out.ar
out.extinction[3]= 0.758*out.ar
out.extinction[4]= 0.537*out.ar
out.nobs[0]= nobs_u
out.nobs[1]= nobs_g
out.nobs[2]= nobs_r
out.nobs[3]= nobs_i
out.nobs[4]= nobs_z
out.mmed[0]= mmed_u
out.mmed[1]= mmed_g
out.mmed[2]= mmed_r
out.mmed[3]= mmed_i
out.mmed[4]= mmed_z
out.mmu[0]= mmu_u
out.mmu[1]= mmu_g
out.mmu[2]= mmu_r
out.mmu[3]= mmu_i
out.mmu[4]= mmu_z
out.msig[0]= msig_u
out.msig[1]= msig_g
out.msig[2]= msig_r
out.msig[3]= msig_i
out.msig[4]= msig_z
out.mrms[0]= mrms_u
out.mrms[1]= mrms_g
out.mrms[2]= mrms_r
out.mrms[3]= mrms_i
out.mrms[4]= mrms_z
out.chi2[0]= chi2_u
out.chi2[1]= chi2_g
out.chi2[2]= chi2_r
out.chi2[3]= chi2_i
out.chi2[4]= chi2_z
return, out
end
