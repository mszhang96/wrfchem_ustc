begin

 fi=addfile("surfdata.nc","r")
 pft0=fi->PCT_PFT(:,:,:)

 fo=addfile("wrfinput_d01","w")
 pct_pft=fo->PCT_PFT(:,:,:,:)
 pct_pft(0,:,:,:) = (/doubletofloat(pft0(:,:,:))/)

 fo->PCT_PFT=pct_pft
end

