; Organize CESM2-LE data 

 beyear = 1850
 edyear = 2014
 term   = edyear-beyear+1
 
 ;------ temperature -------------------
 ;------ Create mask for 85N-60S, and cover continent --------------
 flsm     = addfile("mask_2_5_degree_rg2CESM2.nc", "r") ;  [lat | 192] x [lon | 288]
 lsm_temp = flsm->lsm
 lsm      = toint(lsm_temp)
 lsm_size = dimsizes(lsm)
 lat_2_5  = lsm_temp&lat
 lon_2_5  = lsm_temp&lon

 mask_2_5_degree          = mask(lsm,((lsm.eq.1).or.(lsm.eq.2)),True)
 mask_2_5_degree(0:31,:)    = -127   ;60S
 mask_2_5_degree(186:191,:) = -127   ;85N
 mask_2_5_degree          = where(mask_2_5_degree.eq.2,1,mask_2_5_degree)
 mask_2_5_degree!0        = "lat"
 mask_2_5_degree!1        = "lon"

 mask_2_5_degree = mask_2_5_degree(::-1,:)

; ******define constant value, read in prepared obs data *********************
 beindex = (beyear-1850)*12
 edindex = (edyear-1850+1)*12-1

;TREFHT(time 120, lat 192, lon 288)
 diri      = "/glade/campaign/cgd/cesm/CESM2-LE/timeseries/atm/proc/tseries/month_1/TREFHT/"
 slide     = "199001-199912"
 slidename = stringtochar(systemfunc("csh -c 'cd " + diri + " ; ls *"+slide+"*.nc'"))
 ;print(slidename(0,28:36))

 ensName   = slidename(:,28:36)
 szName    = dimsizes(ensName) ; get ensemble members' name

 sampID = addfile(diri+"b.e21.BHISTcmip6.f09_g17.LE2-1301.010.cam.h0.TREFHT.199001-199912.nc","r")
 samp   = sampID->TREFHT
 szsamp = dimsizes(samp) ; [time | 120] x [lat | 192] x [lon | 288]

 tas_ann_nE = new((/term, szsamp(1),szsamp(2),szName(0)/),typeof(samp),samp@_FillValue)

 do i = 0, szName(0)-1
     ensFile = systemfunc("csh -c 'cd " + diri + " ; ls *"+ensName(i,:)+"*.nc'")
     fileID  = addfiles(diri+ensFile, "r")
     ListSetType(fileID, "cat")

     TREFHT_temp  = fileID[:]->TREFHT
     tas_ann_nE(:,:,:,i)  = month_to_annual(TREFHT_temp(beindex:edindex,:,:),1)-273.15
	 ;tas_ann_nE(:,:,:,i) = runave_n(TREFHT(:,:,:),10,0,0)
	 
     delete(ensFile)
     delete(fileID)
     delete(TREFHT_temp)
     ;delete(TREFHT)
 end do

 system("/bin/rm -f Ann_tas_CESM2_LE_global_1850_2014_trend.nc")
 mod_ID = addfile("Ann_tas_CESM2_LE_global_1850_2014_trend.nc","c")

 mod_ID ->tas_ann_nE = tas_ann_nE


 mod_ID          = addfile("Ann_tas_CESM2_LE_global_1901_2014_trend.nc","r")
 tas_ann_nE_temp = mod_ID->tas_ann_nE 

 tas_ann_nE      = tas_ann_nE_temp(50:113,:,:,:)
 
 system("/bin/rm -f Ann_tas_CESM2_LE_global_1951_2014_trend.nc")
 mod_ID = addfile("Ann_tas_CESM2_LE_global_1951_2014_trend.nc","c")

 mod_ID ->tas_ann_nE = tas_ann_nE

 ;-------- precipitation ---------------------
 beyear = 1951
 edyear = 2014
 term   = edyear-beyear+1

 beindex = (beyear-1850)*12
 edindex = (edyear-1850+1)*12-1

 ucf = 60.0*60.0*24.0*365.0*(1.0/12.0)

 ; PRECT(time, lat 192, lon 288)
 diri      = "/glade/campaign/cgd/cesm/CESM2-LE/timeseries/atm/proc/tseries/month_1/PRECT/"
 slide     = "199001-199912"
 slidename = stringtochar(systemfunc("csh -c 'cd " + diri + " ; ls *"+slide+"*.nc'"))
 ;print(slidename(0,28:36))

 ensName   = slidename(:,28:36)
 szName    = dimsizes(ensName) ; get ensemble members' name

 sampID = addfile(diri+"b.e21.BHISTcmip6.f09_g17.LE2-1301.010.cam.h0.PRECT.199001-199912.nc","r")
 samp   = sampID->PRECT
 szsamp = dimsizes(samp) ; [time | 120] x [lat | 192] x [lon | 288]

 prcp_ann_nE = new((/term, szsamp(1),szsamp(2),szName(0)/),typeof(samp),samp@_FillValue)

 do i = 0, szName(0)-1
     ensFile = systemfunc("csh -c 'cd " + diri + " ; ls *"+ensName(i,:)+"*.nc'")
     fileID  = addfiles(diri+ensFile, "r")
     ListSetType(fileID, "cat")

     PRECT_temp  = fileID[:]->PRECT
     PRECT       = month_to_annual(PRECT_temp(beindex:edindex,:,:),1)
     prcp_ann_nE(:,:,:,i) = runave_n(PRECT(:,:,:),5,0,0)
	 
     delete(ensFile)
     delete(fileID)
     delete(PRECT_temp)
     delete(PRECT)
 end do
 prcp_ann_nE = prcp_ann_nE*ucf

 system("/bin/rm -f Ann_prcp_CESM2_LE_global_1951_2014_trend.nc")
 mod_ID = addfile("Ann_prcp_CESM2_LE_global_1951_2014_trend.nc","c")
 mod_ID ->prcp_ann_nE = prcp_ann_nE







