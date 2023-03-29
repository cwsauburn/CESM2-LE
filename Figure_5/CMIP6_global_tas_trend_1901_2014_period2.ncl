; calculate global continent weighted average evapotranspiration

; ***** read in CMIP6 temperature data ***********************
 ; size is lon 144,lat 72, year 114, ensemble 227
 ann_tas1 = addfile("~/Ann_tas_CMIP6_1_1901_2014.nc","r")
 ann_tas2 = addfile("~/Ann_tas_CMIP6_2_1901_2014.nc","r")
 ann_tas3 = addfile("~/Ann_tas_CMIP6_3_1901_2014.nc","r")

 tas_ann_nE1 = ann_tas1->CMIP6
 tas_ann_nE2_temp = ann_tas2->CMIP6
 tas_ann_nE3 = ann_tas3->CMIP6

 TANT_size   = dimsizes(tas_ann_nE2_temp)
 tas_ann_nE2 = new((/TANT_size(0),TANT_size(1),TANT_size(2),TANT_size(3)-1/),typeof(tas_ann_nE2_temp),tas_ann_nE2_temp@_FillValue)
 tas_ann_nE2(:,:,:,0:6) = tas_ann_nE2_temp(:,:,:,0:6)
 tas_ann_nE2(:,:,:,7:TANT_size(3)-2) = tas_ann_nE2_temp(:,:,:,8:TANT_size(3)-1)

 tas_ann_nE_temp = new((/TANT_size(0),TANT_size(1),2014-1901+1,220/),typeof(tas_ann_nE1),tas_ann_nE1@_FillValue)

 tas_ann_nE_temp(:,:,:,0:131)   = tas_ann_nE1
 tas_ann_nE_temp(:,:,:,132:182) = tas_ann_nE2
 tas_ann_nE_temp(:,:,:,183:219) = tas_ann_nE3

 tas_ann_nE_temp!0 = "lon"
 tas_ann_nE_temp!1 = "lat"
 tas_ann_nE_temp!2 = "year"
 tas_ann_nE_temp!3 = "ensemble"

 tas_ann_nE_temp2  = tas_ann_nE_temp(year|:,lat|:,lon|:,ensemble|:)
 
 tas_ann_nE        = tas_ann_nE_temp2(59:89,:,:,:)
 
 szsamp = dimsizes(tas_ann_nE)
 TSA    = new((/2,szsamp(1),szsamp(2),szsamp(3)/),typeof(tas_ann_nE),tas_ann_nE@_FillValue)
 
 do ensj = 0, 219
	TSA(:,:,:,ensj) = trend_manken(rm_single_dims(tas_ann_nE(:,:,:,ensj)), False, 0)
 end do

; ****** write it into nc file *****************************	   
 system("/bin/rm -f TSA_CMIP6_global_1901_2014_period2.nc")
 RPC_ID = addfile("TSA_CMIP6_global_1901_2014_period2.nc","c")

 RPC_ID ->TSA = TSA	   




	   
