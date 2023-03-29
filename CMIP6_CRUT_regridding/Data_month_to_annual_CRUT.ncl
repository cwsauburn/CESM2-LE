; Read in observation, month to annual, two time periods are 1901-2014, 1951-2014
  
;  ye_be = 1901
; ******** Add Observation, it is global gridded data **************************

; Add cru_ts4.03.1901.2018.pre.dat.nc 0.5*0.5 regridded data
;; short tempanomaly time 1416 x lat 360 x lon 720
;; Time period(190101-201812)

  CRUT_f_ID = addfile("cru_ts4.03.1901.2018.pre.dat.nc","r")
  pre       = CRUT_f_ID->pre
  pre_size  = dimsizes(pre)

; lat 72*lon 144  
  mask_f_ID = addfile("mask_2_5_degree.nc","r")
  mask_lat  = mask_f_ID->lat
  mask_lon  = mask_f_ID->lon
  lsm       = mask_f_ID->lsm
  dimlsm    = dimsizes(lsm)
  posilat   = new(dimlsm(0),"integer")
  posilon   = new(dimlsm(1),"integer")
    
  regrid_pre      = new((/pre_size(0),dimlsm(0),dimlsm(1)/),typeof(pre),pre@_FillValue) 
  regrid_pre!0    = "time"
  regrid_pre&time = pre&time
  regrid_pre!1    = "lat"
  regrid_pre&lat  = lsm&lat 
  regrid_pre!2    = "lon"
  regrid_pre&lon  = lsm&lon  
  
  ; negative latitude
  do lati = 0,35
     if ((mask_lat(lati)-ceil(mask_lat(lati))).gt.0.5) then
         posilat(lati) = floattointeger((ceil(mask_lat(lati))+89)*2)
     else     
         posilat(lati) = doubletoint((ceil(mask_lat(lati))+89)*2+1)
     end if
  end do  
   
  ; positive latitude   
  do lati = 36,dimlsm(0)-1
     if ((mask_lat(lati)-floor(mask_lat(lati))).gt.0.5) then
         posilat(lati) = doubletoint(359-abs(floor(mask_lat(lati))-89)*2)
     else     
         posilat(lati) = doubletoint(359-abs(floor(mask_lat(lati))-89)*2-1)
     end if  
  end do
    
  ; positive longitude
  do lonj = 0, 71
     if ((mask_lon(lonj)-floor(mask_lon(lonj))).gt.0.5) then
        posilon(lonj) = doubletoint((floor(mask_lon(lonj))-0)*2+361)
     else
        posilon(lonj) = doubletoint((floor(mask_lon(lonj))-0)*2+360)
     end if
  end do
        
  ; positive longitude      
  do lonj = 72, dimlsm(1)-1
     if (((360-mask_lon(lonj))-floor(360-mask_lon(lonj))).gt.0.5) then
        posilon(lonj) = doubletoint((ceil(mask_lon(lonj)-360)+179)*2)
     else
        posilon(lonj) = doubletoint((ceil(mask_lon(lonj)-360)+179)*2+1)
     end if
  end do
   
  do lati = 0,dimlsm(0)-1
     do lonj = 0, dimlsm(1)-1
        regrid_pre(:,lati,lonj) = pre(:,posilat(lati),posilon(lonj))
     end do  
  end do
   
  CRUT_ann_all         = month_to_annual(regrid_pre((ye_be-1901)*12:(2014-1901+1)*12-1,:,:),1)
  CRUT_ann_all_reorder = CRUT_ann_all(lon|:,lat|:,year|:)
  
  CRUT_ID = addfile("Ann_prcp_CRUT_1901_2014.nc","c")
  CRUT_ID->CRUT = CRUT_ann_all_reorder   
