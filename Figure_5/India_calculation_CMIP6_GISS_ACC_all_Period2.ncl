; calculate global continent weighted average temperature 
; India defination {9.6:25.3},{72.8:90} lat lon
; 39:45,29:35

; ******** Read in mask file, it is global gridded data ***************************
; size is lat 72 * lon 144
 flsm = addfile("./mask_2_5_degree.nc", "r") ;  [lat | 72] x [lon | 144]
 lsm  = flsm->lsm  
 
; ******** Read in weight, it is global gridded data ******************************
; size is lat 72 * lon 144
 WID        = addfile("./area_2_5_degree.nc","r")
 carea_temp = doubletofloat(WID->CAREA)
 
; ****** read in CMIP6, calculate anomaly and mask by ocean ********************
 ; read in all data together, lon, lat, year, ensemble 144*72*114*(132+52+37=221)
 ann_tas1 = addfile("./Ann_tas_CMIP6_1_1901_2014.nc","r")
 ann_tas2 = addfile("./Ann_tas_CMIP6_2_1901_2014.nc","r")
 ann_tas3 = addfile("./Ann_tas_CMIP6_3_1901_2014.nc","r")
  
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
 
 tas_ann_nE        = tas_ann_nE_temp(year|:,ensemble|:,lat|:,lon|:)

 tas_ann_nE_cli    = rm_single_dims(dim_avg_n(tas_ann_nE(1951-1901:1980-1901,:,:,:),0))

 CMIP6_T_temp      = tas_ann_nE(:,:,39:45,29:35)
 tas_ann_nE_cli_T  = tas_ann_nE_cli(:,39:45,29:35)
 
 do yeari = 0, 113
    do ensj = 0, 219
	   CMIP6_T_temp(yeari,ensj,:,:) = where(lsm(39:45,29:35).eq.1,rm_single_dims(CMIP6_T_temp(yeari,ensj,:,:)),CMIP6_T_temp@_FillValue)
	   CMIP6_T_temp(yeari,ensj,:,:) = rm_single_dims(CMIP6_T_temp(yeari,ensj,:,:))-rm_single_dims(tas_ann_nE_cli_T(ensj,:,:))
	end do
 end do	
 
 CMIP6_T_temp = runave_n(CMIP6_T_temp,5,0,0)
 
; ******** Add Observation, it is global gridded data ******************************
 GISS_f_ID     = addfile("./GISS_global_2.5.nc","r")
 GISS_ann_temp_temp = GISS_f_ID->GISS

 GISS_ann_temp = month_to_annual(GISS_ann_temp_temp(0:(2014-1880+1)*12-1,:,:),1)
 GISS_ann_all  = GISS_ann_temp(lon|:,lat|:,year|:)
 GISS_ann      = runave_n(GISS_ann_all(:,:,(1901-1880):(2014-1880)),5,0,2)

 obs_T_temp = GISS_ann(29:35,39:45,:)*0.01
 
; ******* Calculate area average for CMIP6 and observation ********************** 
 careaT = carea_temp(39:45,29:35)
 careaT = where(lsm(39:45,29:35).eq.1,careaT,careaT@_FillValue)
 sumT = sum(careaT)

 ; time 114, lat , lon , ens 220
 ST = dimsizes(CMIP6_T_temp)
  
 CMIP6_T_temp2 = CMIP6_T_temp
 CMIP6_T_temp2(:,:,:,:) = CMIP6_T_temp@_FillValue
 
 obs_T_temp2 = obs_T_temp
 obs_T_temp2(:,:,:) = CMIP6_T_temp@_FillValue

 do loni = 0, ST(3)-1
    do latj = 0, ST(2)-1
	   CMIP6_T_temp2(:,:,latj,loni) = CMIP6_T_temp(:,:,latj,loni)*careaT(latj,loni)/sumT
	   obs_T_temp2(loni,latj,:)     = obs_T_temp(loni,latj,:)*careaT(latj,loni)/sumT
	end do
 end do  
	   
 CMIP6_T = rm_single_dims(dim_sum_n(dim_sum_n(CMIP6_T_temp2, 3), 2)) 
 obs_T_temp3   = rm_single_dims(dim_sum_n(dim_sum_n(obs_T_temp2,1),0))	 
 obs_T   = obs_T_temp3(0:59)

 CMIP6_T!0 = "year_1901_2014"
 CMIP6_T!1 = "ensemble"
 
 obs_T!0 = "year_1901_2014"

; *********** add CESM2-LE tas *******************************************
 ; year_1901_2014 = 114 ; ensemble = 90 ;
 cesm2id = addfile("Ann_India_CESM2_LE_GISS_tas_1901_2014.nc","r")
 CESM2_T = cesm2id->CESM2_AMV

 comb_T_temp = new((/114,220+90/),typeof(CMIP6_T))
 comb_T_temp(:,0:219) = CMIP6_T
 comb_T_temp(:,220:309) = CESM2_T
 
 comb_T = comb_T_temp(0:59,:) 

;  ************* calculate RPS *******************************************
; calculate RPS following the method described in Smith et al 2020
; find the top 10 ensemble members which have the similar IPO with CMIP6-LE ensemble mean
 
 RPS = variance(obs_T)/variance(dim_avg_n(comb_T,1))
 
 adjust_AMV = dim_avg_n(comb_T,1)*RPS
 
 RMSE_ens = new((/220+90/),"float")
 
 do ensi = 0, 219+90
    ; RMSE_ens(ensi) = sqrt(avg((CMIP6_T(:, ensi)-adjust_AMV)^2))
    RMSE_ens(ensi) = escorc(comb_T(:, ensi),adjust_AMV)
 end do
  
 CMIP6_AMV_ord = dim_pqsort(RMSE_ens,1)
 
 ;print(CMIP6_IPO_ord)
 
; ****** write it into nc file *****************************	   
 system("/bin/rm -f Ann_India_CMIP6_CESM2_LE_GISS_tas_1901_1960.nc")
 RPC_ID = addfile("Ann_India_CMIP6_CESM2_LE_GISS_tas_1901_1960.nc","c")

 RPC_ID ->CMIP6_AMV = comb_T_temp	   
 RPC_ID ->GISS_AMV  = obs_T_temp3

 RPC_ID ->CMIP6_AMV_ord = CMIP6_AMV_ord


	   
