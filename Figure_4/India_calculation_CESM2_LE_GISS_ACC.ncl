; calculate India continent weighted average temperature
; India defination {9.6:25.3},{72.8:90}  106:123, 58:71 

; ******** Read in mask file, it is global gridded data ***************************
; size is lat 192 * lon 288, lsm=0, ocean
 maskID = addfile("mask_2_5_degree_rg2CESM2.nc","r")
 lsm    = maskID->lsm
 lsm    = lsm(::-1,:)
 
; ******** Read in weight, it is global gridded data ******************************
; size is lat 192 * lon 288
 WID        = addfile("~/CLM4.0CN.h0.fixed.nc","r")
 carea_temp = WID->area
 carea_temp@_FillValue = 1e+36

; ****** read in CESM2-LE, calculate anomaly and mask by ocean ********************
 CESMID = addfile("Ann_tas_CESM2_LE_global_1850_2014_trend.nc","r")
 tas_ann_nE_temp0 = CESMID->tas_ann_nE
 tas_ann_nE_temp  = tas_ann_nE_temp0(1901-1850:2014-1850,:,:,:)
 
 tas_ann_nE_temp!0 = "time"
 tas_ann_nE_temp!1 = "lat"
 tas_ann_nE_temp!2 = "lon"
 tas_ann_nE_temp!3 = "ens" 

 tas_ann_nE       = tas_ann_nE_temp(time|:,ens|:,lat|:,lon|:)

 tas_ann_nE_cli   = rm_single_dims(dim_avg_n(tas_ann_nE(1951-1901:1980-1901,:,:,:),0))
 
 CESM2_T_temp     = tas_ann_nE(:,:,106:123,58:71)
 tas_ann_nE_cli_T = tas_ann_nE_cli(:,106:123,58:71)
  
 do yeari = 0, 113
    do ensj = 0, 89
	   CESM2_T_temp(yeari,ensj,:,:) = where(lsm(106:123,58:71).eq.1,rm_single_dims(CESM2_T_temp(yeari,ensj,:,:)),CESM2_T_temp@_FillValue)
	   CESM2_T_temp(yeari,ensj,:,:) = rm_single_dims(CESM2_T_temp(yeari,ensj,:,:))-rm_single_dims(tas_ann_nE_cli_T(ensj,:,:))
    end do
 end do
 
 CESM2_T_temp = runave_n(CESM2_T_temp,5,0,0)
 
; ******** Add Observation, it is global gridded data ******************************
 ; lon 288, lat 192, year 114
 ; base: 1951-1980
 GISS_f_ID     = addfile("~/Ann_tas_GISS_1901_2014_1_degree.nc","r")
 GISS_ann_temp = GISS_f_ID->tas
 
 GISS_ann = runave_n(GISS_ann_temp(:, :, :),5,0,2)

 obs_T_temp = GISS_ann(58:71,106:123,:)

; ******* Calculate area average for CESM2_LE and observation ********************** 
 careaT = carea_temp(106:123,58:71)
 careaT = where(lsm(106:123,58:71).eq.1,careaT,careaT@_FillValue)
 sumT = sum(careaT)

 ; time 114, lat , lon , ens 90
 ST = dimsizes(CESM2_T_temp)
  
 CESM2_T_temp2 = CESM2_T_temp
 CESM2_T_temp2(:,:,:,:) = CESM2_T_temp@_FillValue
 
 obs_T_temp2 = obs_T_temp
 obs_T_temp2(:,:,:) = CESM2_T_temp@_FillValue

 do loni = 0, ST(3)-1
    do latj = 0, ST(2)-1
	   CESM2_T_temp2(:,:,latj,loni) = CESM2_T_temp(:,:,latj,loni)*careaT(latj,loni)/sumT
	   obs_T_temp2(loni,latj,:)     = obs_T_temp(loni,latj,:)*careaT(latj,loni)/sumT
	end do
 end do  
	   
 CESM2_T = rm_single_dims(dim_sum_n(dim_sum_n(CESM2_T_temp2, 3), 2)) 
 obs_T   = rm_single_dims(dim_sum_n(dim_sum_n(obs_T_temp2,1),0))	 

 CESM2_T!0 = "year_1901_2014"
 CESM2_T!1 = "ensemble"
 
 obs_T!0 = "year_1901_2014"

; ************* calculate RPS *******************************************
; calculate RPS following the method described in Smith et al 2020
; find the top 10 ensemble members which have the similar IPO with CESM2-LE ensemble mean
 
 RPS = variance(obs_T)/variance(dim_avg_n(CESM2_T,1))
 
 adjust_AMV = dim_avg_n(CESM2_T,1)*RPS
 
 RMSE_ens = new((/90/),"float")
 
 do ensi = 0, 89
    ;RMSE_ens(ensi) = sqrt(avg((CESM2_T(:, ensi)-adjust_AMV)^2))
    RMSE_ens(ensi) = escorc(CESM2_T(:, ensi),adjust_AMV)
 end do
  
 CESM2_AMV_ord = dim_pqsort(RMSE_ens,1)
 
 ;print(CESM2_IPO_ord)
 
; ****** write it into nc file *****************************	   
 system("/bin/rm -f Ann_India_CESM2_LE_GISS_tas_1901_2014.nc")
 RPC_ID = addfile("Ann_India_CESM2_LE_GISS_tas_1901_2014.nc","c")

 RPC_ID ->CESM2_AMV = CESM2_T	   
 RPC_ID ->GISS_AMV  = obs_T

 RPC_ID ->CESM2_AMV_ord = CESM2_AMV_ord


	   
