; calculate global continent weighted average precipitation

; Central eastern US defination {36.9:44.0},{269.7:290.3} lat lon
 lat1 = 50
 lat2 = 53

 lon1 = 108
 lon2 = 115

/;
; UK defination {51.2,58.9},{-6:1.3} lat lon
 lat1 = 56
 lat2 = 59
 
 lon1 = 142
 lon2 = 143

; Western Australia defination {-12.5:-34.6},{113.3:134.5} lat lon
 lat1 = 22
 lat2 = 30

 lon1 = 45
 lon2 = 53
;/

 ;1901-1950
 ;period = 49
 ;period = 59
 ;period = 69
 period = 79


; Inorder to make 5 years moving average, we need to extend the calculating period to 1900-2015
 beginy = 1901
 endy   = 2014
 term   = endy-beginy+1
 ensnum = 227
 
; ******** Read in mask file, it is global gridded data ***************************
; size is lat 72 * lon 144
 flsm = addfile("/hosted/fores/CMIP6/Trend_ans/temperature_data/mask_2_5_degree.nc", "r") ;  [lat | 72] x [lon | 144]
 lsm  = flsm->lsm  
 
; ******** Read in weight, it is global gridded data ******************************
; size is lat 72 * lon 144
 WID        = addfile("/hosted/fores/CMIP6/Trend_ans/temperature_data/RPC/area_2_5_degree.nc","r")
 carea_temp = doubletofloat(WID->CAREA)
 
; ****** read in CMIP6, calculate anomaly and mask by ocean ******************** 
 ; read in all data together, lon, lat, year, ensemble 144*72*114*(228-1 for 1901-2014, 229-1 for 1951-2014)
  ann_prcp1 = addfile("/hosted/fores/CMIP6/Trend_ans/precipitation_data/Ann_prcp_CMIP6_1_1901_2014.nc","r")
  ann_prcp2 = addfile("/hosted/fores/CMIP6/Trend_ans/precipitation_data/Ann_prcp_CMIP6_2_1901_2014.nc","r")
  ann_prcp3 = addfile("/hosted/fores/CMIP6/Trend_ans/precipitation_data/Ann_prcp_CMIP6_3_1901_2014.nc","r")
  
  prcp_ann_nE1 = ann_prcp1->CMIP6
  prcp_ann_nE2_temp = ann_prcp2->CMIP6
  prcp_ann_nE3 = ann_prcp3->CMIP6
  
  TANT_size   = dimsizes(prcp_ann_nE2_temp)

  prcp_ann_nE2 = new((/TANT_size(0),TANT_size(1),TANT_size(2),TANT_size(3)-1/),typeof(prcp_ann_nE2_temp),prcp_ann_nE2_temp@_FillValue)
  prcp_ann_nE2(:,:,:,0:14)                = prcp_ann_nE2_temp(:,:,:,0:14)
  prcp_ann_nE2(:,:,:,15:TANT_size(3)-2)   = prcp_ann_nE2_temp(:,:,:,16:TANT_size(3)-1)
  
  prcp_ann_nE_temp = new((/TANT_size(0),TANT_size(1),term,ensnum/),typeof(prcp_ann_nE1),prcp_ann_nE1@_FillValue)
  
  prcp_ann_nE_temp(:,:,:,0:130)   = prcp_ann_nE1 
  prcp_ann_nE_temp(:,:,:,131:188) = prcp_ann_nE2 
  prcp_ann_nE_temp(:,:,:,189:226) = prcp_ann_nE3

  prcp_ann_nE_temp!0 = "lon"
  prcp_ann_nE_temp!1 = "lat"
  prcp_ann_nE_temp!2 = "year"
  prcp_ann_nE_temp!3 = "ensemble" 
 
  prcp_ann_nE        = prcp_ann_nE_temp(year|:,ensemble|:,lat|:,lon|:)
 
  prcp_ann_nE_cli    = rm_single_dims(dim_avg_n(prcp_ann_nE(1951-1901:1980-1901,:,:,:),0))

  CMIP6_P_temp        = prcp_ann_nE(:,:,lat1:lat2,lon1:lon2)
  prcp_ann_nE_cli_P   = prcp_ann_nE_cli(:,lat1:lat2,lon1:lon2)
 
  do yeari = 0, term-1
     do ensj = 0, 226
	CMIP6_P_temp(yeari,ensj,:,:) = where(lsm(lat1:lat2,lon1:lon2).eq.1,rm_single_dims(CMIP6_P_temp(yeari,ensj,:,:)),CMIP6_P_temp@_FillValue)
        CMIP6_P_temp(yeari,ensj,:,:) = rm_single_dims(CMIP6_P_temp(yeari,ensj,:,:))-rm_single_dims(prcp_ann_nE_cli_P(ensj,:,:))
     end do
  end do

  CMIP6_P_temp = runave_n(CMIP6_P_temp,5,0,0)
 
; ******** Add Observation, it is global gridded data **************************
; precipitation anomaly lon = 144*lat = 72*time = 114
  HadCRUT_f_ID    = addfile("/hosted/fores/CMIP6/Trend_ans/precipitation_data/Ann_prcp_CRUT_1901_2014.nc","r")
  HadCRUT_ann_all = HadCRUT_f_ID->CRUT

  HadCRUT_ann_all_cli = rm_single_dims(dim_avg_n(HadCRUT_ann_all(:,:,1951-1901:1980-1901),2))
  
  do yeari = 0, 2014-1901
      HadCRUT_ann_all(:,:,yeari) = HadCRUT_ann_all(:,:,yeari)-HadCRUT_ann_all_cli
  end do
  
  HadCRUT_ann = runave_n(HadCRUT_ann_all,5,0,2)
  obs_P_temp = HadCRUT_ann(lon1:lon2,lat1:lat2,:)
  
; ******* Calculate area average for CMIP6 and observation ********************** 
  careaT = carea_temp(lat1:lat2,lon1:lon2)
  careaT = where(lsm(lat1:lat2,lon1:lon2).eq.1,careaT,careaT@_FillValue)
  sumT = sum(careaT)

  ; time 114, ens 228, lat , lon 
  ST = dimsizes(CMIP6_P_temp)
  
  CMIP6_P_temp2 = CMIP6_P_temp
  CMIP6_P_temp2(:,:,:,:) = CMIP6_P_temp@_FillValue
 
  obs_P_temp2 = obs_P_temp
  obs_P_temp2(:,:,:) = CMIP6_P_temp@_FillValue

  do loni = 0, ST(3)-1
     do latj = 0, ST(2)-1
       CMIP6_P_temp2(:,:,latj,loni) = CMIP6_P_temp(:,:,latj,loni)*careaT(latj,loni)/sumT
       obs_P_temp2(loni,latj,:)     = obs_P_temp(loni,latj,:)*careaT(latj,loni)/sumT
     end do
  end do  
	   
  CMIP6_P = rm_single_dims(dim_sum_n(dim_sum_n(CMIP6_P_temp2, 3), 2)) 
  obs_P   = rm_single_dims(dim_sum_n(dim_sum_n(obs_P_temp2,1),0))	 

  CMIP6_P!0 = "year_1901_2014"
  CMIP6_P!1 = "ensemble"
 
  obs_P!0 = "year_1901_2014"
  
; *********** add CESM2-LE tas *******************************************
 ; year_1901_2014 = 114 ; ensemble = 90 ;

 cesm2id = addfile("Ann_Central_eastern_US_CESM2_LE_CRUT_prcp_1901_2014.nc","r") 
 ;cesm2id = addfile("Ann_UK_CESM2_LE_CRUT_prcp_1901_2014.nc","r")
 ;cesm2id = addfile("Ann_western_Australia_CESM2_LE_CRUT_prcp_1901_2014.nc","r")

 CESM2_P = cesm2id->CESM2_AMV

 comb_P = new((/114,227+90/),typeof(CMIP6_P))
 comb_P(:,0:226) = CMIP6_P
 comb_P(:,227:316) = CESM2_P

;  ************* calculate RPS *******************************************
; calculate RPS following the method described in Smith et al 2020
; find the top 10 ensemble members which have the similar IPO with CMIP6-LE ensemble mean
 
 RPS = variance(obs_P(0:period))/variance(dim_avg_n(comb_P(0:period,:),1))
 
 adjust_AMV = dim_avg_n(comb_P(0:period,:),1)*RPS
 
 RMSE_ens = new((/227+90/),"float")
 
 do ensi = 0, 226+90
    ; RMSE_ens(ensi) = sqrt(avg((CMIP6_T(:, ensi)-adjust_AMV)^2))
    RMSE_ens(ensi) = escorc(comb_P(0:period, ensi),adjust_AMV)
 end do
  
 CMIP6_AMV_ord = dim_pqsort(RMSE_ens,1)
 
 ;print(CMIP6_IPO_ord)
 
; ****** write it into nc file *****************************	   
 ;Outfile = "Ann_Central_eastern_US_CMIP6_CESM2_LE_CRUT_prcp_1901_1950.nc"
 ;Outfile = "Ann_Central_eastern_US_CMIP6_CESM2_LE_CRUT_prcp_1901_1960.nc"
 ;Outfile = "Ann_Central_eastern_US_CMIP6_CESM2_LE_CRUT_prcp_1901_1970.nc"
 Outfile = "Ann_Central_eastern_US_CMIP6_CESM2_LE_CRUT_prcp_1901_1980.nc"

 ;Outfile = "Ann_UK_CMIP6_CESM2_LE_CRUT_prcp_1901_1950.nc"
 ;Outfile = "Ann_UK_CMIP6_CESM2_LE_CRUT_prcp_1901_1960.nc"
 ;Outfile = "Ann_UK_CMIP6_CESM2_LE_CRUT_prcp_1901_1970.nc"
 ;Outfile = "Ann_UK_CMIP6_CESM2_LE_CRUT_prcp_1901_1980.nc"

 ;Outfile = "Ann_western_Australia_CMIP6_CESM2_LE_CRUT_prcp_1901_1950.nc"
 ;Outfile = "Ann_western_Australia_CMIP6_CESM2_LE_CRUT_prcp_1901_1960.nc"
 ;Outfile = "Ann_western_Australia_CMIP6_CESM2_LE_CRUT_prcp_1901_1970.nc"
 ;Outfile = "Ann_western_Australia_CMIP6_CESM2_LE_CRUT_prcp_1901_1980.nc"

 system("/bin/rm -f "+Outfile)
 RPC_ID = addfile(Outfile,"c")

 RPC_ID ->CMIP6_AMV = comb_P	   
 RPC_ID ->GISS_AMV  = obs_P

 RPC_ID ->CMIP6_AMV_ord = CMIP6_AMV_ord  
  

  
  
  
  
  
