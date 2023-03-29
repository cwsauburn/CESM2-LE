; Read in CMIP6 and Obs, calculate RCP for global continent region
; plot figures for global continent region 
; Method is described in Guo et al 2011 and Rosie et al 2014

; Inorder to make 5 years moving average, we need to extend the calculating period to 1900-2015
 beginy = 1901
 endy   = 2014
 term = endy-beginy+1
 
 ;------ Create mask for 85N-60S, and cover continent --------------
 flsm     = addfile("/hosted/fores/CMIP6/Trend_ans/temperature_data/mask_2_5_degree.nc", "r") ;  [lat | 72] x [lon | 144]
 lsm_temp = flsm->lsm  
 lsm      = toint(flsm->lsm)
 lsm_size = dimsizes(lsm)
 
 lat_2_5  = flsm->lat
 lon_2_5  = flsm->lon
 
 mask_2_5_degree = mask(lsm,((lsm.eq.1).or.(lsm.eq.2)),True)
 mask_2_5_degree(0:11,:)  = -127
 mask_2_5_degree(70:71,:) = -127

 mask_2_5_degree = where(mask_2_5_degree.eq.2,1,mask_2_5_degree)  

 mask_2_5_degree!0 = "lat"
 mask_2_5_degree!1 = "lon"          
 ;mask_2_5_degree_reorder = mask_2_5_degree(lon|:,lat|:) 
 
; ******define constant value, read in prepared obs data *********************
  
  nMod = 23
  nEns = (/3, 3, 11, 3, 17, 5, 50, 2, 10, 10, 32, 10, 5, 8, 2, 4, 1, 4, 3, 5, 30, 1, 1/)
  
  /;
  norL = 1980
  norY = norL/12 
  time = ispan(1901,2013,1)
  ;/
  
  mseries = (/2, 4, 6, 8, 9, 10, 11, 20/)
   
  mod_name = new((/8/),"string")
  mod_name(0) = "CESM2"
  mod_name(1) = "CNRM-CM6-1"
  mod_name(2) = "CanESM5"
  mod_name(3) = "GISS-E2-1-G"
  mod_name(4) = "GISS-E2-1-H"
  mod_name(5) = "IPSL-CM6A-LR"
  mod_name(6) = "MIROC6"
  mod_name(7) = "NorCPM1" 

  cmipid          = addfile("CMIP6_GISS_CRUT_TAS_1901_2014_global.nc","r")
  tas_ann_nE_temp = cmipid->tas_ann_nE
  ATS             = dimsizes(tas_ann_nE_temp)
 
  GISS_ann    = cmipid->GISS_ann
  HadCRUT_ann = cmipid->HadCRUT_ann 

  /;
; temperature anomaly lon = 144*lat = 72*time = 114
  GISS_f_ID    = addfile("Ann_tas_GISS_1901_2014.nc","r")
  GISS_ann_all = GISS_f_ID->GISS 
  
; temperature anomaly lon = 144*lat = 72*time = 114
  HadCRUT_f_ID    = addfile("Ann_tas_CRUT_1901_2014.nc","r")
  HadCRUT_ann_all = HadCRUT_f_ID->CRUT
  
  GISS_ann    = runave_n(GISS_ann_all,5,0,0)
  HadCRUT_ann = runave_n(HadCRUT_ann_all,5,0,0)

  do mperiod = 0, endy-beginy;
     GISS_ann(mperiod, :, :)    = mask(rm_single_dims(GISS_ann_all(mperiod+GISS_be, :, :)), .not.ismissing(HadCRUT_ann_all(mperiod+Had_be, :, :)), True)
     HadCRUT_ann(mperiod, :, :) = mask(rm_single_dims(HadCRUT_ann_all(mperiod+Had_be, :, :)), .not.ismissing(GISS_ann_all(mperiod+GISS_be, :, :)), True)
  end do

  GISS_ann     = mask(GISS_ann,mask_all.eq.1,True) 
  HadCRUT_ann  = mask(HadCRUT_ann,mask_all.eq.1,True) 
  
  min_mask = mask_2_5_degree
  
  do mperiod = 0, endy-beginy;
     min_mask = mask(min_mask,.not.ismissing(HadCRUT_ann_all(mperiod+Had_be, :, :)), True) 
	 min_mask = mask(min_mask,.not.ismissing(GISS_ann_all(mperiod+GISS_be, :, :)), True) 
  end do	

  RCP_ID ->min_mask = min_mask
  ;/
  
  do mse_ind = 0,7
  
  enindex    = sum(nEns(0:mseries(mse_ind)))-1
  beindex    = enindex-nEns(mseries(mse_ind))+1
    
  tas_ann_nE = tas_ann_nE_temp(:,:,:,beindex:enindex) 
  ATS        = dimsizes(tas_ann_nE)
  
; ********calculate ratio of predictable components (RPC)*************************
  
  RPC_Had  = new((/ATS(0),ATS(1)/),"float",tas_ann_nE@_FillValue)
  RPC_GISS = new((/ATS(0),ATS(1)/),"float",tas_ann_nE@_FillValue)
  
  do loni = 0, ATS(0)-1
     do latj = 0, ATS(1)-1
       if (.not.ismissing(mask_2_5_degree(latj,loni))) then 
          tas_mean = rm_single_dims(dim_avg_n(tas_ann_nE(loni,latj,:,:),1))
	  Tmean    = avg(tas_mean)
	  Var_sm   = variance(tas_mean)

          P_corre0 = escorc(tas_mean,rm_single_dims(HadCRUT_ann(:,latj,loni)))
          P_corre1 = escorc(tas_mean,rm_single_dims(GISS_ann(:,latj,loni)))
		   
          Sum_nm = new((/ATS(3)/),"float",tas_ann_nE@_FillValue)
		   
          do ensk = 0, ATS(3)-1 
            Sum_nm(ensk) = avg((rm_single_dims(tas_ann_nE(loni,latj,:,ensk))-Tmean)^2)
          end do
		   
	  Var_tm = avg(Sum_nm)
		   
	  RPC_Had(loni,latj)  = P_corre0/sqrt(Var_sm/Var_tm)   
	  RPC_GISS(loni,latj) = P_corre1/sqrt(Var_sm/Var_tm) 
       end if
    end do
  end do

  out_file = "RPC_global_continent_1901_2014_"+mod_name(mse_ind)+".nc"  
  
  system("/bin/rm -f "+out_file)
  RPC_ID = addfile(out_file,"c")
  
  RPC_ID ->RPC_Had  = RPC_Had 
  RPC_ID ->RPC_GISS = RPC_GISS 
  
  delete(tas_ann_nE)
  delete(RPC_Had)
  delete(RPC_GISS)
  delete(tas_mean)
  delete(Sum_nm)

 end do 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
