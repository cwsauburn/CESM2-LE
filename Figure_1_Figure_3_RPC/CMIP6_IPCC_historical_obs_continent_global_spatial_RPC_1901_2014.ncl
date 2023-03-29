; Read in CMIP6 and Obs, calculate RCP for global region
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
  
; ****** read in CMIP6, GISS, CRUT temperature anomaly *********************
 cmipid     = addfile("CMIP6_GISS_CRUT_TAS_1901_2014_global.nc","r")
 tas_ann_nE = cmipid->tas_ann_nE
 ATS        = dimsizes(tas_ann_nE)
 
 GISS_ann    = cmipid->GISS_ann
 HadCRUT_ann = cmipid->HadCRUT_ann 
  
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
  
  system("/bin/rm -f RPC_global_continent_1901_2014.nc")
  RPC_ID = addfile("RPC_global_continent_1901_2014.nc","c")
  
  RPC_ID ->RPC_Had  = RPC_Had 
  RPC_ID ->RPC_GISS = RPC_GISS 
  
  
  
  
  

  
  
  
  
  
  
  
  
  
  
