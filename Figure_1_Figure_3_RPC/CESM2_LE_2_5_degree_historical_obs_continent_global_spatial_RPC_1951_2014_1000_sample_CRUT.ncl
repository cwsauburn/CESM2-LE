; Read in CMIP6 and Obs, calculate RCP for global continent region
; plot figures for global continent region 
; Method is described in Guo et al 2011 and Rosie et al 2014

; Inorder to make 5 years moving average, we need to extend the calculating period to 1900-2015
 beyear = 1951
 edyear = 2014
 term   = edyear-beyear+1
 
 ;------ Create mask for 85N-60S, and cover continent --------------
 flsm     = addfile("mask_2_5_degree.nc", "r") ;  [lat | 72] x [lon | 144]
 lsm_temp = flsm->lsm
 lsm      = toint(flsm->lsm)
 lsm_size = dimsizes(lsm)

 lat_2_5  = flsm->lat
 lon_2_5  = flsm->lon

 mask_2_5_degree          = mask(lsm,((lsm.eq.1).or.(lsm.eq.2)),True)
 mask_2_5_degree(0:11,:)  = -127
 mask_2_5_degree(70:71,:) = -127

 mask_2_5_degree          = where(mask_2_5_degree.eq.2,1,mask_2_5_degree)

 mask_2_5_degree!0 = "lat"
 mask_2_5_degree!1 = "lon"  
 
 mask_2_5_degree_reorder = mask_2_5_degree(lon|:,lat|:) 
  
;------ Read in CESM2-LE and interpolate it into 2.5 degree ------------------------ 
 ; size is lat 72, lon 144, year 114, ens 90
 NID              = addfile("Ann_tas_CESM2_LE_global_1901_2014_2_5_degree.nc","r")
 tas_ann_nE_temp  = NID->tas_ann_nE(:,:,beyear-1901:edyear-1901,:)
 tas_ann_nE_temp2 = tas_ann_nE_temp(ncl2|:,lat|:,lon|:,ncl3|:)
 tas_ann_nE       = runave_n(tas_ann_nE_temp2,5,0,0)

 ; time 114, lat , lon , ens 90
 ATS = dimsizes(tas_ann_nE) 

; ******define constant value, read in prepared obs data ************************ 
; ******** Add Observation, it is global gridded data *************************** 
  NIDT            = addfile("Ann_tas_HadCRUT_global_1901_2018_2_5_degree.nc","r")
  HadCRUT_ann_all = NIDT->HadCRUT_ann_all
  
  CRUT_ann        = new((/edyear-beyear+1,ATS(1),ATS(2)/),"float")

  ; HadCRUT4 begin with 1901-1850=51
  Had_be   = beyear-1901
  ; time*lat*lon
  CRUT_ann = runave_n(HadCRUT_ann_all(Had_be:(Had_be+edyear-beyear), :, :),5,0,0)
  
  mask_final      =  mask_2_5_degree_reorder
  mask_final(:,:) = 0

  do loni = 0, ATS(2)-1
    do latj = 0, ATS(1)-1
      obs_temp = rm_single_dims(CRUT_ann(:,latj,loni))
      if ((.not.ismissing(mask_2_5_degree_reorder(loni,latj))).and.(any(obs_temp.ne.obs_temp@_FillValue))) then
         mask_final(loni,latj) = 1
      end if
    end do
  end do	

; ******** Read in weight, it is global gridded data *****************************
; size is lat 72 * lon 144
  WID        = addfile("~/area_2_5_degree.nc","r")
  carea_temp = WID->CAREA
  carea      = doubletofloat(carea_temp(lon|:,lat|:))
  carea      = where(mask_final.eq.1,carea,0)
  carea_tot  = sum(carea)

; ********calculate ratio of predictable components (RPC)*************************
 ;kind  = 0, it is CRUT, kind =1, it is GISS
 ; size is lon lat
 RPC_ID = addfile("Ann_tas_CESM2_LE_2_5_degree_RPC_global_continent_1951_2014.nc","r")
 RPCO   = RPC_ID->RPC_Had  ;----<<------<<------<<------<<-----<<------<<-----------------
  
 RPC    = new((/144,72,1000/),"float",tas_ann_nE@_FillValue)
 Siglv  = new((/144,72/),"float",tas_ann_nE@_FillValue)
  
  do loni = 0, 144-1
    do latj = 0, 72-1
      if (.not.ismissing(mask_2_5_degree(latj,loni))) then 
           
       ens_pool         = new((/64,91/),"float",tas_ann_nE@_FillValue)
       ens_pool(:,0:89) = rm_single_dims(tas_ann_nE(:,latj,loni,:))
       ens_pool(:,90)   = (/rm_single_dims(CRUT_ann(:,latj,loni))/)  ;----<<------<<------<<------<<-----<<------<<-----------------
		   
       tas_mean = rm_single_dims(dim_avg_n(ens_pool(:,0:89),1))	  
       Tmean    = avg(tas_mean)	
       Var_sm   = variance(tas_mean)  

       Sum_nm = new((/90/),"float",tas_ann_nE@_FillValue)
       do ensk = 0, 89
         Sum_nm(ensk) = avg((rm_single_dims(ens_pool(:,ensk))-Tmean)^2)
       end do

       Var_tm = avg(Sum_nm)
	    
       do rantime = 0, 999
        ens_ind = generate_sample_indices(91,0)
        P_corre = escorc(tas_mean,rm_single_dims(ens_pool(:,ens_ind(90))))
			  
        RPC(loni,latj,rantime)   = P_corre/sqrt(Var_sm/Var_tm)
       end do 

       sortRPC = RPC(loni,latj,:)		   
       qsort(sortRPC)

       CV    = new((/5/),typeof(sortRPC),sortRPC@_FillValue) 
       CV(0) = sortRPC(round(0.025*1000,3)-1) 
       CV(1) = sortRPC(round(0.05*1000,3)-1)  
       CV(2) = sortRPC(round(0.95*1000,3)-1) 
       CV(3) = sortRPC(round(0.975*1000,3)-1)  
       CV(4) = RPCO(loni,latj)
		   
       if (.not.any(ismissing(CV))) then
		   
         if (CV(4).le.CV(0)) then
            Siglv(loni,latj) = 0.025
         else if (CV(4).le.CV(1)) then			 
            Siglv(loni,latj) = 0.05
         else if (CV(4).ge.CV(3)) then
	    Siglv(loni,latj) = 0.025
	 else if (CV(4).ge.CV(2)) then
	    Siglv(loni,latj) = 0.05  
	 end if 
         end if
	 end if
         end if
         end if
         delete(CV)	   
       end if
    end do
  end do
  
  STP_LTP_Had   = sum(mask(carea, Siglv,0.025))
  print(STP_LTP_Had/carea_tot)  
	   
  system("/bin/rm -f CESM2_LE_2_5_degree_RPC_global_continent_1951_2014_1000_sample_CRUT.nc")
  RPC_ID = addfile("CESM2_LE_2_5_degree_RPC_global_continent_1951_2014_1000_sample_CRUT.nc","c")
  
  RPC_ID ->RPC_Had  = RPC   
  RPC_ID ->Siglv_Had  = Siglv 

  
  
  
  
  
  
  
  
  
  
  
  
  
  
