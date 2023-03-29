; Read in CMIP6 and Obs, calculate RPC for global continent region
; plot figures for global continent region 
; Method is described in Guo et al 2011 and Rosie et al 2014

; Inorder to make 5 years moving average, we need to extend the calculating period to 1900-2015
 beginy = 1951
 endy   = 2014
 term   = endy-beginy+1
 
 ;------ Create mask for 85N-60S, and cover continent --------------
 flsm     = addfile("/hosted/fores/CMIP6/Trend_ans/temperature_data/mask_2_5_degree.nc", "r") ;  [lat | 72] x [lon | 144]
 lsm_temp = flsm->lsm  
 lsm      = toint(flsm->lsm)
 lsm_size = dimsizes(lsm)
 
 lat_2_5  = flsm->lat
 lon_2_5  = flsm->lon
 
 mask_2_5_degree          = mask(lsm,((lsm.eq.1).or.(lsm.eq.2)),True)
 mask_2_5_degree(0:11,:)  = -127
 mask_2_5_degree(70:71,:) = -127

 mask_2_5_degree   = where(mask_2_5_degree.eq.2,1,mask_2_5_degree)  

 mask_2_5_degree!0 = "lat"
 mask_2_5_degree!1 = "lon"          

 mask_2_5_degree_reorder = mask_2_5_degree(lon|:,lat|:)
/;
 mask_all = new((/endy-beginy+1,lsm_size(0),lsm_size(1)/),typeof(mask_2_5_degree),mask_2_5_degree@_FillValue)
 do i = 0,endy-beginy
    mask_all(i,:,:) = mask_2_5_degree
 end do
 delete(i)
;/
; ****** read in CMIP6, GISS, CRUT temperature anomaly *********************
 cmipid     = addfile("CMIP6_GISS_CRUT_TAS_1901_2014_global.nc","r")
 tas_ann_nE = cmipid->tas_ann_nE(:,:,beginy-1901:endy-1901,:)
 ATS        = dimsizes(tas_ann_nE)
 
 GISS_ann    = cmipid->GISS_ann(beginy-1901:endy-1901,:,:)
 HadCRUT_ann = cmipid->HadCRUT_ann(beginy-1901:endy-1901,:,:) 

 mask_final =  mask_2_5_degree_reorder
 mask_final(:,:) = 0

 do loni = 0, ATS(0)-1
     do latj = 0, ATS(1)-1
        obs_temp = rm_single_dims(HadCRUT_ann(:,latj,loni))
        if ((.not.ismissing(mask_2_5_degree_reorder(loni,latj))).and.(any(obs_temp.ne.obs_temp@_FillValue))) then
            mask_final(loni,latj) = 1
        end if
     end do
 end do

 WID        = addfile("/hosted/fores/CMIP6/Trend_ans/temperature_data/RPC/area_2_5_degree.nc","r")
 carea_temp = WID->CAREA
 carea      = doubletofloat(carea_temp(lon|:,lat|:))
 carea      = where(mask_final.eq.1,carea,0)
 carea_tot  = sum(carea)
  
; ********calculate ratio of predictable components (RPC)*************************
 ;kind  = 0, it is CRUT, kind =1, it is GISS
 ; size is lon lat
 RPC_ID = addfile("/hosted/fores/CMIP6/Trend_ans/temperature_data/RPC/RPC_global_continent_1951_2014.nc","r")
 RPCO   = RPC_ID->RPC_Had  ;----<<------<<------<<------<<-----<<------<<-----------------
  
 RPC    = new((/ATS(0),ATS(1),1000/),"float",tas_ann_nE@_FillValue)
 Siglv  = new((/ATS(0),ATS(1)/),"float",tas_ann_nE@_FillValue)
  
  do loni = 0, ATS(0)-1
    do latj = 0, ATS(1)-1
      if (.not.ismissing(mask_2_5_degree(latj,loni))) then 
           
       ens_pool          = new((/ATS(2),221/),"float",tas_ann_nE@_FillValue)
       ens_pool(:,0:219) = rm_single_dims(tas_ann_nE(loni,latj,:,:))
       ens_pool(:,220)   = rm_single_dims(HadCRUT_ann(:,latj,loni))  ;----<<------<<------<<------<<-----<<------<<-----------------
		   
       tas_mean = rm_single_dims(dim_avg_n(ens_pool(:,0:219),1))	  
       Tmean    = avg(tas_mean)
       Var_sm   = variance(tas_mean) 		  

       Sum_nm = new((/220/),"float",tas_ann_nE@_FillValue) 
       do ensk = 0, 220-1
        Sum_nm(ensk) = avg((rm_single_dims(ens_pool(:,ensk))-Tmean)^2)
       end do

       Var_tm = avg(Sum_nm)

       do rantime = 0, 999
        ens_ind = generate_sample_indices(221,0)
        P_corre = escorc(tas_mean,rm_single_dims(ens_pool(:,ens_ind(220))))
			  
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

 
  system("/bin/rm -f RPC_global_continent_1951_2014_1000_sample_CRUT.nc")
  RPC_ID = addfile("RPC_global_continent_1951_2014_1000_sample_CRUT.nc","c")
  
  RPC_ID ->RPC_Had  = RPC   
  RPC_ID ->Siglv_Had  = Siglv 

  
  
  
  
  
  
  
  
  
  
  
  
  
  
