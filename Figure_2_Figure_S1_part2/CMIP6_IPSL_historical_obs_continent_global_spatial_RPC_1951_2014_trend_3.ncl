; Read in CMIP6 and Obs, calculate RCP for global continent region
; plot figures for global continent region 
; Method is described in Guo et al 2011 and Rosie et al 2014

; Inorder to make 5 years moving average, we need to extend the calculating period to 1900-2015
 beginy = 1951
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
 mask_2_5_degree_reorder = mask_2_5_degree(lon|:,lat|:) 

 /;
 mask_all = new((/endy-beginy+1,lsm_size(0),lsm_size(1)/),typeof(mask_2_5_degree),mask_2_5_degree@_FillValue)
 do i = 0,endy-beginy
    mask_all(i,:,:) = mask_2_5_degree
 end do
 delete(i)
 ;/

; ******define constant value, read in prepared obs data *********************
  
  nMod = 23
  nEns = (/3, 3, 11, 3, 17, 5, 50, 2, 10, 10, 32, 10, 5, 8, 2, 4, 1, 4, 3, 5, 30, 1, 1/)

  mseries = (/6, 10, 20/)
  mse_ind = 1

  enindex    = sum(nEns(0:mseries(mse_ind)))-1
  beindex    = enindex-nEns(mseries(mse_ind))+1

  fileID      = addfile("/hosted/fores/CMIP6/Trend_ans/temperature_data/RPC/CMIP6_GISS_CRUT_TAS_1901_2014_global.nc","r")
  tas_ann_nE  = fileID->tas_ann_nE(:,:,beginy-1901:endy-1901,beindex:enindex)
  HadCRUT_ann = fileID->HadCRUT_ann(beginy-1901:endy-1901,:,:)

  ATS         = dimsizes(tas_ann_nE)
  
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

; ******** Read in weight, it is global gridded data *****************************
; size is lat 72 * lon 144
  WID        = addfile("/hosted/fores/CMIP6/Trend_ans/temperature_data/RPC/area_2_5_degree.nc","r")
  carea_temp = WID->CAREA
  carea      = doubletofloat(carea_temp(lon|:,lat|:)) 
  carea      = where(mask_final.eq.1,carea,0)
  carea_tot  = sum(carea)
    
; ********calculate ratio of predictable components (RPC)*************************

  P_corre0_temp = new((/ATS(0),ATS(1)/),"float",tas_ann_nE@_FillValue)
  denomin_temp  = new((/ATS(0),ATS(1)/),"float",tas_ann_nE@_FillValue)
  RPC_Had_temp  = new((/ATS(0),ATS(1)/),"float",tas_ann_nE@_FillValue)

  step = ispan(2,32,2)

  ;RPC_Had  = new((/ATS(0),ATS(1)/),"float",tas_ann_nE@_FillValue)
  P_corre0 = new((/dimsizes(step),500/),"float",tas_ann_nE@_FillValue)
  denomin  = new((/dimsizes(step),500/),"float",tas_ann_nE@_FillValue)
  RPC      = new((/dimsizes(step),500/),"float",tas_ann_nE@_FillValue)  
  
  do iterk = 0, 499
    print(iterk)
    do ensi_index = 0, 9   
      ensi    = step(ensi_index)  
      pool = generate_sample_indices(32,1)
      unf  = pool(0:ensi-1) 
		
       do loni = 0, ATS(0)-1
         do latj = 0, ATS(1)-1
	    if (mask_final(loni,latj).eq.1) then 

              tas_ann_nE_temp = new((/ATS(2),ensi/),"float",tas_ann_nE@_FillValue)
	      tas_ann_nE_temp(:,:) = 0
				 
	      do unfi = 0, ensi-1
	       tas_ann_nE_temp(:,unfi) = rm_single_dims(tas_ann_nE(loni,latj,:,unf(unfi)))
	      end do
		  
	      tas_mean = rm_single_dims(dim_avg_n(tas_ann_nE_temp,1))
	      delete(tas_ann_nE_temp)
			  
	      Tmean    = avg(tas_mean)
	      Var_sm   = variance(tas_mean)

              P_corre0_temp(loni,latj) = escorc(tas_mean,rm_single_dims(HadCRUT_ann(:,latj,loni)))

              Sum_nm = new((/ensi/),"float",tas_ann_nE@_FillValue)
		   
	      do ensk = 0, ensi-1 
                 Sum_nm(ensk) = avg((rm_single_dims(tas_ann_nE(loni,latj,:,unf(ensk)))-Tmean)^2)
	      end do
		   
	      Var_tm = avg(Sum_nm)
	      delete(Sum_nm)
				 
	      denomin_temp(loni,latj) = sqrt(Var_sm/Var_tm)
		   
	      RPC_Had_temp(loni,latj)  = P_corre0_temp(loni,latj)/denomin_temp(loni,latj) 

              P_corre0_temp(loni,latj) = P_corre0_temp(loni,latj)*(carea(loni,latj)/carea_tot)
	      denomin_temp(loni,latj)  = denomin_temp(loni,latj)*(carea(loni,latj)/carea_tot)
	      RPC_Had_temp(loni,latj)  = RPC_Had_temp(loni,latj)*(carea(loni,latj)/carea_tot)
				 
	    end if	
        end do
       end do 
       P_corre0(ensi_index,iterk) = sum(P_corre0_temp)
       denomin(ensi_index,iterk)  = sum(denomin_temp)
       RPC(ensi_index,iterk)      = sum(RPC_Had_temp)

       P_corre0_temp(:,:) = tas_ann_nE@_FillValue
       denomin_temp(:,:)  = tas_ann_nE@_FillValue
       RPC_Had_temp(:,:)  = tas_ann_nE@_FillValue	
  
       delete(pool)
       delete(unf)
	   
   end do
  end do
  
  system("/bin/rm -f Ann_tas_CMIP6_IPSL_RPC_global_continent_1951_2014_trend_3.nc")
  RPC_ID = addfile("Ann_tas_CMIP6_IPSL_RPC_global_continent_1951_2014_trend_3.nc","c")
  
  RPC_ID ->P_corre0  = P_corre0 
  RPC_ID ->denomin   = denomin
  RPC_ID ->RPC       = RPC
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
