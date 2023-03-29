; Read in CMIP6 and Obs, calculate RCP for global continent region
; plot figures for global continent region 
; Method is described in Guo et al 2011 and Rosie et al 2014

; Inorder to make 5 years moving average, we need to extend the calculating period to 1900-2015
 beyear = 1901
 edyear = 2014
 term   = edyear-beyear+1
 
 ;------ Create mask for 85N-60S, and cover continent --------------
 flsm     = addfile("mask_2_5_degree_rg2CESM2.nc", "r") ;  [lat | 192] x [lon | 288]
 lsm_temp = flsm->lsm
 lsm      = toint(lsm_temp)
 lsm_size = dimsizes(lsm)
 lat_2_5  = lsm_temp&lat
 lon_2_5  = lsm_temp&lon

 mask_2_5_degree          = mask(lsm,((lsm.eq.1).or.(lsm.eq.2)),True)
 mask_2_5_degree(0:31,:)    = -127   ;60S
 mask_2_5_degree(186:191,:) = -127   ;85N
 mask_2_5_degree          = where(mask_2_5_degree.eq.2,1,mask_2_5_degree)
 mask_2_5_degree!0        = "lat"
 mask_2_5_degree!1        = "lon"
 
 mask_2_5_degree = mask_2_5_degree(::-1,:)
 
; ******define constant value, read in prepared obs data *********************
 CESMID     = addfile("Ann_tas_CESM2_LE_global_1901_2014_anomaly.nc","r") 
 tas_ann_nE = CESMID->tas_ann_nE
   
 ; time 114, lat 192, lon 288, ens 90
 ATS = dimsizes(tas_ann_nE) 
 
; ******** Add Observation, it is global gridded data **************************
 ; lon 288, lat 192, year 114
 CRUT_f_ID = addfile("Ann_tas_CRUT_1901_2014_1_degree.nc","r")
 CRUT_ann0 = CRUT_f_ID->tas 
 HadCRUT_ann = runave_n(CRUT_ann0(:, :, :),5,0,2)
  
 mask_final = mask_2_5_degree
 mask_final(:,:) = 0
  
 do loni = 0, ATS(2)-1
    do latj = 0, ATS(1)-1
       obs_temp = rm_single_dims(HadCRUT_ann(loni,latj,:))
	   if ((.not.ismissing(mask_2_5_degree(latj,loni))).and.(any(obs_temp.ne.obs_temp@_FillValue))) then 
		   mask_final(latj,loni) = 1
	   end if
	end do
 end do	 		
    
; ******** Read in weight, it is global gridded data *****************************
; size is lat 192 * lon 288
  WID        = addfile("CLM4.0CN.h0.fixed.nc","r")
  carea_temp = WID->area
  carea_temp@_FillValue = 1e+36
  carea_temp = where(mask_final.eq.1,carea_temp,carea_temp@_FillValue)
  carea      = doubletofloat(carea_temp(lon|:,lat|:)) 
  
  carea_tot  = sum(carea)

; ********calculate ratio of predictable components (RPC)*************************
  
  ;RPC_Had  = new((/ATS(2),ATS(1),90/),"float",tas_ann_nE@_FillValue)
  P_corre0 = new((/90,500/),"float",tas_ann_nE@_FillValue)
  denomin  = new((/90,500/),"float",tas_ann_nE@_FillValue)
  RPC      = new((/90,500/),"float",tas_ann_nE@_FillValue)
  
  P_corre0_temp = new((/ATS(2),ATS(1)/),"float",tas_ann_nE@_FillValue)
  denomin_temp  = new((/ATS(2),ATS(1)/),"float",tas_ann_nE@_FillValue)
  RPC_Had_temp  = new((/ATS(2),ATS(1)/),"float",tas_ann_nE@_FillValue)
  
  do iterk = 0, 499
    do ensi = 84, 84
       print(ensi)
	   pool = generate_sample_indices(90,1)
	   unf  = pool(0:ensi-1)
		
		do loni = 0, ATS(2)-1
		   do latj = 0, ATS(1)-1 
		   
		      obs_temp = rm_single_dims(HadCRUT_ann(loni,latj,:))
		   
			   if (mask_final(latj,loni).eq.1) then 

			     tas_ann_nE_temp = new((/ATS(0),ensi/),"float",tas_ann_nE@_FillValue)
			     tas_ann_nE_temp(:,:) = 0
				 
			     do unfi = 0, ensi-1
			        tas_ann_nE_temp(:,unfi) = rm_single_dims(tas_ann_nE(:,latj,loni,unf(unfi)))
			     end do

		         tas_mean = rm_single_dims(dim_avg_n(tas_ann_nE_temp,1))
			     delete(tas_ann_nE_temp)
			  
		         Tmean    = avg(tas_mean)
		         Var_sm   = variance(tas_mean)

                 P_corre0_temp(loni,latj) = escorc(tas_mean,obs_temp)

		         Sum_nm = new((/ensi/),"float",tas_ann_nE@_FillValue)
		   
		         do ensk = 0, ensi-1 
                    Sum_nm(ensk) = avg((rm_single_dims(tas_ann_nE(:,latj,loni,unf(ensk)))-Tmean)^2)
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
	   P_corre0(ensi-1,iterk) = sum(P_corre0_temp)
	   denomin(ensi-1,iterk)  = sum(denomin_temp)
	   RPC(ensi-1,iterk)      = sum(RPC_Had_temp)
       	
       P_corre0_temp(:,:) = tas_ann_nE@_FillValue
       denomin_temp(:,:)  = tas_ann_nE@_FillValue
       RPC_Had_temp(:,:)  = tas_ann_nE@_FillValue	
	   	   
	   delete(pool)
	   delete(unf)
	   
	end do
  end do
  
  system("/bin/rm -f Ann_tas_CESM2_LE_RPC_global_continent_1901_2014_trend_18.nc")
  RPC_ID = addfile("Ann_tas_CESM2_LE_RPC_global_continent_1901_2014_trend_18.nc","c")
  
  RPC_ID ->P_corre0 = P_corre0 
  RPC_ID ->denomin  = denomin
  RPC_ID ->RPC      = RPC

 
  
  
  
  
  
  
