; Read in CMIP6 and Obs, calculate RPC for global continent region
; plot figures for global continent region 
; Method is described in Guo et al 2011 and Rosie et al 2014

; Inorder to make 5 years moving average, we need to extend the calculating period to 1900-2015
 beginy = 1901
 endy   = 2014
 term   = endy-beginy+1
 
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
 
 mask_all = new((/endy-beginy+1,lsm_size(0),lsm_size(1)/),typeof(mask_2_5_degree),mask_2_5_degree@_FillValue)
 do i = 0,endy-beginy
    mask_all(i,:,:) = mask_2_5_degree
 end do
 delete(i)
 
; ******define constant value, read in prepared obs data *********************

 nMod = 23
 nEns = (/3, 3, 11, 3, 17, 5, 50, 2, 10, 10, 32, 10, 5, 8, 2, 4, 1, 4, 3, 5, 30, 1, 1/)
  
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
  
  do mse_ind = 0,7

  enindex    = sum(nEns(0:mseries(mse_ind)))-1
  beindex    = enindex-nEns(mseries(mse_ind))+1
    
  tas_ann_nE = tas_ann_nE_temp(:,:,:,beindex:enindex)
  ATS        = dimsizes(tas_ann_nE)

; ********calculate ratio of predictable components (RPC)*************************
  ;kind  = 0, it is CRUT, kind =1, it is GISS
  ;lon lat
  RPC_file = "/hosted/fores/CMIP6/Trend_ans/temperature_data/RPC/RPC_global_continent_1901_2014_"+mod_name(mse_ind)+".nc"
  
  RPC_ID = addfile(RPC_file,"r")
  RPCO   = RPC_ID->RPC_Had  ;----<<------<<------<<------<<-----<<------<<-----------------
  ;RPCO  = RPC_ID->RPC_GISS
  
  RPC    = new((/ATS(0),ATS(1),1000/),"float",tas_ann_nE@_FillValue)
  Siglv  = new((/ATS(0),ATS(1)/),"float",tas_ann_nE@_FillValue)
 
  do loni = 0, ATS(0)-1
     do latj = 0, ATS(1)-1
      if (mask_final(loni,latj).eq.1) then 
           
        ens_pool = new((/ATS(2),nEns(mseries(mse_ind))+1/),"float",tas_ann_nE@_FillValue)
        ens_pool(:,0:nEns(mseries(mse_ind))-1) = rm_single_dims(tas_ann_nE(loni,latj,:,:))
        ens_pool(:,nEns(mseries(mse_ind)))     = rm_single_dims(HadCRUT_ann(:,latj,loni))  ;----<<------<<------<<------<<-----<<------<<-----------------
        ;ens_pool(:,nEns(mseries(mse_ind)))   = rm_single_dims(GISS_ann(:,latj,loni))
		   
	tas_mean = rm_single_dims(dim_avg_n(ens_pool(:,0:nEns(mseries(mse_ind))-1),1))
        Tmean    = avg(tas_mean)
 	Var_sm   = variance(tas_mean)
	   
        Sum_nm = new((/nEns(mseries(mse_ind))/),"float",tas_ann_nE@_FillValue)
        do ensk = 0, nEns(mseries(mse_ind))-1
         Sum_nm(ensk) = avg((rm_single_dims(ens_pool(:,ensk))-Tmean)^2)
        end do

        Var_tm = avg(Sum_nm)

        do rantime = 0, 999
	  ens_ind = generate_sample_indices(nEns(mseries(mse_ind))+1,0)
          P_corre = escorc(tas_mean,rm_single_dims(ens_pool(:,ens_ind(nEns(mseries(mse_ind))))))
			  
          RPC(loni,latj,rantime) = P_corre/sqrt(Var_sm/Var_tm)   
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

  out_file = "RPC_global_continent_1901_2014_1000_sample_CRUT_"+mod_name(mse_ind)+".nc"
  system("/bin/rm -f "+out_file)
  RPC_ID = addfile(out_file,"c")
  
  RPC_ID ->RPC_Had  = RPC 
  ;RPC_ID ->RPC_GISS = RPC
  
  RPC_ID ->Siglv_Had  = Siglv 
  ;RPC_ID ->Siglv_GISS = Siglv 
  
  delete(tas_ann_nE)
  delete(RPC_file)
  delete(RPCO)
  delete(RPC)
  delete(Siglv) 
  delete(ens_pool)
  delete(tas_mean)
  delete(ens_ind)
  delete(out_file)
  delete(Sum_nm)
 end do 
  
  
  
  
  
  
  
  
  
  
  
  
