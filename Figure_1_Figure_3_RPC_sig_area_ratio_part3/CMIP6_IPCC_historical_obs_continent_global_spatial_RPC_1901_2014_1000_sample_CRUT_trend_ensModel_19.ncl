; Read in CMIP6 and Obs, calculate RPC for global continent region
; plot figures for global continent region 
; Method is described in Guo et al 2011 and Rosie et al 2014

; Inorder to make 5 years moving average, we need to extend the calculating period to 1900-2015
 beginy = 1901
 endy   = 2014
 term = endy-beginy+1
 
 ;------ Create mask for 85N-60S, and cover continent --------------
 flsm     = addfile("~/mask_2_5_degree.nc", "r") ;  [lat | 72] x [lon | 144]
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
  
  mseries = (/6, 10/)
  mse_ind = 0
  
  ; read in all data together, lon, lat, year, ensemble 144*72*114*(132+52+37=221)
  ann_tas1 = addfile("~/Ann_tas_CMIP6_1_1901_2014.nc","r")
  ann_tas2 = addfile("~/Ann_tas_CMIP6_2_1901_2014.nc","r")
  ann_tas3 = addfile("~/Ann_tas_CMIP6_3_1901_2014.nc","r")
  
  tas_ann_nE1      = ann_tas1->CMIP6
  tas_ann_nE2_temp = ann_tas2->CMIP6
  tas_ann_nE3      = ann_tas3->CMIP6
  
  TANT_size               = dimsizes(tas_ann_nE2_temp)
  tas_ann_nE2             = new((/TANT_size(0),TANT_size(1),TANT_size(2),TANT_size(3)-1/),typeof(tas_ann_nE2_temp),tas_ann_nE2_temp@_FillValue)
  tas_ann_nE2(:,:,:,0:6)  = tas_ann_nE2_temp(:,:,:,0:6) 
  tas_ann_nE2(:,:,:,7:TANT_size(3)-2)   = tas_ann_nE2_temp(:,:,:,8:TANT_size(3)-1)
    
  tas_ann_nE_temp = new((/TANT_size(0),TANT_size(1),term,220/),typeof(tas_ann_nE1),tas_ann_nE1@_FillValue)
  
  tas_ann_nE_temp(:,:,:,0:131)   = runave_n(tas_ann_nE1,5,0,2) 
  tas_ann_nE_temp(:,:,:,132:182) = runave_n(tas_ann_nE2,5,0,2) 
  tas_ann_nE_temp(:,:,:,183:219) = runave_n(tas_ann_nE3,5,0,2)
  
  enindex    = sum(nEns(0:mseries(mse_ind)))-1
  beindex    = enindex-nEns(mseries(mse_ind))+1
    
  tas_ann_nE = tas_ann_nE_temp(:,:,:,beindex:enindex)
  ATS        = dimsizes(tas_ann_nE)

; ******** Add Observation, it is global gridded data **************************

; Add GISS 2.5*2.5 regridded data
; short tempanomaly [time | 1671] x [lat | 72] x [lon | 144] 
; base_period = '1951-1980'
; Time period(188001-201903)

  GISS_f_ID = addfile("~/GISS_global_2.5.nc","r")
  GISS      = GISS_f_ID->GISS 

 ; Add cru_ts4.03.1901.2018.tmp.dat.nc 0.5*0.5 regridded data
 ; short tempanomaly time 1416 x lat 360 x lon 720
 ; Time period(190101-201812)
		
  HadCRUT_f_ID = addfile("~/cru_ts4.03.1901.2018.tmp.dat.nc","r")
  tmp1 = (HadCRUT_f_ID->tmp)
  tmp  = tmp1
  tmp(:,:,0:359)   = tmp1(:,:,360:719)
  tmp(:,:,360:719) = tmp1(:,:,0:359)
  
  lat_CRU = tmp&lat
  lon_CRU = tmp&lon
  lon_CRU = where(lon_CRU.lt.0, lon_CRU+360, lon_CRU)
  
  gwt  = latRegWgt(tmp&lat, typeof(tmp), 0) 
  opt1 = True
  opt1@critpc = 75
  HadCRUT = area_hi2lores_Wrap(lon_CRU, lat_CRU, tmp, True, gwt, lon_2_5, lat_2_5, opt1)  

  GISS_ann_all    = new((/139,lsm_size(0),lsm_size(1)/),"float")
  HadCRUT_ann_all = new((/118,lsm_size(0),lsm_size(1)/),"float")
  
  GISS_ann_all    = month_to_annual(GISS(0:(2019-1880)*12-1,:,:),1)
  HadCRUT_ann_all = month_to_annual(HadCRUT(:,:,:),1)
  
  GISS_ann    = new((/endy-beginy+1,lsm_size(0),lsm_size(1)/),"float")
  HadCRUT_ann = new((/endy-beginy+1,lsm_size(0),lsm_size(1)/),"float")  
  
  ; GISS begin with 1901-1880=21, 
  GISS_be = beginy-1880
  
  ; HadCRUT4 begin with 1901-1850=51
  Had_be = beginy-1901
  
  GISS_ann    = runave_n(GISS_ann_all(GISS_be:(GISS_be+endy-beginy), :, :),5,0,0)
  HadCRUT_ann = runave_n(HadCRUT_ann_all(Had_be:(Had_be+endy-beginy), :, :),5,0,0)
  
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
  WID        = addfile("~/area_2_5_degree.nc","r")
  carea_temp = WID->CAREA
  carea      = doubletofloat(carea_temp(lon|:,lat|:)) 
  carea      = where(mask_final.eq.1,carea,0)
  carea_tot  = sum(carea)
      
; ********calculate ratio of predictable components (RPC)*************************

  RPC    = new((/ATS(0),ATS(1),1000/),"float",tas_ann_nE@_FillValue)
  Siglv  = new((/ATS(0),ATS(1)/),"float",tas_ann_nE@_FillValue)

  Siglv_final  = new((/50,500/),"float",tas_ann_nE@_FillValue)
  
  flag = 0;
  do iterk = 0, 499
     do ensi = 46, 47
		 ens_sel = generate_sample_indices(50,1)	
		 
		 flag = flag +1
		 print(flag) 
		  
         do loni = 0, ATS(0)-1
             do latj = 0, ATS(1)-1
                if (mask_final(loni,latj).eq.1) then 

                   ens_pool = new((/ATS(2),ensi+1/),"float",tas_ann_nE@_FillValue)
                   ens_pool(:,0:ensi-1) = rm_single_dims(tas_ann_nE(loni,latj,:,ens_sel(0:ensi-1)))
                   ens_pool(:,ensi)     = rm_single_dims(HadCRUT_ann(:,latj,loni))  ;----<<------<<------<<------<<-----<<------<<-----------------

                   tas_mean = rm_single_dims(dim_avg_n(ens_pool(:,0:ensi-1),1))
                   Tmean    = avg(tas_mean)
                   Var_sm   = variance(tas_mean)	

                   Sum_nm = new((/ensi/),"float",tas_ann_nE@_FillValue)
                   do ensk = 0, ensi-1
                      Sum_nm(ensk) = avg((rm_single_dims(ens_pool(:,ensk))-Tmean)^2)
                   end do

                   Var_tm = avg(Sum_nm)
                   delete(Sum_nm)				   

                   do rantime = 0, 999
                      ens_ind = generate_sample_indices(ensi+1,1)
                      P_corre = escorc(tas_mean,rm_single_dims(ens_pool(:,ens_ind(ensi))))
                      delete(ens_ind)
					  
                      RPC(loni,latj,rantime)   = P_corre/sqrt(Var_sm/Var_tm)   
                   end do 
                   sortRPC = RPC(loni,latj,:)		   
                   qsort(sortRPC)

                   CV  = new((/5/),typeof(sortRPC),sortRPC@_FillValue) 
                   CV(0) = sortRPC(round(0.025*1000,3)-1) 
                   CV(1) = sortRPC(round(0.05*1000,3)-1)  
                   CV(2) = sortRPC(round(0.95*1000,3)-1) 
                   CV(3) = sortRPC(round(0.975*1000,3)-1)  
				   
				   P_corre0 = escorc(tas_mean,rm_single_dims(HadCRUT_ann(:,latj,loni)))
				   RPCO     = P_corre0/sqrt(Var_sm/Var_tm)
                   CV(4)    = RPCO
				   
                   delete(tas_mean) 
				   delete(ens_pool)
				   
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
          Siglv_final(ensi-1,iterk) = STP_LTP_Had/carea_tot          

		  Siglv(:,:) = tas_ann_nE@_FillValue
          RPC(:,:,:) = tas_ann_nE@_FillValue  		  
      end do
  end do
 
  system("/bin/rm -f RPC_global_continent_1901_2014_1000_sample_CRUT_CanESM5_trend_19.nc")
  RPC_ID = addfile("RPC_global_continent_1901_2014_1000_sample_CRUT_CanESM5_trend_19.nc","c")
  RPC_ID ->Siglv_final_Had  = Siglv_final 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
