; Read in CESM-LE and Obs, calculate RPC for global continent region
; plot figures for global continent region 
; Method is described in Guo et al 2011 and Rosie et al 2014

 beyear = 1901
 edyear = 2014
 term   = edyear-beyear+1
 
 ;------ Create mask for 85N-60S, and cover continent -----------------------------
 flsm     = addfile("mask_2_5_degree.nc", "r") ;  [lat | 72] x [lon | 144]
 lsm_temp = flsm->lsm
 lsm      = toint(flsm->lsm)
 lsm_size = dimsizes(lsm)
 
 lat_2_5  = flsm->lat
 lon_2_5  = flsm->lon

 mask_2_5_degree          = mask(lsm,((lsm.eq.1).or.(lsm.eq.2)),True)
 mask_2_5_degree(0:11,:)  = -127   ;60S
 mask_2_5_degree(70:71,:) = -127   ;85N
 
 mask_2_5_degree          = where(mask_2_5_degree.eq.2,1,mask_2_5_degree)
 mask_2_5_degree!0        = "lat"
 mask_2_5_degree!1        = "lon"
  
 mask_2_5_degree_reorder = mask_2_5_degree(lon|:,lat|:) 
 
;------ Read in CESM2-LE and interpolate it into 2.5 degree ------------------------ 
 ; size is 114 yr*192 lat*288 lon*90 ens
 ; ncl0*ncl1*ncl2*ncl3
 CESM2ID = addfile("Ann_tas_CESM2_LE_global_1901_2014_trend.nc","r")
 
 tas_ann_nE_temp  = CESM2ID->tas_ann_nE
 tas_ann_nE_temp2 = tas_ann_nE_temp(ncl0|:,ncl3|:,ncl1|:,ncl2|:)
 
 SID      = addfile("b.e21.BHISTcmip6.f09_g17.LE2-1301.010.cam.h0.TREFHT.199001-199912.nc","r")
 CESM_lat = SID->lat
 CESM_lon = SID->lon
 
 gwt         = latRegWgt(CESM_lat, typeof(tas_ann_nE_temp2), 0)
 opt1        = True
 opt1@critpc = 75
 tas_ann_nE_temp3 = new((/72,144,114,90/),"float")

 do yr = 0, 113
   do ens = 0, 89
     tas_ann_nE_temp3(:,:,yr,ens) = area_hi2lores_Wrap(CESM_lon, CESM_lat, rm_single_dims(tas_ann_nE_temp2(yr,ens,:,:)), True, gwt, lon_2_5, lat_2_5, opt1)
     ;tas_ann_nE_temp3(:,:,yr,ens) = area_conserve_remap_Wrap(CESM_lon, CESM_lat, rm_single_dims(tas_ann_nE_temp2(yr,ens,:,:)), lon_2_5, lat_2_5, False)
   end do
 end do

 system("/bin/rm -f Ann_tas_CESM2_LE_global_1901_2014_2_5_degree.nc")
 NID = addfile("Ann_tas_CESM2_LE_global_1901_2014_2_5_degree.nc","c")
 NID->tas_ann_nE = tas_ann_nE_temp3
 
; ******define constant value, read in prepared obs data ************************ 
; ******** Add Observation, it is global gridded data ***************************
 HadCRUT_f_ID = addfile("~/cru_ts4.03.1901.2018.tmp.dat.nc","r")
 tmp1 = (HadCRUT_f_ID->tmp)
 tmp  = tmp1
 tmp(:,:,0:359)   = tmp1(:,:,360:719)
 tmp(:,:,360:719) = tmp1(:,:,0:359)

 lat_CRU = tmp&lat
 lon_CRU = tmp&lon
 lon_CRU = where(lon_CRU.lt.0, lon_CRU+360, lon_CRU)

 ;gwt  = latRegWgt(tmp&lat, typeof(tmp), 0)
 ;opt1 = True
 ;opt1@critpc = 75
 ;HadCRUT = area_hi2lores_Wrap(lon_CRU, lat_CRU, tmp, True, gwt, lon_2_5, lat_2_5, opt1)
 HadCRUT = area_conserve_remap_Wrap(lon_CRU, lat_CRU, tmp, lon_2_5, lat_2_5, False) 

 HadCRUT_ann_all = new((/118,lsm_size(0),lsm_size(1)/),"float")
 HadCRUT_ann_all = month_to_annual(HadCRUT(:,:,:),1)
  
 NID = addfile("Ann_tas_HadCRUT_global_1901_2018_2_5_degree.nc","c")
 NID->HadCRUT_ann_all = HadCRUT_ann_all
 
 

  
  
  
