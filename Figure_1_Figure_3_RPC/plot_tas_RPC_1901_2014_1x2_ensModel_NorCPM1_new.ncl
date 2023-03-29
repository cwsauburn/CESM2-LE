
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  

 ;------ Read mask for 85N-60S, and cover continent -----------------------
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
 mask_2_5_degree!0      = "lat"
 mask_2_5_degree!1      = "lon"

 ;------- read in all data together, lon, lat, 144*72 ----------
 RPC_ID = addfile("RPC_global_continent_1951_2014_NorCPM1.nc","r")
 
 ; long 144*lat 72*record 1
 Trend_Obs_ID = addfile("~/tas_Ann_GISS_1901-2014_trend_hurst.nc","r") 
 TSA_GISS     = rm_single_dims(Trend_Obs_ID->TSA)
 TSA_GISS@_FillValue = -9999 
 
 RPC_Had  = RPC_ID ->RPC_Had  
 RPC_GISS = RPC_ID ->RPC_GISS

 RPC_GISS = mask(RPC_GISS, .not.ismissing(TSA_GISS), True)

 RPC_Had!0      = "lon"
 RPC_Had!1      = "lat"
 RPC_Had&lon    = lsm_temp&lon 
 RPC_Had&lat    = lsm_temp&lat  

 RPC_GISS!0      = "lon"
 RPC_GISS!1      = "lat"
 RPC_GISS&lon    = lsm_temp&lon 
 RPC_GISS&lat    = lsm_temp&lat  
  
 sample = new((/lsm_size(0),lsm_size(1)/),"float",-9999)
 sample!0   = "lat"
 sample&lat = lsm_temp&lat
 sample!1   = "lon"
 sample&lon = lsm_temp&lon 
 
 RPC_Had_results  = sample 
 RPC_GISS_results = sample
 
 RPC_Had_results  = RPC_Had(lat|:,lon|:)
 RPC_GISS_results = RPC_GISS(lat|:,lon|:)

 ;-------- Read in significant area ---------------------
 RPC_IDH    = addfile("RPC_global_continent_1951_2014_1000_sample_CRUT_NorCPM1.nc","r")
 ; long 144*lat 72*record 120
 Siglv_Had  = RPC_IDH->Siglv_Had
 
 RPC_IDG    = addfile("RPC_global_continent_1901_2014_1000_sample_GISS.nc","r")
 Siglv_GISS = RPC_IDG->Siglv_GISS
 
 Had_sig_area = new(dimsizes(Siglv_Had),"integer")
 Had_sig_area = where(Siglv_Had.le.0.025, 1, 0)
 Had_sig_area = where(ismissing(Had_sig_area), 0, Had_sig_area)
 Had_sig_area = where(ismissing(RPC_Had), 0, Had_sig_area)

 GISS_sig_area = new(dimsizes(Siglv_GISS),"integer")
 GISS_sig_area = where(Siglv_GISS.le.0.025, 1, 0)
 GISS_sig_area = where(ismissing(GISS_sig_area), 0, GISS_sig_area)
 GISS_sig_area = where(ismissing(RPC_GISS), 0, GISS_sig_area)
  
 Had_sig_area!0   = "lon" 
 Had_sig_area!1   = "lat" 

 Had_sig_area&lon = lon_2_5
 Had_sig_area&lat = lat_2_5
 
 Had_Sarea  = sample 
 Had_Sarea  = Had_sig_area(lat|:,lon|:)
   
 GISS_sig_area!0   = "lon" 
 GISS_sig_area!1   = "lat" 

 GISS_sig_area&lon = lon_2_5
 GISS_sig_area&lat = lat_2_5
 
 GISS_Sarea  = sample 
 GISS_Sarea  = GISS_sig_area(lat|:,lon|:)
 
 ;-------- Read in AREA MASK, and calculate total area -------------
 ; [lat | 72] x [lon | 144]
 farea = addfile("area_2_5_degree.nc", "r")
 carea = farea->CAREA
 carea@_FillValue = -9999      ; no fill value in original data, only used for masked data
 carea      = mask(carea, mask_2_5_degree.eq.1, True)
 tot_area   = sum(carea) 
   
 STP_LTP_Had = sum(mask(carea, Had_Sarea,1))
 print(STP_LTP_Had/tot_area)

 STP_LTP_GISS = sum(mask(carea, GISS_Sarea,1))
 print(STP_LTP_GISS/tot_area)
 print(tot_area) 

 ;-------- Plot the three panels as figure 3-----------------------------
 ;******** displying the time mean data ******************************
 wks = gsn_open_wks("png", "TAS_RPC_1901_2014_spatial_global_conus_NorCPM1")

 gsn_define_colormap(wks, "BlueDarkRed18") 
 res = True

 res@gsnDraw            = False             ; don't draw yet
 res@gsnFrame           = False             ; don't advance frame yet
 res@lbLabelsOn         = False
 res@cnInfoLabelOn      = False
 res@gsnMaximize        = True

 res@mpLimitMode   = "LatLon"
 res@mpOutlineBoundarySets  = "Geophysical"
 res@mpGeophysicalLineColor = "black"
 res@mpUSStateLineColor     = "black"

 res@pmTickMarkDisplayMode  = "Always"

 res@cnFillOn   = True
 res@cnLinesOn  = False
 res@cnFillMode = "RasterFill"
 res@cnLineLabelsOn = False

 res@gsnSpreadColors = True

 res@lbOrientation = "Vertical"

 res@lbBoxLinesOn         = False           ; Turn off labelbar box lines
 res@lbLabelAutoStride    = True            ; Control labelbar labels.

 res@lbLabelBarOn         = False 
 res@cnLevelSelectionMode = "ManualLevels"

 res@mpMinLatF    = -60.0
 res@mpMaxLatF    = 85.0
 res@mpProjection = "Robinson"
 res@gsnMaximize  = True
 res@cnMissingValFillColor = "gray"

 res@tmXBOn = False
 res@tmYLOn = False
 res@tmXTOn = False
 res@tmYROn = False

 ; *** the decide which maximum and minimum value you want to assign *****
 ; **** accordingly decide the min, max and stide for manual cn level selction *****
 res@cnLevelSelectionMode = "ExplicitLevels"
 res@cnLevels             = (/-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2/)

 ;******************************************

 plots = new(1, graphic)
 plots(0) = gsn_csm_contour_map(wks, RPC_Had_results, res)
 
 ;----------- plot significant area -------------------
 res2 = True
 res2@gsnDraw            = False             ; don't draw yet
 res2@gsnFrame           = False             ; don't advance fra
 
 res2@gsnAddCyclic = True
 res2@gsnMaximize = True
 
 res2@cnInfoLabelOn       = False       ; turn off info label
 res2@cnLinesOn           = False       ; do not draw contour lines
 res2@cnLineLabelsOn      = False       ; do not draw contour labels
 res2@cnConstFEnableFill = True         ; Fill even scalar field is constant
 res2@cnFillScaleF        = 0.5         ; add extra density
 ;res2@cnFillMode = "RasterFill"
 
 res2@cnLevelSelectionMode = "ManualLevels"
 res2@cnMinLevelValF  =  1.2
 res2@cnMaxLevelValF  =  0.8
 res2@cnLevelSpacingF =  0.1
 
 opt = True
 opt@gsnShadeFillType = "pattern"       ; pattern fill
 opt@gsnShadeHigh = 3                   ; use pattern #2
 opt@gsnShadeLow  = 3
 res2@cnConstFEnableFill = True        ; Fill even scalar field is constant 
 ;------------------------------------------------------- 

 plots2 = new(1, graphic)

 plots2(0) = gsn_csm_contour(wks, Had_Sarea, res2)
 plots2(0) = gsn_contour_shade(plots2(0), -0.2, 0.8, opt)  ; shade all areas greater than 95 CI

 overlay(plots(0), plots2(0))
 
 panres1                   = True                   ; modify the panel plot
 panres1@gsnFrame          = False                  ; don't advance frame yet
 panres1@gsnPanelLabelBar  = True                   ; add common colorbar
 panres1@lbOrientation     = "horizontal"
 panres1@gsnMaximize    = True
 panres1@gsnFrame       = False

 txres = True
 txres@txFontHeightF = 0.050

 panres1@lbOrientation     = "horizontal"
  
 panres1@gsnPanelMainString      = "TAS RPC 1901-2014"
 panres1@gsnPanelMainFontHeightF = 0.015

 panres1@amJust = "BottomLeft"
 panres1@amOrthogonalPosF = 0.48
 panres1@amParallelPosF   = -0.1
 panres1@gsnPanelFigureStringsPerimOn = False
 panres1@gsnPanelFigureStringsBackgroundFillColor = -1
 panres1@gsnPanelFigureStringsFontHeightF = 0.012
  

 gsn_panel(wks, plots, (/1, 1/), panres1)  
 frame(wks) 
 
 
