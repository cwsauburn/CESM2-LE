
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  

 ;------ Read mask for 85N-60S, and cover continent -----------------------
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
 mask_2_5_degree!0        = "lat"
 mask_2_5_degree!1        = "lon"

 mask_2_5_degree_reorder = mask_2_5_degree(lon|:,lat|:)

 ;------- read in all data together, lon, lat, 144*72 ----------
 RPC_ID  = addfile("Ann_tas_CESM2_LE_2_5_degree_RPC_global_continent_1901_2014.nc","r")
 RPC_Had = RPC_ID ->RPC_Had  
 
 RPC_Had!0      = "lon"
 RPC_Had!1      = "lat"
 RPC_Had&lon    = lsm_temp&lon 
 RPC_Had&lat    = lsm_temp&lat  
 
 sample     = new((/lsm_size(0),lsm_size(1)/),"float",-9999)
 sample!0   = "lat"
 sample&lat = lsm_temp&lat
 sample!1   = "lon"
 sample&lon = lsm_temp&lon 
 
 RPC_Had_results  = sample 
 RPC_Had_results  = RPC_Had(lat|:,lon|:)

 ;-------- Read in significant area ---------------------
 RPC_IDH    = addfile("CESM2_LE_2_5_degree_RPC_global_continent_1901_2014_1000_sample_CRUT.nc","r")
 ; long 144*lat 72*record 120
 Siglv_Had  = RPC_IDH->Siglv_Had
  
 Had_sig_area = new(dimsizes(Siglv_Had),"integer")
 Had_sig_area = where(Siglv_Had.le.0.025, 1, 0)
 Had_sig_area = where(ismissing(Had_sig_area), 0, Had_sig_area)
 Had_sig_area = where(ismissing(RPC_Had), 0, Had_sig_area)

 Had_sig_area!0   = "lon" 
 Had_sig_area!1   = "lat" 

 Had_sig_area&lon = lon_2_5
 Had_sig_area&lat = lat_2_5
 
 Had_Sarea  = sample 
 Had_Sarea  = Had_sig_area(lat|:,lon|:)
    
 ;-------- Plot the three panels as figure 3-----------------------------
 ;******** displying the time mean data ******************************
 wks = gsn_open_wks("png", "CESM2_LE_2_5_degree_TAS_RPC_1901_2014_spatial_global_conus")

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
 plots = new(2, graphic)
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
 ;opt@gsnShadeFillScalesF = 0.5
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
  
 panres1@gsnPanelMainString      = "CESM2-LE 2_5_degree TAS RPC 1901-2014 sigratio 0.39"
 panres1@gsnPanelMainFontHeightF = 0.015

 panres1@amJust = "BottomLeft"
 panres1@amOrthogonalPosF = 0.48
 panres1@amParallelPosF   = -0.1
 panres1@gsnPanelFigureStringsPerimOn = False
 panres1@gsnPanelFigureStringsBackgroundFillColor = -1
 panres1@gsnPanelFigureStringsFontHeightF = 0.012
  

  gsn_panel(wks, plots, (/1, 1/), panres1)  
  frame(wks) 
 
 
