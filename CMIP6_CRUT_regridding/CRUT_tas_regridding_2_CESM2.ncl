
 ;IPCC_tas_regridding_All.ncl

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;load "./taylor_diagram.ncl"

 inDir  = "~/"
 outDir = "~/"

 varTypName = "ATM"
 interval = "month"
 varName0 = "TAS"
 
 varName1 = "tas_Amon"
 varName2 = "tas" 
 
 mod_name = (/"ACCESS-CM2","ACCESS-ESM1-5","AWI-CM-1-1-MR","BCC-CSM2-MR","CAMS-CSM1-0","CanESM5","CanESM5-CanOE","CESM2","CESM2-WACCM","CIESM","CMCC-CM2-SR5","CNRM-CM6-1","CNRM-CM6-1-HR","CNRM-ESM2-1","EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L","FGOALS-g3","FIO-ESM-2-0","GFDL-ESM4","GISS-E2-1-G","HadGEM3-GC31-LL","HadGEM3-GC31-MM","IITM-ESM","INM-CM4-8","INM-CM5-0","IPSL-CM6A-LR","KACE-1-0-G","KIOST-ESM","MCM-UA-1-0","MIROC6","MIROC-ES2L","MPI-ESM1-2-HR","MPI-ESM1-2-LR","MRI-ESM2-0","NESM3","NorESM2-LM","NorESM2-MM","TaiESM1","UKESM1-0-LL"/) 
 mod_ensem = new((/40/),"integer")
 
;*************Now let us regrid and then mask the data***********
;****************Adding land mask file***************************
 ; Add cru_ts4.03.1901.2018.tmp.dat.nc 0.5*0.5 regridded data
 ; short tempanomaly time 1416 x lat 360 x lon 720
 ; Time period(190101-201812)
 CRU_ID = addfile("~/cru_ts4.03.1901.2018.tmp.dat.nc", "r")

 tmp1 = CRU_ID->tmp
 tmp  = tmp1
 tmp(:,:,0:359)   = tmp1(:,:,360:719)
 tmp(:,:,360:719) = tmp1(:,:,0:359)
 pre_size  = dimsizes(tmp)
  
 lat_CRU = tmp&lat
 lon_CRU = tmp&lon
 lon_CRU = where(lon_CRU.lt.0, lon_CRU+360, lon_CRU)
 
 CMIP_ID = addfile("~/tas_Amon_CESM2_historical_r10i1p1f1_gn_185001-189912.nc", "r")
 var     = CMIP_ID->tas
 
 ;var_rg = area_conserve_remap_Wrap(lon_CRU,lat_CRU,tmp,var&lon,var&lat,False)
  
 gwt  = latRegWgt(tmp&lat, typeof(tmp), 0)
 opt1 = True
 opt1@critpc = 75
 var_rg = area_hi2lores_Wrap(lon_CRU,lat_CRU,tmp, True, gwt, var&lon,var&lat, opt1)
  
 system("/bin/rm -f cru_ts4.03.1901.2018.tmp.dat_rg.nc")
 fout1 = addfile("cru_ts4.03.1901.2018.tmp.dat_rg.nc","c")  	  
 fout1->$varName2$ = var_rg
	  
;---- read in regridded CRUT and calculate annual TAS ------------------
  
  CRUT_f_ID  = addfile("cru_ts4.03.1901.2018.tmp.dat_rg.nc","r")
  CRUT_tas   = CRUT_f_ID->$varName2$
 
  CRUT_ann_1 = month_to_annual(CRUT_tas((1901-1901)*12:(2014-1901+1)*12-1,:,:),1) 
  CRUT_ann_2 = month_to_annual(CRUT_tas((1951-1901)*12:(2014-1901+1)*12-1,:,:),1)
  
  CRUT_ann_reorder_1 = CRUT_ann_1(lon|:,lat|:,year|:)
  CRUT_ann_reorder_2 = CRUT_ann_2(lon|:,lat|:,year|:)  
  
  system("/bin/rm -f "+"~/Ann_tas_CRUT_1901_2014_1_degree.nc")
  fout1 = addfile("~/Ann_tas_CRUT_1901_2014_1_degree.nc","c")  	  
  fout1->$varName2$ = CRUT_ann_reorder_1

  system("/bin/rm -f "+"~/Ann_tas_CRUT_1951_2014_1_degree.nc")
  fout1 = addfile("~/Ann_tas_CRUT_1951_2014_1_degree.nc","c")  	  
  fout1->$varName2$ = CRUT_ann_reorder_2
  
;---- read in CESM and calculate annual TAS ------------------  

  sampeID = addfile(inDir+"/"+modName+"/"+varTypName+"/"+interval+"/"+varName0+"/"+varName0+"/"+"tas_Amon_CESM2_historical_r11i1p1f1_gn_185001-201412.nc","r")
  tas_11  = sampeID->$varName2$
  tas_sz  = dimsizes(tas_11)
  
  CESM_1 = new((/2014-1901+1,tas_sz(1),tas_sz(2),11/),"float") 
  CESM_2 = new((/2014-1951+1,tas_sz(1),tas_sz(2),11/),"float")

  CESM_1!0   = "year"
  CESM_1!1   = "lat"
  CESM_1!2   = "lon"
  CESM_1!3   = "ens"
	  
  CESM_2!0   = "year"
  CESM_2!1   = "lat"
  CESM_2!2   = "lon"
  CESM_2!3   = "ens"	 
   
  betime1 = (1901-1850)*12
  betime2 = (1951-1850)*12
  
 do ensi = 1,11

   EnsName = "r"+tostring(ensi)+"i1p1f1"
   
   ; make sure we are dealing with one ensemble member each time
   fils = systemfunc("ls "+inDir+"/"+modName+"/"+varTypName+"/"+interval+"/"+varName0+"/"+varName0+"/"+"*"+EnsName+"*_185001-201412.nc")    
   
   file_ID = addfile(fils, "r")
   var = file_ID->$varName2$ 
   
   if (ensi.gt.6) then
      var = var - 273.15
      var@units = "degree C"
   end if
   
   CESM_1_temp = month_to_annual(var(betime1:tas_sz(0)-1,:,:),1) 
   CESM_2_temp = month_to_annual(var(betime2:tas_sz(0)-1,:,:),1)
 
   CESM_1(:,:,:,ensi-1) = CESM_1_temp
   CESM_2(:,:,:,ensi-1) = CESM_2_temp
   
   delete(EnsName)  
   delete(fils)    
   delete(file_ID)
   delete(var)
 end do

 
 CESM_1 = month_to_annual(tas_all(betime1:tas_sz(0)-1,:,:,:),1) 
 CESM_2 = month_to_annual(tas_all(betime2:tas_sz(0)-1,:,:,:),1)
         
 CMIP6_01_reorder = CESM_1(lon|:,lat|:,year|:,ensemble|:)   
 CMIP6_02_reorder = CESM_2(lon|:,lat|:,year|:,ensemble|:) 
 
 system("/bin/rm -f Ann_tas_CMIP6_CESM2_1901_2014.nc")
 system("/bin/rm -f Ann_tas_CMIP6_CESM2_1951_2014.nc")
 
 CMIP6_ID_01 = addfile("~/Ann_tas_CMIP6_CESM2_1901_2014.nc","c")
 CMIP6_ID_01->CMIP6 = CMIP6_01_reorder 
 
 CMIP6_ID_02 = addfile("~/Ann_tas_CMIP6_CESM2_1951_2014.nc","c")
 CMIP6_ID_02->CMIP6 = CMIP6_02_reorder 
  
  
  
  
  
  
  