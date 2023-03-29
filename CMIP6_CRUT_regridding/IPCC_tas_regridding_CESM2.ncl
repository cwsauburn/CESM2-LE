
;IPCC_tas_regridding_CESM2.ncl

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

;************Now let us regrid and then mask the data**************************************
;****************Adding land mask file***************************

 flsm = addfile("~/mask_2_5_degree.nc", "r")

 lsm = flsm->lsm

 lat_2_5 = lsm&lat
 lon_2_5 = lsm&lon

 lat_2_5!0 = "lat"
 lat_2_5&lat = lsm&lat

;****Edit this section to get the data ********************
;****

 inDir = "~/"
 outDir = "~/"

 modName = "CESM2"
 
 modInName = "CESM2"
 
 varName0 = "TAS"
 
 varName1 = "tas_Amon"
 
 varName2 = "tas"
 
 varTypName = "ATM"

 interval = "month"
  
 nEns = 11

;********Now Read the data*****************************************

do i = 10, 10

   EnsName = "r"+tostring(i+1)+"i1p1f1"
   
   fils = systemfunc("ls "+inDir+"/"+modName+"/"+varTypName+"/"+interval+"/"+varName0+"/"+varName0+"/"+varName1+"_"+modName+"_historical_"+EnsName+"_"+"*.nc")
     
   f = addfiles(fils, "r")
   print(fils)

   ListSetType(f, "cat")
   var = f[:]->$varName2$
   
   var = var - 273.15
   var@units = "degree C"

   printVarSummary(var) 

   var_rg = area_conserve_remap_Wrap(var&lon, var&lat, var, lon_2_5, lat_2_5, False)

   system("/bin/rm -f "+outDir+"/"+modName+"/"+varTypName+"/"+interval+"/"+varName0+"/"+varName0+"/"+varName1+"_"+modName+"_historical_"+EnsName+"_185001-201412_rg.nc")
   
   fout1 = addfile(outDir+"/"+modName+"/"+varTypName+"/"+interval+"/"+varName0+"/"+varName0+"/"+varName1+"_"+modName+"_historical_"+EnsName+"_185001-201412_rg.nc", "c")    
   
   fout1->tas = var_rg
   delete(var)
   delete(fils)
   delete(f)
   delete(EnsName)
   
end do

;**************************************************************



