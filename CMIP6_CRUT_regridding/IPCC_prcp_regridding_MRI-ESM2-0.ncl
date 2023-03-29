
;IPCC_tas_regridding_MRI-ESM2-0.ncl

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

 modName = "MRI-ESM2-0"
 
 modInName = "MRI-ESM2-0"
 
 varName0 = "PRCP"
 
 varName1 = "pr_Amon"
 
 varName2 = "pr"
 
 varTypName = "ATM"

 interval = "month"
  
 nEns = 5
 
;********Now Read the data*****************************************

ucf = 60.0*60.0*24.0*365.0*(1.0/12.0)


do i = 0, nEns-1
  
   EnsName = "r"+tostring(i+1)+"i1p1f1"
   
   fils = systemfunc("ls "+inDir+"/"+modName+"/"+varTypName+"/"+interval+"/"+varName0+"/"+varName1+"_"+modName+"_historical_"+EnsName+"_"+"*.nc")  
    
   f = addfiles(fils, "r")
   print(fils)

   ListSetType(f, "cat")
   var1 = f[:]->$varName2$
   var = var1(:, :, :)
   delete(var1)

   var = var*ucf
   var@units = "mm/month"

   printVarSummary(var) 

   var_rg = area_conserve_remap_Wrap(var&lon, var&lat, var, lon_2_5, lat_2_5, False)

   system("/bin/rm -f "+outDir+"/"+modName+"/"+varTypName+"/"+interval+"/"+varName0+"/"+varName1+"_"+modName+"_historical_"+EnsName+"_185001-201412_rg.nc")
       
   fout1 = addfile(outDir+"/"+modName+"/"+varTypName+"/"+interval+"/"+varName0+"/"+varName1+"_"+modName+"_historical_"+EnsName+"_185001-201412_rg.nc", "c") 
 
   fout1->pr = var_rg
   delete(var)
   delete(fils)
   delete(f)
   delete(EnsName)
   
end do

;**************************************************************




