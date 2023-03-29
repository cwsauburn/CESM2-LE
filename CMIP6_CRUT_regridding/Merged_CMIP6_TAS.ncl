;****Edit this section to get the data ********************
;**********************************************************
 inDir  = "~/"
 outDir = "~/"

 varTypName = "ATM"
 interval = "month"
 varName0 = "TAS"
 
 varName1 = "tas_Amon"
 varName2 = "tas" 
 
 mod_name = (/"ACCESS-CM2","ACCESS-ESM1-5","AWI-CM-1-1-MR","BCC-CSM2-MR","CAMS-CSM1-0","CanESM5","CanESM5-CanOE","CESM2","CESM2-WACCM","CIESM","CMCC-CM2-SR5","CNRM-CM6-1","CNRM-CM6-1-HR","CNRM-ESM2-1","EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L","FGOALS-g3","FIO-ESM-2-0","GFDL-ESM4","GISS-E2-1-G","HadGEM3-GC31-LL","HadGEM3-GC31-MM","IITM-ESM","INM-CM4-8","INM-CM5-0","IPSL-CM6A-LR","KACE-1-0-G","KIOST-ESM","MCM-UA-1-0","MIROC6","MIROC-ES2L","MPI-ESM1-2-HR","MPI-ESM1-2-LR","MRI-ESM2-0","NESM3","NorESM2-LM","NorESM2-MM","TaiESM1","UKESM1-0-LL"/) 
 mod_ensem = new((/40/),"integer")


;********Now Read the data*****************************************

;ucf = 60.0*60.0*24.0*365.0*(1.0/12.0)
mod_num = 7
modName  = mod_name(mod_num)


do i = 7,11

   EnsName = "r"+tostring(i)+"i1p1f1"
   
   ; make sure we are dealing with one ensemble member each time
   fils = systemfunc("ls "+inDir+"/"+modName+"/"+varTypName+"/"+interval+"/"+varName0+"/"+varName0+"/"+"*"+EnsName+"*.nc")    
   
   f = addfiles(fils, "r")
   ListSetType(f, "cat")
   var = f[:]->$varName2$   
 
   var = var - 273.15
   var@units = "degree C"
   
   ; connect all files together  
   system("/bin/rm -f "+inDir+"/"+modName+"/"+varTypName+"/"+interval+"/"+varName0+"/"+varName0+"/"+"tas_Amon_historical_"+EnsName+"_gn_185001-201412.nc")   
   fout1 = addfile(inDir+"/"+modName+"/"+varTypName+"/"+interval+"/"+varName0+"/"+varName0+"/"+"tas_Amon_historical_"+EnsName+"_gn_185001-201412.nc","c")  	  
   fout1->$varName2$ = var   
   
   delete(fils)    
   delete(f)
   delete(var)
   delete(fout1)
end do
