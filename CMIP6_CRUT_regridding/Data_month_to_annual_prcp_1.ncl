
; Read in CMIP6, month to annual, two time periods are 1901-2014, 1951-2014

; *********Add CMIP6 14 + 12 model temperature 1901-2014 ************************************************ 
; base period: 1951-1980 degrees Celsius (deg-C).
; dimension  time = 1980 (185001-201412, 1980 months), lat = 72 ,lon = 144
  
  time1 = 1901
  time2 = 1951
  time3 = 2014
  lat   = 72
  lon   = 144
  modName  = (/"BCC-CSM2-MR","BCC-ESM1","CESM2","CESM2-WACCM","CNRM-CM6-1","CNRM-ESM2-1","CanESM5","EC-Earth3","GISS-E2-1-G","GISS-E2-1-H","IPSL-CM6A-LR","MIROC6","MRI-ESM2-0","UKESM1-0-LL"/) 
  
  varName0 = "PRCP"
  varName1 = "pr_Amon"
  varName2 = "pr"
  varTypName = "ATM" 
  interval = "month" 
   
  inDir = "~/"  
  
  nMod = 14 
  nEns = (/3, 3, 9, 3, 10, 5, 25, 2, 10, 10, 31, 10, 5, 5/) ;sum=131
    
  norL = 1980
  norY = norL/12
   
  CMIP6_1 =  new((/time3-time1+1,lat,lon,sum(nEns)/),"float")
  CMIP6_2 =  new((/time3-time2+1,lat,lon,sum(nEns)/),"float")
  CMIP6_1!0  = "year"  
  CMIP6_1!1  = "lat"
  CMIP6_1!2  = "lon"
  CMIP6_1!3  = "ensemble"
  CMIP6_2!0  = "year"  
  CMIP6_2!1  = "lat"
  CMIP6_2!2  = "lon"
  CMIP6_2!3  = "ensemble"  
  GISS_f_ID  = addfile("GISS_global_2.5.nc","r")
  CMIP6_1&lat= GISS_f_ID->lat
  CMIP6_1&lon= GISS_f_ID->lon
  CMIP6_2&lat= GISS_f_ID->lat
  CMIP6_2&lon= GISS_f_ID->lon

  count = -1
  do nM6 = 0, nMod-1
      
     fils = systemfunc("ls "+inDir+"/"+modName(nM6)+"/"+varTypName+"/"+interval+"/"+varName0+"/"+varName1+"_"+modName(nM6)+"_historical_"+"r*_"+"*.nc")  
         
     do nE6 = 0, nEns(nM6)-1
        count = count+1
        
        flist = addfiles(fils(nE6), "r")
        namestr = tochar(fils(nE6))
        listsize = dimsizes(namestr)
        betime1 = (time1-stringtointeger(namestr(listsize-19:listsize-16)))*12-(stringtointeger(namestr(listsize-15:listsize-14))-1)  
        betime2 = (time2-stringtointeger(namestr(listsize-19:listsize-16)))*12-(stringtointeger(namestr(listsize-15:listsize-14))-1)
        edtime1 = (time3-time1+1)*12
        edtime2 = (time3-time2+1)*12
        
        delete(namestr)
        delete(listsize)
        
        pure_data_temp = flist[:]->$varName2$
        pure_data1 = pure_data_temp(betime1:betime1+edtime1-1,:,:)
        pure_data2 = pure_data_temp(betime2:betime2+edtime2-1,:,:)
        delete(betime1)
        delete(betime2)
        delete(pure_data_temp)
     
        D1 = dimsizes(pure_data1)
        D2 = dimsizes(pure_data2)
        
        CMIP6_1(:(D1(0)/12)-1,:,:,count) = month_to_annual(pure_data1(0:(D1(0)/12)*12-1,:,:),1) 
        CMIP6_2(:(D2(0)/12)-1,:,:,count) = month_to_annual(pure_data2(0:(D2(0)/12)*12-1,:,:),1) 
        
        delete(D1)
        delete(D2)
        delete(flist)
        delete(pure_data1)
        delete(pure_data2)
     end do
     delete(fils)
  end do

  system("/bin/rm -f Ann_prcp_CMIP6_1_1901_2014.nc")
  system("/bin/rm -f Ann_prcp_CMIP6_1_1951_2014.nc")  
          
  CMIP6_01_reorder = CMIP6_1(lon|:,lat|:,year|:,ensemble|:)   
  CMIP6_02_reorder = CMIP6_2(lon|:,lat|:,year|:,ensemble|:) 
  
  CMIP6_ID_01 = addfile("Ann_prcp_CMIP6_1_1901_2014.nc","c")
  CMIP6_ID_01->CMIP6 = CMIP6_01_reorder 
  CMIP6_ID_02 = addfile("Ann_prcp_CMIP6_1_1951_2014.nc","c")
  CMIP6_ID_02->CMIP6 = CMIP6_02_reorder 
  
    

   