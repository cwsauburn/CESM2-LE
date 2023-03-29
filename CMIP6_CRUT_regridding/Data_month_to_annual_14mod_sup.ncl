; Read in supplement old 14 model added ensemble members, and calculate month to annual, two time period is 1901-2014, 1951-2014
; Please find Table excel file to find the added ensemble members  

; Read in CMIP6, month to annual, two time period is 1901-2015, 1951-2015

; *********Add CMIP6 5 model temperature************************************************ 
; base period: 1951-1980 degrees Celsius (deg-C).
; dimension  time = 1980 (185001-201412, 1980 months), lat = 72 ,lon = 144
  
  time1 = 1901
  time2 = 1951
  time3 = 2014
  lat   = 72
  lon   = 144
  modName  = (/"CESM2","CNRM-CM6-1","CanESM5","IPSL-CM6A-LR","UKESM1-0-LL"/) 
  expName = (/(/"r11i1p1f1","*********","*********","*********","*********","*********","*********"/),(/"r11i1p1f2","r12i1p1f2","r13i1p1f2","r14i1p1f2","r15i1p1f2","r16i1p1f2","r17i1p1f2"/),(/"*********","*********","*********","*********","*********","*********","*********"/),(/"r32i1p1f1","*********","*********","*********","*********","*********","*********"/),(/"r5i1p1f3*","r6i1p1f3*","r7i1p1f3*","r9i1p1f2*","*********","*********","*********"/)/)
  
  varName0 = "PRCP"
  varName1 = "pr_Amon"
  varName2 = "pr"
  varTypName = "ATM" 
  interval = "month" 
     
  inDir = "~"  
  
  nMod = 5 
  nEns = (/1,7,25,1,4/) ;sum=38
    
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
        
     do nE6 = 0, nEns(nM6)-1
        count = count+1
        
        if (nM6.ne.2) then       
           fils    = systemfunc("ls "+inDir+"/"+modName(nM6)+"/"+varTypName+"/"+interval+"/"+varName0+"/"+varName1+"_"+modName(nM6)+"_historical_"+expName(nM6,nE6)+"*.nc") 
           flist   = addfiles(fils, "r")
           namestr = tochar(fils)
        else
           fils    = systemfunc("ls "+inDir+"/"+modName(nM6)+"/"+varTypName+"/"+interval+"/"+varName0+"/"+varName1+"_"+modName(nM6)+"_historical_"+"r*i1p2f1_"+"*.nc") 
           flist   = addfiles(fils(nE6), "r")
           namestr = tochar(fils(nE6))             
        end if

        listsize = dimsizes(namestr)
        betime1 = (time1-stringtointeger(namestr(listsize-19:listsize-16)))*12-(stringtointeger(namestr(listsize-15:listsize-14))-1)  
        betime2 = (time2-stringtointeger(namestr(listsize-19:listsize-16)))*12-(stringtointeger(namestr(listsize-15:listsize-14))-1)
        delete(namestr)
        delete(listsize)
        
        pure_data_temp = flist[:]->$varName2$
        pure_data1 = pure_data_temp(betime1:,:,:)
        pure_data2 = pure_data_temp(betime2:,:,:)
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
        delete(fils)
     end do 
  end do
  
  system("/bin/rm -f Ann_prcp_CMIP6_3_1901_2014.nc")
  system("/bin/rm -f Ann_prcp_CMIP6_3_1951_2014.nc")  
          
  CMIP6_01_reorder = CMIP6_1(lon|:,lat|:,year|:,ensemble|:)   
  CMIP6_02_reorder = CMIP6_2(lon|:,lat|:,year|:,ensemble|:) 
  CMIP6_ID_01 = addfile("Ann_prcp_CMIP6_3_1901_2014.nc","c")
  CMIP6_ID_01->CMIP6 = CMIP6_01_reorder 
  CMIP6_ID_02 = addfile("Ann_prcp_CMIP6_3_1951_2014.nc","c")
  CMIP6_ID_02->CMIP6 = CMIP6_02_reorder 
  
    
  