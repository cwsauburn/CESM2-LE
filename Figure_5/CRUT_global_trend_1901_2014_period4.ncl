 ; read in CESM2-LE in NCAR folder, calculate 
 ; Mann-Kendall non-parametric test for monotonic trend
 ; (time, lat, lon) (,192,288)
  
;------ list of files & read TAS -------------
    
  diri   = "~/"
  
  ; size 288lon*192lat*114year
  fileID = addfile(diri+"Ann_tas_CRUT_1901_2014_1_degree.nc","r")
  TREFHT_temp = fileID->tas
  
  TREFHT = TREFHT_temp(:,:,79:109)
 
;------ calculate trend ----------------------     
  TSA    = trend_manken(TREFHT, False, 2)
	 	 
;------ write results ------------------------
 
  filename = "TSA_CRUT_global_1901_2014_period4.nc"
  system("/bin/rm -f "+filename)
 
  outID    = addfile(filename,"c")
  outID->TSA = TSA



