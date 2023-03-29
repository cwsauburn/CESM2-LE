clc,clear

filename = 'Ann_India_CMIP6_CESM2_LE_GISS_tas_1901_2014.nc';

sel_year = 5:109;

CMIP6_IPO = ncread(filename,'CMIP6_AMV');
CMIP6_IPO = CMIP6_IPO(:,sel_year);

cli_mod = nanmean(CMIP6_IPO,1);

GISS_IPO  = ncread(filename,'GISS_AMV');
GISS_IPO  = GISS_IPO(sel_year);

COEALL1   = zeros(2,310);
for i = 1:310
 COEALL1(:,i) = polyfit(sel_year,CMIP6_IPO(i,:),1)*114;
end
COEALL2 = polyfit(sel_year,cli_mod,1)*114;
pd2     = fitdist(COEALL1(1,:)','Normal');
range2  = 1.96*0.4868/sqrt(310);

COEALL3 = polyfit(sel_year,GISS_IPO,1)*114;

%--------------------------------------------------------
CMIP6_IPO_ord_temp = ncread(filename,'CMIP6_AMV_ord');	
CMIP6_IPO_ord = flip(CMIP6_IPO_ord_temp,1); 

CMIP6_mean1 = nanmean(CMIP6_IPO(CMIP6_IPO_ord(1:31)+1,:),1);

COEALL4 = COEALL1(:,CMIP6_IPO_ord(1:31)+1);
pd4     = fitdist(COEALL4(1,:)','Normal');
range4  = 1.96*0.165/sqrt(31); 
COEALL5 = polyfit(sel_year,CMIP6_mean1,1)*114;




