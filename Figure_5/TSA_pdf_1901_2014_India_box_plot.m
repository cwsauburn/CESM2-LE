% read in CMIP6 2.5 degree TSA and CESM-LE 1.0 degree TSA
% calculate probability distribution function
% find the ensemble member with Simth et al 2020 method (anology analysis)

clc,clear

CRUT_ave_all = zeros(1,3);
com_all_all  = zeros(3,310);
com_sel_all  = zeros(3,31);

for period = 1:4

% read in ensemble member ranked dataset ----------------------------------
 filename = ['Ann_India_CMIP6_CESM2_LE_GISS_tas_1901_1950.nc';...
 'Ann_India_CMIP6_CESM2_LE_GISS_tas_1901_1960.nc';...
 'Ann_India_CMIP6_CESM2_LE_GISS_tas_1901_1970.nc';...
 'Ann_India_CMIP6_CESM2_LE_GISS_tas_1901_1980.nc';];

% read in mask to mask global trend for continent field -------------------
% size is lon 144*lat 72 mask= 1 or 2
% India 29:35, 39:45 lon*lat

CMIP_mask_temp0 = ncread('mask_2_5_degree.nc','lsm');
CMIP_mask_temp1 = int8(or(CMIP_mask_temp0==1,CMIP_mask_temp0==2));

CMIP_mask_temp2 = CMIP_mask_temp1;
CMIP_mask_temp2(:,:)  = 0;
CMIP_mask_temp2(29:35,39:45) = CMIP_mask_temp1(29:35,39:45);

% size is lat 72*lon 144, mask= 1
%CMIP_mask     = permute(CMIP_mask_temp2,[2 1]);
       
% size is lon 288*lat 192 mask= 1 or 2
% India defination  lon{72.8:90}, lat {9.6:25.3}, 58:71, 106:123, 
CESM_mask_temp0 = ncread('mask_2_5_degree_rg2CESM2.nc','lsm');
CESM_mask_temp1 = flip(CESM_mask_temp0,2);
CESM_mask_temp2 = int8(or(CESM_mask_temp1==1,CESM_mask_temp1==2));

CESM_mask       = CESM_mask_temp2;
CESM_mask(:,:)  = 0;
CESM_mask(58:71,106:123) = CESM_mask_temp2(58:71,106:123);

% CESM_mask_temp3 = permute(CESM_mask,[2 1]);

% ensemble order ----------------------------------------------------------
CMIP6_IPO_ord_temp = ncread(filename(period,:),'CMIP6_AMV_ord');	
CMIP6_IPO_ord = flip(CMIP6_IPO_ord_temp,1); 

% read in TSA file and mask them ------------------------------------------
% 220 ensx 144lon x 72lat x 2
CMIP_ID = ['TSA_CMIP6_global_1901_2014_period',num2str(period),'.nc'];
CMIP_TSA = ncread(CMIP_ID,'TSA');

CMIP_ave = zeros(1,220)+NaN;
for ensi = 1:220
    CMIP_temp = squeeze(CMIP_TSA(ensi,:,:,2));
    CMIP_ave(ensi) = mean(CMIP_temp(CMIP_mask_temp2==1))*30.0;
end

% 90 ens x 288lon x 192lat x 2
CESM_ID = ['TSA_CESM2_LE_global_1901_2014_90ens_period',num2str(period),'.nc'];

CESM_TSA_temp1 = ncread(CESM_ID,'TSA');

CESM_ave = zeros(1,90)+NaN;
for ensi = 1:90
    CESM_temp = squeeze(CESM_TSA_temp1(ensi,:,:,2));
    CESM_ave(ensi) = mean(CESM_temp(CESM_mask==1))*30.0;
end

com_all = [CMIP_ave,CESM_ave];
com_sel = com_all(CMIP6_IPO_ord(1:31)+1);

% read in observation trend (CRU TS 4.0.3) --------------------------------
% lat,lon,record 192x288x2
CRUT_ID = ['TSA_CRUT_global_1901_2014_period',num2str(period),'.nc'];
CRUT_tsa_temp1 = ncread(CRUT_ID,'TSA');
CRUT_tsa_temp2 = squeeze(CRUT_tsa_temp1(:,:,2)); 
CRUT_ave       = nanmean(CRUT_tsa_temp2(CESM_mask==1))*30.0;

CRUT_ave_all(period) = CRUT_ave;
com_all_all(period,:) = com_all;
com_sel_all(period,:) = com_sel;

%a(period,1) = sqrt(mean((com_all-CRUT_ave).^2));
%a(period,2) = sqrt(mean((com_sel-CRUT_ave).^2));
end

%%
figure 
ylim([-0.3 0.5])
xlim([0 5])

%TTAS = 

hold on
for period = 1:4
   boxplot(com_all_all(period,:)','positions',[period-0.2],'Colors',[150/225, 75/225, 0], 'Symbol','')
   %A1 = quantile(com_all_all(period,:)',[0.05 0.95]);
   %plot([period-0.2],A1(1),'s','MarkerSize',9,'MarkerEdgeColor',[150/225, 75/225, 0],'MarkerFaceColor','w')
   %plot([period-0.2],A1(2),'s','MarkerSize',9,'MarkerEdgeColor',[150/225, 75/225, 0],'MarkerFaceColor','w')

   boxplot(com_sel_all(period,:)','positions',[period+0.2],'Colors','b', 'Symbol','')
   %B1 = quantile(com_sel_all(period,:)',[0.05 0.95]);
   %plot([period+0.2],B1(1),'s','MarkerSize',9,'MarkerEdgeColor','b','MarkerFaceColor','w')
   %plot([period+0.2],B1(2),'s','MarkerSize',9,'MarkerEdgeColor','b','MarkerFaceColor','w')
end
plot([1,2,3,4],CRUT_ave_all,'ko', 'markersize',7)

ylabel('Temperature Trend (°C/30yr)')
%title(tit_name(kind))

set(gca,'FontSize',12)
set(gca,'fontname','Times') 

n = findobj(gca,'tag','Outliers');
%n.MarkerEdgeColor = rgb('white');

xticks([1 2 3 4])
xticklabels({'1951-1980','1961-1990','1971-2000','1981-2010'})




