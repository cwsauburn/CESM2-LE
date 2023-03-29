% read in CMIP6 2.5 degree TSA and CESM-LE 1.0 degree TSA
% calculate probability distribution function
% find the ensemble member with Simth et al 2020 method (anology analysis)

clc,clear

CRUT_ave_all = zeros(1,3);
com_all_all  = zeros(3,317);
com_sel_all  = zeros(3,32);

% read in ensemble member ranked dataset ----------------------------------
 filename = ["Ann_Central_eastern_US_CMIP6_CESM2_LE_CRUT_prcp_1901_1980.nc";...
     "Ann_western_Australia_CMIP6_CESM2_LE_CRUT_prcp_1901_1980.nc";...
     "Ann_UK_CMIP6_CESM2_LE_CRUT_prcp_1901_1980.nc"];

 CESM2filename = ["TSA_CESM2_LE_PRCP_Central_eastern_US_1981_2010_ncl_func.nc";...
     "TSA_CESM2_LE_PRCP_Western_Australia_1981_2010_ncl_func.nc";...
     "TSA_CESM2_LE_PRCP_UK_1981_2010_ncl_func.nc"];

 CMIP6filename = ["TSA_CMIP6_PRCP_Central_eastern_US_1981_2010_ncl_func.nc";...
     "TSA_CMIP6_PRCP_Western_Australia_1981_2010_ncl_func.nc";...
     "TSA_CMIP6_PRCP_UK_1981_2010_ncl_func.nc"];
     
 CRUTfilename = ["TSA_CRUT_PRCP_Central_eastern_US_1981_2010_ncl_func.nc";...
     "TSA_CRUT_PRCP_Western_Australia_1981_2010_ncl_func.nc";...
     "TSA_CRUT_PRCP_UK_1981_2010_ncl_func.nc"];
     

for period = 1:3
    % ensemble order ----------------------------------------------------------
    CMIP6_IPO_ord_temp = ncread(filename(period,:),'CMIP6_AMV_ord');	
    CMIP6_IPO_ord = flip(CMIP6_IPO_ord_temp,1); 
    
    % read in TSA file and mask them ------------------------------------------
    % 220 ensx 144lon x 72lat x 2
    
    CMIP_TSA = ncread(CMIP6filename(period),'TSA');
    CMIP_TSA(CMIP_TSA==0) = NaN;

    CMIP_ave = zeros(1,227)+NaN;
    for ensi = 1:227
        CMIP_temp = squeeze(CMIP_TSA(ensi,:,:,2));
        CMIP_ave(ensi) = nanmean(CMIP_temp,'all')*30.0*12;
    end
    
    % 90 ens x 288lon x 192lat x 2
    CESM_TSA_temp1 = ncread(CESM2filename(period),'TSA');
    CESM_TSA_temp1(CESM_TSA_temp1==0) = NaN;

    CESM_ave = zeros(1,90)+NaN;
    for ensi = 1:90
        CESM_temp = squeeze(CESM_TSA_temp1(ensi,:,:,2));
        CESM_ave(ensi) = nanmean(CESM_temp,'all')*30.0*12;
    end
    
    com_all = [CMIP_ave,CESM_ave];
    com_sel = com_all(CMIP6_IPO_ord(1:32)+1);
    
    % read in observation trend (CRU TS 4.0.3) --------------------------------
    % lat,lon,record 192x288x2
    CRUT_tsa_temp1 = ncread(CRUTfilename(period),'TSA');
    CRUT_tsa_temp1(CRUT_tsa_temp1==0) = NaN;

    CRUT_tsa_temp2 = squeeze(CRUT_tsa_temp1(:,:,2)); 
    CRUT_ave       = nanmean(CRUT_tsa_temp2,'all')*30.0*12;
    
    CRUT_ave_all(period) = CRUT_ave;
    com_all_all(period,:) = com_all;
    com_sel_all(period,:) = com_sel;
    
    %a(period,1) = sqrt(mean((com_all-CRUT_ave).^2));
    %a(period,2) = sqrt(mean((com_sel-CRUT_ave).^2));
end

%%
figure 

xlim([0 4])

hold on
for period = 1:3
   boxplot(com_all_all(period,:)','positions',[period-0.2],'Colors',[150/225, 75/225, 0], 'Symbol','')
   %boxplot(com_all_all(period,:)','positions',[period-0.2], 'Symbol','')
   %A1 = quantile(com_all_all(period,:)',[0.05 0.95]);
   %plot([period-0.2],A1(1),'s','MarkerSize',9,'MarkerEdgeColor',[150/225, 75/225, 0],'MarkerFaceColor','w')
   %plot([period-0.2],A1(2),'s','MarkerSize',9,'MarkerEdgeColor',[150/225, 75/225, 0],'MarkerFaceColor','w')

   boxplot(com_sel_all(period,:)','positions',[period+0.2],'Colors','b', 'Symbol','')
   %boxplot(com_sel_all(period,:)','positions',[period+0.2], 'Symbol','')
   %B1 = quantile(com_sel_all(period,:)',[0.05 0.95]);
   %plot([period+0.2],B1(1),'s','MarkerSize',9,'MarkerEdgeColor','b','MarkerFaceColor','w')
   %plot([period+0.2],B1(2),'s','MarkerSize',9,'MarkerEdgeColor','b','MarkerFaceColor','w')
end
plot([1,2,3],CRUT_ave_all,'ko', 'markersize',7)



ylabel('Precipitatioin Trend (mm/30 yr)')
%title(tit_name(kind))

set(gca,'FontSize',12)
set(gca,'fontname','Times') 

n = findobj(gca,'tag','Outliers');
%n.MarkerEdgeColor = rgb('white');
ylim([-200 300])
xticks([1 2 3 ])
xticklabels({'NE US','Western Australia','UK'})


