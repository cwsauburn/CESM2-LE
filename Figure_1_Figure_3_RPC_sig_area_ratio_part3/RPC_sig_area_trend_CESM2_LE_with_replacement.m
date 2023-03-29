% combine RPC significant area trend and plot figure

clc,clear
%%
i = 21;

filename = ['RPC_global_continent_1901_2014_1000_sample_CRUT_trend_CanESM5_',num2str(i),'.nc'];
Siglv_Can = ncread(filename,'Siglv_final_Had');

copyfile RPC_global_continent_1901_2014_1000_sample_CRUT_trend_CanESM5_21.nc RPC_global_continent_1901_2014_1000_sample_CRUT_trend_CanESM5.nc 

for i = 1:20
    
    posi1 = (i*2-2)+10;
    posi2 = (i*2-1)+10;  
    
    filename = ['RPC_global_continent_1901_2014_1000_sample_CRUT_trend_CanESM5_',num2str(i),'.nc'];
    Siglv_Can_temp = ncread(filename,'Siglv_final_Had');
    
    Siglv_Can(:,posi1) = Siglv_Can_temp(:,posi1);
    Siglv_Can(:,posi2) = Siglv_Can_temp(:,posi2);

end

ncwrite('RPC_global_continent_1901_2014_1000_sample_CRUT_trend_CanESM5.nc','Siglv_final_Had',Siglv_Can)

%%
clc,clear
i = 12;

filename = ['RPC_global_continent_1901_2014_1000_sample_CRUT_trend_IPSL_',num2str(i),'.nc'];
Siglv_IPSL = ncread(filename,'Siglv_final_Had');

copyfile RPC_global_continent_1901_2014_1000_sample_CRUT_trend_IPSL_12.nc RPC_global_continent_1901_2014_1000_sample_CRUT_trend_IPSL.nc
for i = 1:11
    
    posi1 = (i*2-2)+10;
    posi2 = (i*2-1)+10;  
    
    filename = ['RPC_global_continent_1901_2014_1000_sample_CRUT_trend_IPSL_',num2str(i),'.nc'];
    Siglv_IPSL_temp = ncread(filename,'Siglv_final_Had');
    
    Siglv_IPSL(:,posi1) = Siglv_IPSL_temp(:,posi1);
    Siglv_IPSL(:,posi2) = Siglv_IPSL_temp(:,posi2);

end

ncwrite('RPC_global_continent_1901_2014_1000_sample_CRUT_trend_IPSL.nc','Siglv_final_Had',Siglv_IPSL)

%%
clc,clear
i = 11;

filename = ['RPC_global_continent_1901_2014_1000_sample_CRUT_trend_NorCPM1_',num2str(i),'.nc'];
Siglv_final_Had = ncread(filename,'Siglv_final_Had');

copyfile RPC_global_continent_1901_2014_1000_sample_CRUT_trend_NorCPM1_11.nc RPC_global_continent_1901_2014_1000_sample_CRUT_trend_NorCPM1.nc
for i = 1:10
    
    posi1 = (i*2-2)+10;
    posi2 = (i*2-1)+10;    
    
    filename = ['RPC_global_continent_1901_2014_1000_sample_CRUT_trend_NorCPM1_',num2str(i),'.nc'];
    Siglv_final_Had_temp = ncread(filename,'Siglv_final_Had');
    
    
    Siglv_final_Had(:,posi1) = Siglv_final_Had_temp(:,posi1);
    Siglv_final_Had(:,posi2) = Siglv_final_Had_temp(:,posi2);
end

ncwrite('RPC_global_continent_1901_2014_1000_sample_CRUT_trend_NorCPM1.nc','Siglv_final_Had',Siglv_final_Had)

%%

i = 1;
% 500*90
filename = ['Ann_tas_CESM2-LE_RPC_global_continent_1901_2014_1000_sample_CRUT_trend_',num2str(i),'.nc'];
Siglv_Can = ncread(filename,'Siglv_final_Had');

copyfile Ann_tas_CESM2-LE_RPC_global_continent_1901_2014_1000_sample_CRUT_trend_1.nc Ann_tas_CESM2-LE_RPC_global_continent_1901_2014_1000_sample_CRUT_trend.nc 

for i = 2:34
    
    posi = 10+(i-1)*2;
    
    filename = ['Ann_tas_CESM2-LE_RPC_global_continent_1901_2014_1000_sample_CRUT_trend_',num2str(i),'.nc'];
    Siglv_Can_temp = ncread(filename,'Siglv_final_Had');
    
    Siglv_Can(:,posi) = Siglv_Can_temp(:,posi);
end

posi = [80, 84, 88, 90];
for i = 35: 38

    filename = ['Ann_tas_CESM2-LE_RPC_global_continent_1901_2014_1000_sample_CRUT_trend_',num2str(i),'.nc'];
    Siglv_Can_temp = ncread(filename,'Siglv_final_Had');
    
    Siglv_Can(:,posi(i-34)) = Siglv_Can_temp(:,posi(i-34));
end

ncwrite('Ann_tas_CESM2-LE_RPC_global_continent_1901_2014_1000_sample_CRUT_trend.nc','Siglv_final_Had',Siglv_Can)

%%
clc,clear
i = 22;

filename        = ['RPC_global_continent_1901_2014_1000_sample_CRUT_trend_',num2str(i),'.nc'];
Siglv_final_Had = ncread(filename,'Siglv_final_Had');

copyfile RPC_global_continent_1901_2014_1000_sample_CRUT_trend_22.nc RPC_global_continent_1901_2014_1000_sample_CRUT_trend.nc
for i = 1:21
    
    posi1 = i*2-1;
    posi2 = i*2;    
    
    filename = ['RPC_global_continent_1901_2014_1000_sample_CRUT_trend_',num2str(i),'.nc'];
    Siglv_final_Had_temp = ncread(filename,'Siglv_final_Had');
    
    
    Siglv_final_Had(:,posi1) = Siglv_final_Had_temp(:,posi1);
    Siglv_final_Had(:,posi2) = Siglv_final_Had_temp(:,posi2);
end

ncwrite('RPC_global_continent_1901_2014_1000_sample_CRUT_trend.nc','Siglv_final_Had',Siglv_final_Had)

%%

i = 1;
% 500*90
filename = ['Ann_tas_CESM2-LE_RPC_global_continent_1901_2014_1000_sample_CRUT_trend_',num2str(i),'.nc'];
Siglv_Can = ncread(filename,'Siglv_final_Had');

copyfile Ann_tas_CESM2-LE_RPC_global_continent_1901_2014_1000_sample_CRUT_trend_1.nc Ann_tas_CESM2-LE_RPC_global_continent_1901_2014_1000_sample_CRUT_trend.nc 

for i = 2:34
    
    posi = 10+(i-1)*2;
    
    filename = ['Ann_tas_CESM2-LE_RPC_global_continent_1901_2014_1000_sample_CRUT_trend_',num2str(i),'.nc'];
    Siglv_Can_temp = ncread(filename,'Siglv_final_Had');
    
    Siglv_Can(:,posi) = Siglv_Can_temp(:,posi);
end

posi = [80, 84, 88, 90];
for i = 35: 38

    filename = ['Ann_tas_CESM2-LE_RPC_global_continent_1901_2014_1000_sample_CRUT_trend_',num2str(i),'.nc'];
    Siglv_Can_temp = ncread(filename,'Siglv_final_Had');
    
    Siglv_Can(:,posi(i-34)) = Siglv_Can_temp(:,posi(i-34));
end

ncwrite('Ann_tas_CESM2-LE_RPC_global_continent_1901_2014_1000_sample_CRUT_trend.nc','Siglv_final_Had',Siglv_Can)

%%

Siglv_Can       = ncread('RPC_global_continent_1901_2014_1000_sample_CRUT_CanESM5_trend.nc','Siglv_final_Had');
Siglv_IPSL      = ncread('RPC_global_continent_1901_2014_1000_sample_CRUT_IPSL_CM6A_LR_trend.nc','Siglv_final_Had');
Siglv_final_Had = ncread('RPC_global_continent_1901_2014_1000_sample_CRUT_trend.nc','Siglv_final_Had');
Siglv_Nor       = ncread('RPC_global_continent_1901_2014_1000_sample_CRUT_NorCPM1_trend.nc','Siglv_final_Had');
Siglv_CESM_temp = ncread('Ann_tas_CESM2-LE_RPC_global_continent_1901_2014_1000_sample_CRUT_trend_2_5_degree.nc','Siglv_final_Had');

Siglv_CESM = zeros(500,38)+NaN;
xCESM = [10:2:76,80, 84, 88, 90];

for i = 1:38
    Siglv_CESM(:,i) = Siglv_CESM_temp(:,xCESM(i));
end

% calculating mean value and averaged band
mean_Can = nanmean(Siglv_Can(:,10:50),1);
Can_bnd1 = mean_Can+2*std(Siglv_Can(:,10:50),0,1,'omitnan');
Can_bnd2 = mean_Can-2*std(Siglv_Can(:,10:50),0,1,'omitnan');

mean_IPSL = nanmean(Siglv_IPSL(:,10:32),1);
IPSL_bnd1 = mean_IPSL+2*std(Siglv_IPSL(:,10:32),0,1,'omitnan');
IPSL_bnd2 = mean_IPSL-2*std(Siglv_IPSL(:,10:32),0,1,'omitnan');

mean_CMIP = nanmean(Siglv_final_Had,1);
CMIP_bnd1 = mean_CMIP+2*std(Siglv_final_Had,0,1,'omitnan');
CMIP_bnd2 = mean_CMIP-2*std(Siglv_final_Had,0,1,'omitnan');

mean_Nor  = nanmean(Siglv_Nor(:,10:30),1);
Nor_bnd1  = mean_Nor+2*std(Siglv_Nor(:,10:30),0,1,'omitnan');
Nor_bnd2  = mean_Nor-2*std(Siglv_Nor(:,10:30),0,1,'omitnan');

mean_CESM = nanmean(Siglv_CESM,1);
CESM_bnd1 = mean_CESM+2*std(Siglv_CESM,0,1,'omitnan');
CESM_bnd2 = mean_CESM-2*std(Siglv_CESM,0,1,'omitnan');

%%

% ploting figure
xCan  = 10:50;
xIPSL = 10:32;
xNor  = 10:30;

x1 = 10:2:50;
x2 = 55:5:100;
x3 = 110:10:220;

xCMIP = [x1,x2,x3];

color(1,1,:) = [0.5843,0.8157,0.9882]; 
color(1,2,:) = [0,0,1];
color(2,1,:) = [1,1,0];
color(2,2,:) = [0.93,0.69,0.13];
color(3,1,:) = [0.9,0.9,0.9];
color(3,2,:) = [0,0,0];
color(4,1,:) = [0.74,0.44,0.92];
color(4,2,:) = [0.72,0.27,1.00];
color(5,1,:) = [0.98,0.73,0.77];
color(5,2,:) = [0.64,0.08,0.18];

figure
box on
hold on
grid on

plot(xNor,mean_Nor,'color',color(4,2,:),'Linewidth',1.8)
plot(xCan, mean_Can,'color', color(1,2,:),'Linewidth',1.8)
plot(xIPSL,mean_IPSL,'color',color(2,2,:),'Linewidth',1.8)
plot(xCMIP,mean_CMIP,'color',color(3,2,:),'Linewidth',1.8)
plot(xCESM,mean_CESM,'color',color(5,2,:),'Linewidth',1.8)

xNorB     = [xNor, fliplr(xNor)];
inBetween = [Nor_bnd2,fliplr(Nor_bnd1)];
h4        = fill(xNorB, inBetween,color(4,1,:));
set(h4,'EdgeColor',color(4,2,:))
alpha(h4,.2)

xCanB     = [xCan, fliplr(xCan)];
inBetween = [Can_bnd2,fliplr(Can_bnd1)];
h1        = fill(xCanB, inBetween,color(1,1,:));
set(h1,'EdgeColor',color(1,2,:))
alpha(h1,.2) 

xIPSLB     = [xIPSL, fliplr(xIPSL)];
inBetween  = [IPSL_bnd2,fliplr(IPSL_bnd1)];
h2         = fill(xIPSLB, inBetween,color(2,1,:));
set(h2,'EdgeColor',color(2,2,:))
alpha(h2,.2)   

xCMIPB     = [xCMIP, fliplr(xCMIP)];
inBetween  = [CMIP_bnd2,fliplr(CMIP_bnd1)];
h3         = fill(xCMIPB, inBetween,color(3,1,:));
set(h3,'EdgeColor',color(3,2,:))
alpha(h3,.2)

xCESMB     = [xCESM, fliplr(xCESM)];
inBetween  = [CESM_bnd2,fliplr(CESM_bnd1)];
h4         = fill(xCESMB, inBetween,color(5,1,:));
set(h4,'EdgeColor',color(5,2,:))
alpha(h4,.2)

%plot(xNor, mean_Nor,'color',color(4,2,:),'Linewidth',1.8)
%plot(xCan, mean_Can,'color', color(1,2,:),'Linewidth',1.8)
%plot(xIPSL,mean_IPSL,'color',color(2,2,:),'Linewidth',1.8)
%plot(xCMIP,mean_CMIP,'color',color(3,2,:),'Linewidth',1.8)
%plot(xCESM,mean_CESM,'color',color(5,2,:),'Linewidth',1.8)

hold off

legend('NorCPM1(30)','CanESM5(50)','IPSL-CM6A-LR(32)','CMIP6(220)','CESM2-LE(90)')
%title('TAS RPC significant area ratio 1901-2014','FontSize',14)
xlabel('ensemble size')
%ylabel('Fraction of global land area where the null hypothesis is rejected')

set(gca,'FontSize',15)
set(gca,'fontname','Times') 
set(gcf,'OuterPosition', [0, 1, 650*1.4, 470]);
xlim([9.5 220.5])
ylim([0.05 0.80])




