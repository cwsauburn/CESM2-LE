% combine RPC trend and plot figure

clc,clear

%%
i = 1;

filename = ['Ann_tas_CMIP6_RPC_global_continent_1951_2014_trend_',num2str(i),'.nc'];
P_corre0 = ncread(filename,'P_corre0');
denomin  = ncread(filename,'denomin');
RPC      = ncread(filename,'RPC');

for i = 2:11
    
    posi1 = i*2+24;
    posi2 = i*2+24+1;  
    
    filename = ['Ann_tas_CMIP6_RPC_global_continent_1951_2014_trend_',num2str(i),'.nc'];
    P_corre0_temp = ncread(filename,'P_corre0');
    denomin_temp  = ncread(filename,'denomin');
    RPC_temp      = ncread(filename,'RPC');
    
    %P_corre0((i*25-24):i*25,:) = P_corre0_temp((i*25-24):i*25,:);
    %P_corre0((i*25-24):i*25,:) = P_corre0_temp(1:25,:);
    
    P_corre0(:,posi1:posi2) = P_corre0_temp(:,posi1:posi2);
    denomin(:,posi1:posi2)  = denomin_temp(:,posi1:posi2);
    RPC(:,posi1:posi2)  = RPC_temp(:,posi1:posi2);
end

i = 12;
filename = ['Ann_tas_CMIP6_RPC_global_continent_1951_2014_trend_',num2str(i),'.nc'];
P_corre0_temp = ncread(filename,'P_corre0');
denomin_temp  = ncread(filename,'denomin');
RPC_temp      = ncread(filename,'RPC');

P_corre0(:,48) = P_corre0_temp(:,48);
denomin(:,48)  = denomin_temp(:,48);
RPC(:,48)  = RPC_temp(:,48);

ncwrite('Ann_tas_CMIP6_RPC_global_continent_1951_2014_trend.nc','P_corre0',P_corre0)
ncwrite('Ann_tas_CMIP6_RPC_global_continent_1951_2014_trend.nc','denomin',denomin)
ncwrite('Ann_tas_CMIP6_RPC_global_continent_1951_2014_trend.nc','RPC',RPC)

%%

% size is 500 iteration *56 ens
P_corre0 = ncread('Ann_tas_CMIP6_RPC_global_continent_1951_2014_trend.nc','P_corre0');
denomin  = ncread('Ann_tas_CMIP6_RPC_global_continent_1951_2014_trend.nc','denomin');
RPC      = ncread('Ann_tas_CMIP6_RPC_global_continent_1951_2014_trend.nc','RPC');

P_corr_bnd1 = nanmean(P_corre0,1)+2*std(P_corre0,0,1,'omitnan');
P_corr_bnd2 = nanmean(P_corre0,1)-2*std(P_corre0,0,1,'omitnan');

denomin_bnd1 = nanmean(denomin,1)+2*std(denomin,0,1,'omitnan');
denomin_bnd2 = nanmean(denomin,1)-2*std(denomin,0,1,'omitnan');

RPC_bnd1 = nanmean(RPC,1)+2*std(RPC,0,1,'omitnan');
RPC_bnd2 = nanmean(RPC,1)-2*std(RPC,0,1,'omitnan');

% calculate signal to total ratio 95% significance level critical value
F_CV = ncread('f_critical_value_0.05.nc','F_CV');
N1 = 2014-1901+1;
N2 = 2014-1951+1;
STR_CV = zeros(size(F_CV))+NaN;

for ensj = 2:228
    SNR1 = F_CV(ensj,1)*(N1-1)/(N1*(ensj-1));
    STR_CV(ensj,1) = SNR1/(SNR1+1);

    SNR2 = F_CV(ensj,2)*(N2-1)/(N2*(ensj-1));
    STR_CV(ensj,2) = SNR2/(SNR2+1);
end

STR_CV = sqrt(STR_CV);

%%
% plot Normalized streamflow RMSD for 4 months ---------------------------% 
color(1,1,:) = [0.5843,0.8157,0.9882]; 
color(1,2,:) = [0,0,1];
color(2,1,:) = [0.925,0.839,0.839];
color(2,2,:) = [1,0,0];

x1 = 1:2:51;
x2 = 55:5:100;
x3 = 110:10:220;

x   = [x1,x2,x3];
xCV = 1:220;

figure
box on
hold on
grid on

plot(x,nanmean(RPC,1),'color','k','Linewidth',1.8)
plot(x,nanmean(P_corre0,1),'color',color(1,2,:),'Linewidth',1.8)
plot(x,nanmean(denomin,1),'color',color(2,2,:),'Linewidth',1.8)
% plot(xCV,STR_CV(1:220,2),'--','color',[0.58,0.58,0.58],'Linewidth',1.8)

x2 = [x, fliplr(x)];
inBetween = [RPC_bnd2,fliplr(RPC_bnd1)];
h2 = fill(x2, inBetween,[0.9,0.9,0.9]);
set(h2,'EdgeColor','k')
alpha(h2,.5) 

x2 = [x, fliplr(x)];
inBetween = [P_corr_bnd2,fliplr(P_corr_bnd1)];
h2 = fill(x2, inBetween,color(1,1,:));
set(h2,'EdgeColor',color(1,2,:))
alpha(h2,.5)   

inBetween = [denomin_bnd2,fliplr(denomin_bnd1)];
h1 = fill(x2, inBetween,color(2,1,:));
set(h1,'EdgeColor',color(2,2,:))
alpha(h1,.5)

plot(x,nanmean(RPC,1),'color','k','Linewidth',1.8)
plot(x,nanmean(P_corre0,1),'color',color(1,2,:),'Linewidth',1.8)
plot(x,nanmean(denomin,1),'color',color(2,2,:),'Linewidth',1.8)
% plot(xCV,STR_CV(1:220,2),'--','color',[0.58,0.58,0.58],'Linewidth',1.8)
hold off

legend('RPC','EM-OBS Correl','Sig./Total  Ratio')
% legend('RPC','EM-OBS Correl','Sig./Total  Ratio','Ratio CV(0.95)')
title('CMIP6 TAS RPC trend 1951-2014','FontSize',14)
xlabel('ensemble size')
ylabel('')

set(gca,'FontSize',15)
set(gca,'fontname','Times') 
set(gcf,'OuterPosition', [0, 1, 650*1.2, 470*1.2]);
xlim([0.5 220.5])
 ylim([0.2 1.1])






