% combine RPC trend and plot figure

clc,clear

%%

P_corre0 = zeros(500,34);
denomin  = zeros(500,34);
RPC      = zeros(500,34);

filename = 'Ann_tas_CMIP6_CanESM5_RPC_global_continent_1951_2014_trend_1.nc';
P_corre0_temp = ncread(filename,'P_corre0');
denomin_temp  = ncread(filename,'denomin');
RPC_temp      = ncread(filename,'RPC');

for i = 1:9
    P_corre0(:,i*2-1) = P_corre0_temp(:,i); 
    denomin(:,i*2-1)  = denomin_temp(:,i);
    RPC(:,i*2-1)      = RPC_temp(:,i);   
end

filename = 'Ann_tas_CMIP6_CanESM5_RPC_global_continent_1951_2014_trend_4.nc';
P_corre0_temp = ncread(filename,'P_corre0');
denomin_temp  = ncread(filename,'denomin');
RPC_temp      = ncread(filename,'RPC');
    
for i = 1:9
    P_corre0(:,i*2) = P_corre0_temp(:,i); 
    denomin(:,i*2)  = denomin_temp(:,i);
    RPC(:,i*2)      = RPC_temp(:,i);   
end

filename = 'Ann_tas_CMIP6_CanESM5_RPC_global_continent_1951_2014_trend_2.nc';
P_corre0_temp = ncread(filename,'P_corre0');
denomin_temp  = ncread(filename,'denomin');
RPC_temp      = ncread(filename,'RPC');
    
P_corre0(:,10+9:18+9) = P_corre0_temp(:,10:18); 
denomin(:,10+9:18+9)  = denomin_temp(:,10:18);
RPC(:,10+9:18+9)      = RPC_temp(:,10:18);   

filename = 'Ann_tas_CMIP6_CanESM5_RPC_global_continent_1951_2014_trend_3.nc';
P_corre0_temp = ncread(filename,'P_corre0');
denomin_temp  = ncread(filename,'denomin');
RPC_temp      = ncread(filename,'RPC');

P_corre0(:,19+9:25+9) = P_corre0_temp(:,19:25); 
denomin(:,19+9:25+9)  = denomin_temp(:,19:25);
RPC(:,19+9:25+9)      = RPC_temp(:,19:25);   

P_corre0(:,26+9) = P_corre0_temp(:,26); 
denomin(:,26+9)  = denomin_temp(:,26);
RPC(:,26+9)      = RPC_temp(:,26);   

%%

posi = [1:18,19:2:49,50];

P_corr_bnd1 = nanmean(P_corre0,1)+2*std(P_corre0,0,1,'omitnan');
P_corr_bnd2 = nanmean(P_corre0,1)-2*std(P_corre0,0,1,'omitnan');

denomin_bnd1 = nanmean(denomin,1)+2*std(denomin,0,1,'omitnan');
denomin_bnd2 = nanmean(denomin,1)-2*std(denomin,0,1,'omitnan');

RPC_bnd1 = nanmean(RPC,1)+2*std(RPC,0,1,'omitnan');
RPC_bnd2 = nanmean(RPC,1)-2*std(RPC,0,1,'omitnan');

%{
% calculate signal to total ratio 95% significance level critical value
F_CV = ncread('f_critical_value_0.05.nc','F_CV');
N1 = 2014-1901+1;
N2 = 2014-1951+1;
STR_CV = zeros(size(F_CV))+NaN;

for ensj = 2:32
    SNR1 = F_CV(ensj,1)*(N1-1)/(N1*(ensj-1));
    STR_CV(ensj,1) = SNR1/(SNR1+1);

    SNR2 = F_CV(ensj,2)*(N2-1)/(N2*(ensj-1));
    STR_CV(ensj,2) = SNR2/(SNR2+1);
end

STR_CV = sqrt(STR_CV);
%}

%%
% plot Normalized streamflow RMSD for 4 months ---------------------------% 
color(1,1,:) = [0.5843,0.8157,0.9882]; 
color(1,2,:) = [0,0,1];
color(2,1,:) = [0.925,0.839,0.839];
color(2,2,:) = [1,0,0];

x3 = 1:90;

figure
box on
hold on
grid on

plot(posi,nanmean(RPC,1),'color','k','Linewidth',1.8)
plot(posi,nanmean(P_corre0,1),'color',color(1,2,:),'Linewidth',1.8)
plot(posi,nanmean(denomin,1),'color',color(2,2,:),'Linewidth',1.8)
%plot(x3,STR_CV(1:90,1),'--','color',[0.58,0.58,0.58],'Linewidth',1.8)

x2 = [posi, fliplr(posi)];
inBetween = [RPC_bnd2,fliplr(RPC_bnd1)];
h2 = fill(x2, inBetween,[0.9,0.9,0.9]);
set(h2,'EdgeColor','k')
alpha(h2,.5) 

x2 = [posi, fliplr(posi)];
inBetween = [P_corr_bnd2,fliplr(P_corr_bnd1)];
h3 = fill(x2, inBetween,color(1,1,:));
set(h3,'EdgeColor',color(1,2,:))
alpha(h3,.5)   

inBetween = [denomin_bnd2,fliplr(denomin_bnd1)];
h1 = fill(x2, inBetween,color(2,1,:));
set(h1,'EdgeColor',color(2,2,:))
alpha(h1,.5)

plot(posi,nanmean(RPC,1),'color','k','Linewidth',1.8)
plot(posi,nanmean(P_corre0,1),'color',color(1,2,:),'Linewidth',1.8)
plot(posi,nanmean(denomin,1),'color',color(2,2,:),'Linewidth',1.8)
%plot(x3,STR_CV(1:90,1),'--','color',[0.58,0.58,0.58],'Linewidth',1.8)

hold off

%legend('RPC','EM-OBS Correl','Sig./Total  Ratio','Ratio CV(0.95)')
legend('RPC','EM-OBS Correl','Sig./Total  Ratio')
title('CanESM5 TAS RPC trend 1951-2014','FontSize',14)
xlabel('ensemble size')
ylabel('')

set(gca,'FontSize',15)
set(gca,'fontname','Times') 
set(gcf,'OuterPosition', [0, 1, 650*1.2, 470*1.2]);
xlim([0.5 50.5])
ylim([0.2 1.1])





