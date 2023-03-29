% combine RPC trend and plot figure

clc,clear

%%

i = 20;

filename = ['Ann_tas_CESM2_LE_RPC_global_continent_1901_2014_trend_',num2str(i),'.nc'];
P_corre0 = ncread(filename,'P_corre0');
denomin  = ncread(filename,'denomin');
RPC      = ncread(filename,'RPC');

posi = [1,2; 3,4; 5,6; 7,8; 9,10; 12,14; 16,18; 20,22; 24,26; 28,32; 36,40; ...
    44,48; 52,56; 60,64; 68,72; 76,NaN; 80,NaN; 84,NaN; 86,NaN; 90,NaN];

copyfile Ann_tas_CESM2_LE_RPC_global_continent_1901_2014_trend_20.nc Ann_tas_CESM2_LE_RPC_global_continent_1901_2014_trend.nc 

for i = 1:20
    
    filename = ['Ann_tas_CESM2_LE_RPC_global_continent_1901_2014_trend_',num2str(i),'.nc'];
    P_corre0_temp = ncread(filename,'P_corre0');
    denomin_temp  = ncread(filename,'denomin');
    RPC_temp      = ncread(filename,'RPC');

    posi1 = posi(i,1);
    
    P_corre0(:,posi1) = P_corre0_temp(:,posi1); 
    denomin(:,posi1)  = denomin_temp(:,posi1);
    RPC(:,posi1)      = RPC_temp(:,posi1);    
    
    if ~isnan(posi(i,2))
       posi2 = posi(i,2);
       
       P_corre0(:,posi2) = P_corre0_temp(:,posi2);
       denomin(:,posi2)  = denomin_temp(:,posi2);    
       RPC(:,posi2)      = RPC_temp(:,posi2);       
    end
end

ncwrite('Ann_tas_CESM2_LE_RPC_global_continent_1901_2014_trend.nc','P_corre0',P_corre0)
ncwrite('Ann_tas_CESM2_LE_RPC_global_continent_1901_2014_trend.nc','denomin',denomin)
ncwrite('Ann_tas_CESM2_LE_RPC_global_continent_1901_2014_trend.nc','RPC',RPC)


%%
% size is 500 iteration *90 ens
P_corre0_temp = ncread('Ann_tas_CESM2_LE_RPC_global_continent_1901_2014_trend.nc','P_corre0');
denomin_temp  = ncread('Ann_tas_CESM2_LE_RPC_global_continent_1901_2014_trend.nc','denomin');
RPC_temp      = ncread('Ann_tas_CESM2_LE_RPC_global_continent_1901_2014_trend.nc','RPC');

posi = [1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,22,24,26,28,32,36,40, ...
    44,48,52,56,60,64,68,72,76,80,84,86,90];

P_corre0 = P_corre0_temp(:,posi);
denomin  = denomin_temp(:,posi);
RPC      = RPC_temp(:,posi);

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

for ensj = 2:90
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
h2 = fill(x2, inBetween,color(1,1,:));
set(h2,'EdgeColor',color(1,2,:))
alpha(h2,.5)   

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
%title('CESM2-LE TAS RPC trend 1901-2014','FontSize',14)
xlabel('ensemble size')
ylabel('')

set(gca,'FontSize',15)
set(gca,'fontname','Times') 
set(gcf,'OuterPosition', [0, 1, 650*1.2, 470*1.2]);
xlim([0.5 90.5])
ylim([0.2 1.1])





