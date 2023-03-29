clc,clear

% plot GISS and CMIP6-LE figure like Smith et al 2020
% ensemble 220+90* year 114

% filename = 'Ann_India_CMIP6_CESM2_LE_GISS_tas_1901_1950.nc';
% filename = 'Ann_India_CMIP6_CESM2_LE_GISS_tas_1901_1960.nc';
% filename = 'Ann_India_CMIP6_CESM2_LE_GISS_tas_1901_1970.nc';
% filename = 'Ann_India_CMIP6_CESM2_LE_GISS_tas_1901_1980.nc';
filename = 'Ann_India_CMIP6_CESM2_LE_GISS_tas_1901_2014.nc';

% sel_year = 50:80;
% sel_year = 60:90;
% sel_year = 70:100;
% sel_year = 80:110;
 sel_year = 5:109;

CMIP6_IPO = ncread(filename,'CMIP6_AMV');
CMIP6_IPO = CMIP6_IPO(:,sel_year);
cli_mod = nanmean(nanmean(CMIP6_IPO,2),1);

% year 114 
GISS_IPO  = ncread(filename,'GISS_AMV');
GISS_IPO  = GISS_IPO(sel_year);

% ensemble order
CMIP6_IPO_ord_temp = ncread(filename,'CMIP6_AMV_ord');	
CMIP6_IPO_ord = flip(CMIP6_IPO_ord_temp,1); 

CMIP6_CI1 = zeros(2,size(CMIP6_IPO,2)) + 0;
CMIP6_CI1(1,:) = nanmean(CMIP6_IPO(CMIP6_IPO_ord(1:31)+1,:),1)+2.086*nanstd(CMIP6_IPO(CMIP6_IPO_ord(1:31)+1,:),1);
CMIP6_CI1(2,:) = nanmean(CMIP6_IPO(CMIP6_IPO_ord(1:31)+1,:),1)-2.086*nanstd(CMIP6_IPO(CMIP6_IPO_ord(1:31)+1,:),1);

CMIP6_CI2 = zeros(2,size(CMIP6_IPO,2)) + 0;
CMIP6_CI2(1,:) = nanmean(CMIP6_IPO,1)+1.66*nanstd(CMIP6_IPO,1);
CMIP6_CI2(2,:) = nanmean(CMIP6_IPO,1)-1.66*nanstd(CMIP6_IPO,1);

CMIP6_mean1 = nanmean(CMIP6_IPO(CMIP6_IPO_ord(1:31)+1,:),1);
CMIP6_mean2 = nanmean(CMIP6_IPO,1);

% RMSD1 = sqrt(nanmean((CMIP6_mean1-GISS_IPO').^2));
% RMSD2 = sqrt(nanmean((CMIP6_mean2-GISS_IPO').^2));

% calculate correlation ACC&RPC ----------------------------
ACC_1 = corr(CMIP6_mean1',GISS_IPO);
ACC_2 = corr(CMIP6_mean2',GISS_IPO);

vars1 = zeros(1,31)+NaN;
for ens = 1:31
    vars1(ens) = nansum((CMIP6_IPO(CMIP6_IPO_ord(ens)+1,:)-cli_mod).^2)/103;
end

vars2 = zeros(1,220+90)+NaN;
for ens = 1:310
    vars2(ens) = nansum((CMIP6_IPO(ens,:)-cli_mod).^2)/103;
end

s2t1 = sqrt(var(CMIP6_mean1)/nanmean(vars1));
s2t2 = sqrt(var(CMIP6_mean2)/nanmean(vars2));

RPC_1 = ACC_1/s2t1;
RPC_2 = ACC_2/s2t2;


%%
xtime = 1905:2009;

% Plot seleted time series plot ----------------------------
figure
hold on
grid on
box on

color(1,:) = [0.925,0.839,0.839];
color(2,:) = [1,0,0];

plot(xtime,GISS_IPO,'k','LineWidth',0.5)
plot(xtime,CMIP6_mean1,'r','LineWidth',0.5)

xtimeflip = [xtime, fliplr(xtime)];
inBetween = [CMIP6_CI1(2,:),fliplr(CMIP6_CI1(1,:))];
h1        = fill(xtimeflip,inBetween,color(1,:));
set(h1,'EdgeColor',color(2,:))
alpha(h1,.5) 

plot(xtime,GISS_IPO,'k','LineWidth',2)
plot(xtime,CMIP6_mean1,'r','LineWidth',2)

COEF1 = polyfit(sel_year,GISS_IPO,1);
COEF2 = polyfit(sel_year,CMIP6_mean1,1);

plot(xtime,COEF1(1)*(xtime-1905)+COEF1(2),'k',linestyle='--',LineWidth=1.2)
plot(xtime,COEF2(1)*(xtime-1905)+COEF2(2),'r',linestyle='--',LineWidth=1.2)

hold off
legend('GISS','CMIP6+CESM2-LE(31)')
xlabel('year')
ylabel('temperature anomaly °C')

set(gca,'FontSize',15)
set(gca,'fontname','Times') 
set(gcf,'OuterPosition', [0, 1, 650*1.5, 450]);
xlim([1905 2009])
%ylim([-1.0 2.5])

% Plot unselected time series plot ---------------------------
figure
hold on
grid on
box on
plot(xtime,GISS_IPO,'k','LineWidth',0.5)
plot(xtime,CMIP6_mean2,'r','LineWidth',0.5)

xtimeflip = [xtime, fliplr(xtime)];
inBetween = [CMIP6_CI2(2,:),fliplr(CMIP6_CI2(1,:))];
h2        = fill(xtimeflip,inBetween,color(1,:));
set(h2,'EdgeColor',color(2,:))
alpha(h2,.5) 

plot(xtime,GISS_IPO,'k','LineWidth',2)
plot(xtime,CMIP6_mean2,'r','LineWidth',2)

COEF3 = polyfit(sel_year,CMIP6_mean2,1);
plot(xtime,COEF1(1)*(xtime-1905)+COEF1(2),'k',linestyle='--',LineWidth=1.2)
plot(xtime,COEF3(1)*(xtime-1905)+COEF3(2),'r',linestyle='--',LineWidth=1.2)

hold off
legend('GISS','CMIP6+CESM2-LE(310)')
xlabel('year')
ylabel('temperature anomaly °C')

set(gca,'FontSize',15)
set(gca,'fontname','Times') 
set(gcf,'OuterPosition', [0, 1, 650*1.5, 450]);
xlim([1905 2009])
%ylim([-1.0 2.5])



