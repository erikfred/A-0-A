% get_Tdep.m
%
% Fiddling around with temperature dependence of the pressure gauges and
% barometer
%

clear; close all

load ../pressure_data/POBS1.mat

% extract sample fields
ttest=data.ta(data.lNormState1);
ptest=data.a1(data.lNormState1);
Ttest=data.Ta1(data.lNormState1);

% trim beginning (large T fluctuation)
ttest(1:4e4)=[]; ptest(1:4e4)=[]; Ttest(1:4e4)=[];

% make a filter
[b,a]=butter(3,2/3600,'high');
ptest2=filtfilt(b,a,ptest);

temp=[];
for i=10^5:10^5:3.1*10^7
    tinv=ttest(i:i+10^5)-ttest(i);
    Gp=[Ttest(i:i+10^5),ones(size(Ttest(i:i+10^5))),tinv];
    mp=inv(Gp'*Gp)*Gp'*ptest(i:i+10^5);
    mtest=Gp*mp;
    temp=[temp;mp(1)];
    % disp(["Temperature Dependence (" num2str(round(i/60/60/24)) " days) = "])
    % disp([num2str(mp(1)) ' hPa/C' newline])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3); clf;
subplot(211); hold on
plot(calInfoAll1.t0,calInfoAll1.pCal,'o','markersize',10,'linewidth',2)
plot(calInfoAll2.t0,calInfoAll2.pCal,'s','markersize',10,'linewidth',2)
plot(barInfoAll.t0p,zeros(size(barInfoAll.pCal))+865,'xk','markersize',7,'linewidth',1)
ylabel('P (hPa)')
yyaxis right
plot(dataDec.ta(3600:end),dataDec.Ta1(3600:end),'color',[0 0.4470 0.7410],'linewidth',0.1)
plot(dataDec.ta(3600:end),dataDec.Ta2(3600:end),'color','r','linewidth',0.1)
ylabel('T (C)')
legend('Gauge 1','Gauge2','Flat Barometer','location','northwest')
datetick('x',6)
title([POBS_dir(i).name ' calibrations'])
set(gca,'fontsize',14)
box on; grid on
subplot(212); hold on
barcor1=calInfoAll1.pCal-zeros(size(barInfoAll.pCal));
barcor2=calInfoAll2.pCal-zeros(size(barInfoAll.pCal));
plot(calInfoAll1.t0,barcor1-barcor1(end),'o','markersize',10,'linewidth',2)
plot(calInfoAll2.t0,barcor2-barcor2(end),'s','markersize',10,'linewidth',2)
ylabel('P (hPa)')
legend('Gauge 1','Gauge2','location','northeast')
datetick('x',6)
title('Flat barometer calibrations')
set(gca,'fontsize',14)
box on; grid on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print(['../figures/calibrations/' POBS_dir(i).name '_cal_summary_nobarometer'],'-dpng','-r300')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ttest=calInfoAll2.t0p(2:end); tinv=ttest-ttest(1);
ptest=calInfoAll2.pCal(2:end);
Ttest=barInfoAll.T(2:end)';

G=[Ttest,tinv,ones(size(tinv))];
m=inv(G'*G)*G'*ptest;

ttest1=calInfoAll1.t0p(2:end); tinv1=ttest1-ttest1(1);
ptest1=calInfoAll1.pCal(2:end);
Ttest1=barInfoAll.T(2:end)';

G1=[Ttest1,tinv1,ones(size(tinv1))];
% add in exponential
lamda_list=linspace(1,180,1000);
for jj=1:1000
    gexp=exp(-tinv1/lamda_list(jj)); gexp(gexp<10^-7)=0;
    G1e=[G1,gexp];
    m1e(:,jj)=inv(G1e'*G1e)*G1e'*ptest1;
    e_fit(:,jj)=G1e*m1e(:,jj);
    stds(jj)=std(ptest1-e_fit(:,jj));
end
[~,imin]=min(stds);
e_cor=ptest1-e_fit(:,imin);
m1=[m1e(:,imin);lamda_list(imin)];