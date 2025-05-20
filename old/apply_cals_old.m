% apply_cals.m
%
% Apply drift calibrations to pressure observations by two methods:
%   1) Generate a smooth model by fitting the calibrations to the
%   functional form y = ax + b + c*exp(-t/d)
%   2) Generate a linearly-interpolated model to fill in gaps between
%   calibration points
%

clear; close all

POBS_dir=dir('../stitched_data/');
POBS_list={POBS_dir.name}';
file_check=cellfun(@(v)v(1),POBS_list);
i_list=find(eq(file_check,'P')); % subselects only files begining with 'P'
n=0;
for i=1:length(i_list)
    k=i_list(i);
    load(['../pressure_data/' POBS_dir(k).name],'calInfoAll1','calInfoAll2','barInfoAll')
    if isempty(calInfoAll1)
        continue
    elseif length(calInfoAll1.i0)<=4
        continue
    end

    load([POBS_dir(k).folder '/' POBS_dir(k).name])

    if i==6 % quick fix for POBS-02 bad barometer
        barInfoAll.pCal=0;
    elseif i==10 % POBS-07 has additional characters we don't need
        POBS_dir(k).name(6:end-4)=[];
    end

    % % trim off first 2 months of data
    % dataf.tf(1:2*30*24)=[]; dataf.p1f(1:2*30*24)=[]; dataf.p2f(1:2*30*24)=[];

    % gather variables
    t=barInfoAll.t0p;
    bp=barInfoAll.pCal;
    bT=barInfoAll.T;
    p1=calInfoAll1.pCal;
    p2=calInfoAll2.pCal;
    T1=calInfoAll1.T;
    T2=calInfoAll2.T;

    % trim as necesary
    if i==1 % offset in last 10 cals
        t(end-10:end)=[]; bp(end-10:end)=[]; bT(end-10:end)=[];
        p1(end-10:end)=[]; p2(end-10:end)=[]; T1(end-10:end)=[]; T2(end-10:end)=[];
    elseif i==4
        t(end-9:end)=[]; bp(end-9:end)=[]; bT(end-9:end)=[];
        p1(end-9:end)=[]; p2(end-9:end)=[]; T1(end-9:end)=[]; T2(end-9:end)=[];
    elseif i==5
        p1(end-5:end)=p2(end-5:end)+16.3; % approximate fix that allows code to
                                            % work until I truncate later
    end

    % correct barometer for temperature as y = ax + bT + c
    Gb=[t-t(1),bT-median(bT),ones(size(t))];
    mb=inv(Gb'*Gb)*Gb'*bp;
    bp=bp-(bT-median(bT))*mb(2);

    % correct pressure gauges for barometer
    p1=p1-bp;
    p2=p2-bp;

    % generate reasonable range of exponential time constants
    lamda_list=linspace(1/max(t-t(1)),180/max(t-t(1)),1000);

    % better time basis for inversions
    tinv=(t-t(1))/max(t-t(1));

    m1=[]; p1_fit=[]; stds1=[];
    m2=[]; p2_fit=[]; stds2=[];
    % assume form y = ax + bT + c + d*exp(-t/f)
    for jj=1:1000
        % construct matrix for inversion
        G1=[tinv,T1-median(T1),ones(size(tinv))];
        G2=[tinv,T2-median(T2),ones(size(tinv))];
        gexp=exp(-tinv/lamda_list(jj)); gexp(gexp<10^-7)=0;
        % Gauge 1
        G1=[G1,gexp];
        m1(:,jj)=inv(G1'*G1)*G1'*p1;
        p1_fit(:,jj)=G1*m1(:,jj);
        stds1(jj)=std(p1-p1_fit(:,jj));
        % Gauge 2
        G2=[G2,gexp];
        m2(:,jj)=inv(G2'*G2)*G2'*p2;
        p2_fit(:,jj)=G2*m2(:,jj);
        stds2(jj)=std(p2-p2_fit(:,jj));
    end
    [~,imin1]=min(stds1);
    [~,imin2]=min(stds2);

    % quality control figures
    figure(13); clf
    subplot(311); hold on
    plot(linspace(1,180,1000),stds1,'.','markersize',5)
    ylabel('RMS (cm)')
    xlabel('time constant (d)')
    title([POBS_dir(k).name(1:end-4) ' Gauge 1'])
    set(gca,'fontsize',14)
    box on; grid on
    subplot(312); hold on
    plot(t,p1,'o','markersize',10,'linewidth',2)
    plot(t,p1_fit(:,imin1),'linewidth',1)
    legend('cals','fit','location','northeast')
    ylabel('P (cm)')
    datetick('x',6)
    set(gca,'fontsize',14)
    box on; grid on
    subplot(313); hold on
    plot(t,p1-p1_fit(:,imin1),'sk','markersize',10,'linewidth',2)
    ylabel('residuals (cm)')
    datetick('x',6)
    set(gca,'fontsize',14)
    box on; grid on

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print(['../figures/cal_fits/fitting/' POBS_dir(k).name(1:end-4) '_G1'],'-dpng','-r100')

    figure(14); clf
    subplot(311); hold on
    plot(linspace(1,180,1000),stds2,'.','markersize',5)
    ylabel('RMS (cm)')
    xlabel('time constant (d)')
    title([POBS_dir(k).name(1:end-4) ' Gauge 2'])
    set(gca,'fontsize',14)
    box on; grid on
    subplot(312); hold on
    plot(t,p2,'o','markersize',10,'linewidth',2)
    plot(t,p2_fit(:,imin2),'linewidth',1)
    legend('cals','fit','location','northeast')
    ylabel('P (cm)')
    datetick('x',6)
    set(gca,'fontsize',14)
    box on; grid on
    subplot(313); hold on
    plot(t,p2-p2_fit(:,imin2),'sk','markersize',10,'linewidth',2)
    ylabel('residuals (cm)')
    datetick('x',6)
    set(gca,'fontsize',14)
    box on; grid on

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print(['../figures/cal_fits/fitting/' POBS_dir(k).name(1:end-4) '_G2'],'-dpng','-r100')

    % intervene as necessary
    % keyboard

    % add fit info into structure
    calInfoAll1.m=[m1(:,imin1);lamda_list(imin1)];
    calInfoAll2.m=[m2(:,imin2);lamda_list(imin2)];

    %% apply modeled drift correction to the pressure records

    % trim as necesary
    if i==1 % offset over last 10 cals
        icut=dataf.tf>=calInfoAll1.t0p(end-10);
        dataf.tf(icut)=[]; dataf.p1f(icut)=[]; dataf.T1f(icut)=[];
        dataf.p2f(icut)=[]; dataf.T2f(icut)=[];
    elseif i==4
        icut=dataf.tf>=calInfoAll1.t0p(end-9);
        dataf.tf(icut)=[]; dataf.p1f(icut)=[]; dataf.T1f(icut)=[];
        dataf.p2f(icut)=[]; dataf.T2f(icut)=[];
    end

    % correct full pressure time series with drift model
    tinv=(dataf.tf-dataf.tf(1))/max(dataf.tf-dataf.tf(1));
    T1inv=interp1(t,T1,dataf.tf,'spline','extrap')-median(T1);
    T2inv=interp1(t,T2,dataf.tf,'spline','extrap')-median(T2);
    p1_m=tinv*calInfoAll1.m(1)+T1inv*calInfoAll1.m(2)+calInfoAll1.m(3)+calInfoAll1.m(4)*...
        exp(-tinv/calInfoAll1.m(5));
    p2_m=tinv*calInfoAll2.m(1)+T2inv*calInfoAll2.m(2)+calInfoAll2.m(3)+calInfoAll2.m(4)*...
        exp(-tinv/calInfoAll2.m(5));

    % store corrected time series
    dataf.p1_mcor=dataf.p1f-p1_m;
    dataf.p2_mcor=dataf.p2f-p2_m;
    if i==5
        inan=find(dataf.tf>=739143); % determined from visual inspection
        p1_m(inan)=NaN;
        dataf.p1f(inan)=NaN;
        dataf.p1_mcor(inan)=NaN;
    end

    % plots to demonstrate correction
    figure(33); clf;
    subplot(211); hold on
    plot(dataf.tf,dataf.p1f-mean(dataf.p1f)+10,'linewidth',1)
    plot(dataf.tf,p1_m-mean(p1_m)+10,'linewidth',1)
    p1_plot=dataf.p1_mcor-mean(dataf.p1_mcor);
    plot(dataf.tf,p1_plot,'k','linewidth',1)
    p1_fit=polyfit(dataf.tf,p1_plot,1);
    text(dataf.tf(end-1000),p1_plot(end)-7.5,...
        [num2str(round(p1_fit(1)*365,1)) ' hPa/yr'],'fontsize',14)
    if i==5
        plot(dataf.tf,dataf.p1f-nanmean(dataf.p1f)+10,'linewidth',1)
        plot(dataf.tf,p1_m-nanmean(p1_m)+10,'r','linewidth',1)
        p1_plot=dataf.p1_mcor-nanmean(dataf.p1_mcor);
        plot(dataf.tf,p1_plot,'k','linewidth',1)
        p1_fit=polyfit(dataf.tf(~isnan(dataf.p1_mcor)),p1_plot(~isnan(dataf.p1_mcor)),1);
        text(dataf.tf(end-1000),p1_plot(inan(1)-1)-7.5,...
            [num2str(round(p1_fit(1)*365,1)) ' hPa/yr'],'fontsize',14)
    end
    ylabel('P (hPa)')
    legend('P','cal fit','corrected','location','northeast')
    datetick('x',6)
    title({POBS_dir(k).name(1:end-4); 'Gauge 1'})
    set(gca,'fontsize',14)
    box on; grid on
    lim1=ylim;
    subplot(212); hold on
    plot(dataf.tf,dataf.p2f-mean(dataf.p2f)+10,'linewidth',1)
    plot(dataf.tf,p2_m-mean(p2_m)+10,'linewidth',1)
    p2_plot=dataf.p2_mcor-mean(dataf.p2_mcor);
    plot(dataf.tf,p2_plot,'k','linewidth',1)
    p2_fit=polyfit(dataf.tf,p2_plot,1);
    text(dataf.tf(end-1000),p2_plot(end)-7.5,...
        [num2str(round(p2_fit(1)*365,1)) ' hPa/yr'],'fontsize',14)
    ylabel('P (hPa)')
    legend('P','cal fit','corrected','location','northeast')
    datetick('x',6)
    title('Gauge 2')
    set(gca,'fontsize',14)
    box on; grid on
    lim2=ylim;

    % equilibrate ylimits
    subplot(211)
    ylim([min(lim1(1),lim2(1)) max(lim1(2),lim2(2))])
    subplot(212)
    ylim([min(lim1(1),lim2(1)) max(lim1(2),lim2(2))])

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print(['../figures/cal_fits/corrections/model/' POBS_dir(k).name(1:end-4) '_driftcor'],'-dpng','-r300')

    p1m_rate(i)=p1_fit(1)*365;
    p2m_rate(i)=p2_fit(1)*365;

    figure(17); clf
    subplot(211); hold on
    plot(dataf.tf,p1_plot,'linewidth',1)
    plot(dataf.tf,p2_plot,'linewidth',1)
    ylabel('P (hPa)')
    legend('Gauge 1','Gauge 2','location','northeast')
    datetick('x',6)
    title([POBS_dir(k).name(1:end-4) ' Gauge Comparisons -- model drift correction'])
    set(gca,'fontsize',14)
    box on; grid on
    subplot(212); hold on
    plot(dataf.tf,p1_plot-p2_plot,'k','linewidth',1)
    ylabel('P (hPa)')
    legend('Gauge 1-Gauge 2','location','northeast')
    datetick('x',6)
    set(gca,'fontsize',14)
    box on; grid on

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print(['../figures/cal_fits/corrections/model/gauge_comp/' POBS_dir(k).name(1:end-4)  '_gaugecomp'],'-dpng','-r300')

    %% apply direct drift correction to the pressure records

    % expand calibrations onto full time domain of pressure data
    p1_m=interp1(t,p1,dataf.tf,'linear','extrap');
    p2_m=interp1(t,p2,dataf.tf,'linear','extrap');
    % spline fit does better job of capturing initial exponent
    p1_ms=interp1(t,p1,dataf.tf,'spline','extrap');
    p2_ms=interp1(t,p2,dataf.tf,'spline','extrap');
    ispl=find(dataf.tf<=t(end));
    p1_m(ispl)=p1_ms(ispl);
    p2_m(ispl)=p2_ms(ispl);
    % apply temperature correction, from above model
    p1_m=p1_m-T1inv*calInfoAll1.m(2);
    p2_m=p2_m-T2inv*calInfoAll2.m(2);

    % store corrected time series
    dataf.p1_dcor=dataf.p1f-p1_m;
    dataf.p2_dcor=dataf.p2f-p2_m;
    if i==5
        inan=find(dataf.tf>=739143); % determined from visual inspection
        p1_m(inan)=NaN;
        dataf.p1f(inan)=NaN;
        dataf.p1_dcor(inan)=NaN;
        p1(end-5:end)=NaN;
    end

    figure(53); clf;
    subplot(211); hold on
    plot(dataf.tf,dataf.p1f-mean(dataf.p1f)+10,'linewidth',1)
    plot(t,p1-mean(p1)+10,'or')
    plot(dataf.tf,p1_m-mean(p1)+10,'r','linewidth',1)
    p1_plot=dataf.p1_dcor-mean(dataf.p1_dcor);
    plot(dataf.tf,p1_plot,'k','linewidth',1)
    % determine linear fit to corrected data
    p1_fit=polyfit(dataf.tf,p1_plot,1);
    text(dataf.tf(end-1000),p1_plot(end)-7.5,...
        [num2str(round(p1_fit(1)*365,1)) ' hPa/yr'],'fontsize',14)
    if i==5
        clf; subplot(211); hold on
        plot(dataf.tf,dataf.p1f-nanmean(dataf.p1f)+10,'linewidth',1)
        plot(t,p1-nanmean(p1)+10,'or')
        plot(dataf.tf,p1_m-nanmean(p1)+10,'r','linewidth',1)
        p1_plot=dataf.p1_dcor-nanmean(dataf.p1_dcor);
        plot(dataf.tf,p1_plot,'k','linewidth',1)
        p1_fit=polyfit(dataf.tf(~isnan(dataf.p1_dcor)),p1_plot(~isnan(dataf.p1_dcor)),1);
        text(dataf.tf(end-1000),p1_plot(inan(1)-1)-7.5,...
            [num2str(round(p1_fit(1)*365,1)) ' hPa/yr'],'fontsize',14)
    end
    ylabel('P (hPa)')
    legend('P','cals','cal extrap','corrected','location','northeast')
    datetick('x',6)
    title({POBS_dir(k).name(1:end-4); 'Gauge 1'})
    set(gca,'fontsize',14)
    box on; grid on
    lim1=ylim;
    subplot(212); hold on
    plot(dataf.tf,dataf.p2f-mean(dataf.p2f)+10,'linewidth',1)
    plot(t,p2-mean(p2)+10,'or')
    plot(dataf.tf,p2_m-mean(p2)+10,'r','linewidth',1)
    p2_plot=dataf.p2_dcor-mean(dataf.p2_dcor);
    plot(dataf.tf,p2_plot,'k','linewidth',1)
    % determine linear fit to corrected data
    p2_fit=polyfit(dataf.tf,p2_plot,1);
    text(dataf.tf(end-1000),p2_plot(end)-7.5,...
        [num2str(round(p2_fit(1)*365,1)) ' hPa/yr'],'fontsize',14)
    ylabel('P (hPa)')
    legend('P','cals','cal extrap','corrected','location','northeast')
    datetick('x',6)
    title('Gauge 2')
    set(gca,'fontsize',14)
    box on; grid on
    lim2=ylim;

    % equilibrate ylimits
    subplot(211)
    ylim([min(lim1(1),lim2(1)) max(lim1(2),lim2(2))])
    subplot(212)
    ylim([min(lim1(1),lim2(1)) max(lim1(2),lim2(2))])

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print(['../figures/cal_fits/corrections/direct/' POBS_dir(k).name(1:end-4) '_driftcor'],'-dpng','-r300')

    % uncorrected linear rates
    p1m_fit=polyfit(dataf.tf,dataf.p1f,1);
    p2m_fit=polyfit(dataf.tf,dataf.p2f,1);
    p1m_rate(i)=p1m_fit(1)*365;
    p2m_rate(i)=p2m_fit(1)*365;

    % corrected linear rates
    p1_rate(i)=p1_fit(1)*365;
    p2_rate(i)=p2_fit(1)*365;

    figure(18); clf
    subplot(211); hold on
    plot(dataf.tf,p1_plot,'linewidth',1)
    plot(dataf.tf,p2_plot,'linewidth',1)
    ylabel('P (hPa)')
    legend('Gauge 1','Gauge 2','location','northeast')
    datetick('x',6)
    title([POBS_dir(k).name(1:end-4) ' Gauge Comparisons -- direct drift correction'])
    set(gca,'fontsize',14)
    box on; grid on
    subplot(212); hold on
    plot(dataf.tf,p1_plot-p2_plot,'k','linewidth',1)
    ylabel('P (hPa)')
    legend('Gauge 1-Gauge 2','location','northeast')
    datetick('x',6)
    set(gca,'fontsize',14)
    box on; grid on

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print(['../figures/cal_fits/corrections/direct/gauge_comp/' POBS_dir(k).name(1:end-4) '_gaugecomp'],'-dpng','-r300')

    % stacked timeseries of (un)corrected data to demonstrate difference
    figure(111); hold on
    plot(dataf.tf,dataf.p2f-mean(dataf.p2f)+10*n,'color',[0 114 189]/255,'linewidth',1)
    plot(dataf.tf,p2m_fit(1)*dataf.tf+p2m_fit(2)-mean(dataf.p2f)+10*n,'color',[0 114 189]/255,'linewidth',1)
    text(dataf.tf(1)-55,p2m_fit(1)*dataf.tf(1)+p2m_fit(2)-mean(dataf.p2f)+10*n,POBS_dir(k).name(1:end-4),'color',[0 114 189]/255,'fontsize',12)

    figure(112); hold on
    plot(dataf.tf,p2_plot+10*n,'r','linewidth',1)
    plot(dataf.tf,p2_fit(1)*dataf.tf+p2_fit(2)+10*n,'r','linewidth',1)
    text(dataf.tf(1)-55,p2_fit(1)*dataf.tf(1)+p2_fit(2)+10*n,POBS_dir(k).name(1:end-4),'color','r','fontsize',12)

    % save corrected data
    save(['../stitched_data/drift_corrected/' POBS_dir(k).name],'dataf','calInfoAll1','calInfoAll2')

    n=n+1;
end

%% MAPS OF INFERED DEFORMATION

addpath('../../hikurangi/code/m_map')
load('../pressure_data/geometry');

% reconcile different station order
stalist = {POBS_dir(i_list).name};
for i=1:length(stalist)
    staID = str2double(stalist{(i)}(5:end-4));
    staID_padded = sprintf('%02d',staID);
    ii(i) = find(strcmp(['POBS-' staID_padded],staname));
end

% setup base map
figure(54); clf; hold on
lat1=-40; lat2=-38.5;
lon1=177.5; lon2=179.5;
m_proj('mercator','longitudes',[lon1 lon2],'latitudes',[lat1 lat2]);

X=linspace(174.5,180,1320);
Y=linspace(-37.5,-42,1080);
Z=imread('../../hikurangi/gshhg-bin-2.3.6/NZ_bathymetry.tiff');
Z(Z>=0)=NaN;
[~,hc]=m_contourf(X,Y,Z,-5000:100:0,'linecolor','none');

m_gshhs_i('patch',[.7 .7 .7]);

m_contour(X,Y,Z,[-150 -150],'k','linewidth',1) % shelf at 150 m?
m_text(177.75,-39.55,'150 m','fontsize',12,'rotation',35)
m_contour(X,Y,Z,[-2750 -2750],'k','linewidth',1) % 2750 m does a decent job visually
m_text(178.4,-39.9,'2750 m','fontsize',12,'rotation',80)

colormap(gca,cmocean('-deep'))
cb1 = colorbar(gca,'eastoutside');
ylabel(cb1,'Depth (m)')
set(cb1,'fontsize',14)
m_grid('xlabeldir','end','fontsize',14);

% station markers
m_plot(stalon,stalat,'^k','markersize',10)
m_text(stalon-0.21,stalat,staname,'fontsize',12)

% scale arrow
m_quiver(177.6,-39.95,0,5/20,'off','color','k','linewidth',2)
m_text(177.65,-39.8,'5 cm/yr','fontsize',12)

%---direct corrections
% Gauge 1
h1=m_quiver(stalon(ii)',stalat(ii)',zeros(1,10),-p1_rate/20,'off','color','r','linewidth',2);
title('Gauge 1 -- direct drift correction')

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../figures/cal_fits/corrections/direct/def_map_G1','-dpng','-r300')
print('../figures/cal_fits/corrections/direct/def_map_G1','-depsc','-vector')

% Gauge 2
delete(h1)
h1=m_quiver(stalon(ii)',stalat(ii)',zeros(1,10),-p2_rate/20,'off','color','r','linewidth',2);
title('Gauge 2 -- direct drift correction')

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../figures/cal_fits/corrections/direct/def_map_G2','-dpng','-r300')
print('../figures/cal_fits/corrections/direct/def_map_G2','-depsc','-vector')

%---smooth model corrections
% Gauge 1
delete(h1)
h1=m_quiver(stalon(ii)',stalat(ii)',zeros(1,10),-p1m_rate/20,'off','color','r','linewidth',2);
title('Gauge 1 -- model drift correction')

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../figures/cal_fits/corrections/model/def_map_G1','-dpng','-r300')
print('../figures/cal_fits/corrections/model/def_map_G1','-depsc','-vector')

% Gauge 2
delete(h1)
h1=m_quiver(stalon(ii)',stalat(ii)',zeros(1,10),-p2m_rate/20,'off','color','r','linewidth',2);
title('Gauge 2 -- model drift correction')

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../figures/cal_fits/corrections/model/def_map_G2','-dpng','-r300')
print('../figures/cal_fits/corrections/model/def_map_G2','-depsc','-vector')