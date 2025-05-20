% model_filter_test.m
%
% Compares A-0-A drift-corrected pressure data against oceanographic
% models at various periods
%

clear; close all

% load HYCOM pressures
load('../../hikurangi/ocean_data/HYCOM/hycom_pressure.mat')
% load GLORYS pressures
load('../../hikurangi/ocean_data/GLORYS/glorys_pressure.mat')
% load satellite altimetry
load('../../hikurangi/ocean_data/SSH/satellite_ssh.mat')
% load ECCO2 pressures
load('../../hikurangi/ocean_data/ECCO2/ecco2_pressure_update.mat')

POBS_dir=dir('../stitched_data/drift_corrected');
POBS_list={POBS_dir.name}';
file_check=cellfun(@(v)v(1),POBS_list);
i_list=find(eq(file_check,'P')); % subselects only files begining with 'P'

% temporarily exclude POBS09 & POBS10
i_list([1,3])=[];

h_rate=nan(10,1);
s_rate=nan(10,1);
g_rate=nan(10,1);
a_rate=nan(10,1);
e_rate=nan(10,1);
for i=i_list'
    load([POBS_dir(i).folder '/' POBS_dir(i).name])

    if strcmp(POBS_dir(i).name,'POBS12.mat') % gauge 1 data were truncated with NaNs, which will cause problems later
        inan=find(isnan(dataf.p1_dcor));
        dataf.p1_dcor(inan)=dataf.p2_dcor(inan)-(dataf.p2_dcor(inan(1))-dataf.p1_dcor(inan(1)-1));
    end

    % extract corrected A-0-A time series
    t1=dataf.tf; t2=t1;
    p1=dataf.p1_dcor; p1=p1-mean(p1);
    p2=dataf.p2_dcor; p2=p2-mean(p2);
    
    %% HYCOM
    % find model index that corresponds to current station
    shortname=POBS_dir(i).name(1:end-4);
    if length(shortname)==5
        shortname=[shortname(1:4) '-0' shortname(5)];
    elseif length(shortname)==6
        shortname=[shortname(1:4) '-' shortname(5:6)];
    else
        shortname=[shortname(1:4) '-0' shortname(5)];
    end
    i_hy=find(strcmp(shortname,p_hy.nearest_sta));

    % trim HYCOM data to match APG
    th=p_hy.t{i_hy};
    ph=p_hy.bpa{i_hy};
    ph=ph(th>=dataf.tf(1) & th<=dataf.tf(end));
    th=th(th>=dataf.tf(1) & th<=dataf.tf(end));

    % linear fits to the time series
    pm=polyfit(th-th(1),ph,1);
    h_rate(i)=pm(1)*365; % store HYCOM rate
    ph_l=polyval(pm,th-th(1));

    % figure(45); clf; hold on
    % subplot(211); hold on
    % plot(dataf.tf,dataf.p1_dcor-mean(dataf.p1_dcor),'linewidth',1)
    % plot(th,ph,'linewidth',1)
    % plot(dataf.tf,p1_l,'color',[0 114 189]/255,'linewidth',1)
    % plot(th,ph_l,'r','linewidth',1)
    % ylabel('P (hPa)')
    % legend('APG','HYCOM','location','northwest')
    % datetick('x',6)
    % title({shortname; 'Gauge 1'})
    % set(gca,'fontsize',14)
    % box on; grid on
    % subplot(212); hold on
    % plot(dataf.tf,dataf.p2_dcor-mean(dataf.p2_dcor),'linewidth',1)
    % plot(th,ph,'linewidth',1)
    % plot(dataf.tf,p2_l,'color',[0 114 189]/255,'linewidth',1)
    % plot(th,ph_l,'r','linewidth',1)
    % ylabel('P (hPa)')
    % legend('APG','HYCOM','location','northwest')
    % datetick('x',6)
    % title('Gauge 2')
    % set(gca,'fontsize',14)
    % box on; grid on
    % 
    % fh=gcf;
    % fh.PaperUnits='inches';
    % fh.PaperPosition=[0 0 8.5 11];
    % print(['../figures/trend_comps/hycom/hycom_' POBS_dir(i).name(1:end-4)],'-dpng','-r300')

    %% SATELLITE SSH
    % find grid index that corresponds to current station
    i_ssh=find(strcmp(shortname,sat_ssh.nearest_sta));

    % trim satellite data to match APG
    ts=sat_ssh.t{i_ssh};
    ps=sat_ssh.sla_int{i_ssh};
    ps=ps(ts>=dataf.tf(1) & ts<=dataf.tf(end));
    ts=ts(ts>=dataf.tf(1) & ts<=dataf.tf(end));

    % linear fits to the time series
    pm=polyfit(ts-ts(1),ps,1);
    s_rate(i)=pm(1)*365; % store satellite rate
    ps_l=polyval(pm,ts-ts(1));

    % figure(46); clf; hold on
    % subplot(211); hold on
    % plot(dataf.tf,dataf.p1_dcor-mean(dataf.p1_dcor),'linewidth',1)
    % plot(ts,ps,'linewidth',1)
    % plot(dataf.tf,p1_l,'color',[0 114 189]/255,'linewidth',1)
    % plot(ts,ps_l,'r','linewidth',1)
    % ylabel('P (hPa)')
    % legend('APG','SSH','location','northwest')
    % datetick('x',6)
    % title({shortname; 'Gauge 1'})
    % set(gca,'fontsize',14)
    % box on; grid on
    % subplot(212); hold on
    % plot(dataf.tf,dataf.p2_dcor-mean(dataf.p2_dcor),'linewidth',1)
    % plot(ts,ps,'linewidth',1)
    % plot(dataf.tf,p2_l,'color',[0 114 189]/255,'linewidth',1)
    % plot(ts,ps_l,'r','linewidth',1)
    % ylabel('P (hPa)')
    % legend('APG','SSH','location','northwest')
    % datetick('x',6)
    % title('Gauge 2')
    % set(gca,'fontsize',14)
    % box on; grid on
    % 
    % fh=gcf;
    % fh.PaperUnits='inches';
    % fh.PaperPosition=[0 0 8.5 11];
    % print(['../figures/trend_comps/satellite/satellite_' POBS_dir(i).name(1:end-4)],'-dpng','-r300')

    %% GLORYS
    % find model index that corresponds to current station
    i_glo=find(strcmp(shortname,p_glo.nearest_sta));

    % trim GLORYS data to match APG
    tg=p_glo.t{i_glo};
    pg=p_glo.bpa{i_glo};
    pg=pg(tg>=dataf.tf(1) & tg<=dataf.tf(end));
    tg=tg(tg>=dataf.tf(1) & tg<=dataf.tf(end));

    % linear fits to the time series
    pm=polyfit(tg-tg(1),pg,1);
    g_rate(i)=pm(1)*365; % store GLORYS rate
    pg_l=polyval(pm,tg-tg(1));
    pm=polyfit(t1-t1(1),p1,1);
    p1_l=polyval(pm,t1-t1(1));
    pm=polyfit(t2-t2(1),p2,1);
    a_rate(i)=pm(1)*365; % store A-0-A rate
    p2_l=polyval(pm,t2-t2(1));

    % figure(47); clf; hold on
    % subplot(211); hold on
    % plot(t1,p1,'linewidth',1)
    % plot(tg,pg,'linewidth',1)
    % plot(t1,p1_l,'color',[0 114 189]/255,'linewidth',1)
    % plot(tg,pg_l,'r','linewidth',1)
    % ylabel('P (hPa)')
    % legend('APG','GLORYS','location','northwest')
    % datetick('x',6)
    % title({shortname; 'Gauge 1'})
    % set(gca,'fontsize',14)
    % box on; grid on
    % subplot(212); hold on
    % plot(t2,p2,'linewidth',1)
    % plot(tg,pg,'linewidth',1)
    % plot(t2,p2_l,'color',[0 114 189]/255,'linewidth',1)
    % plot(tg,pg_l,'r','linewidth',1)
    % ylabel('P (hPa)')
    % legend('APG','GLORYS','location','northwest')
    % datetick('x',6)
    % title('Gauge 2')
    % set(gca,'fontsize',14)
    % box on; grid on
    % 
    % fh=gcf;
    % fh.PaperUnits='inches';
    % fh.PaperPosition=[0 0 8.5 11];
    % print(['../figures/trend_comps/glorys/glorys_' POBS_dir(i).name(1:end-4)],'-dpng','-r300')

    % interpolate GLORYS onto APG time basis
    pgi=interp1(tg,pg,t1,'linear','extrap');

    % filter parameters
    Tc=[30]; % cutoff periods in days

    for j=1:length(Tc)

        % low-pass filter data and model
        [b,a]=butter(4,(2/(Tc(j)*24)),'low');
        p1f=filtfilt(b,a,p1);
        p2f=filtfilt(b,a,p2);
        pmf=filtfilt(b,a,pgi);

        % check1=std(detrend(p1f)); check1=round(check1,1);
        % check2=std(detrend(pmf)); check2=round(check2,1);
        % check3=std(detrend(p1f-pmf)); check3=round(check3,1);

        sfc1=inv(detrend(pmf)'*detrend(pmf))*detrend(pmf)'*detrend(p1f);
        % sfc1=1;
        sfc2=inv(detrend(pmf)'*detrend(pmf))*detrend(pmf)'*detrend(p2f);
        % sfc2=1;

        % linear fits to the filtered time series
        pm=polyfit(t1-t1(1),pmf,1);
        pg_l=polyval(pm,t1-t1(1));
        gf_rate=pm(1)*365; % [hPa/yr]
        pm=polyfit(t1-t1(1),p1f,1);
        p1_l=polyval(pm,t1-t1(1));
        p1f_rate=pm(1)*365; % [hPa/yr]
        pm=polyfit(t2-t2(1),p2f,1);
        p2_l=polyval(pm,t2-t2(1));
        p2f_rate=pm(1)*365; % [hPa/yr]

        figure(48); clf; hold on
        subplot(211); hold on
        plot(t1,p1f,'linewidth',1)
        plot(t1,sfc1*pmf,'linewidth',1)
        plot(t1,p1_l,'color',[0 114 189]/255,'linewidth',1)
        plot(t1,sfc1*pg_l,'r','linewidth',1)
        plot(t1,p1_l-sfc1*pg_l,'k','linewidth',1)
        text(t1(end-500),1,num2str(p1f_rate),'color',[0 114 189]/255)
        text(t1(end-500),0,num2str(sfc1*gf_rate),'color','r')
        text(t1(end-500),-1,num2str(p1f_rate-sfc1*gf_rate),'color','k')
        ylabel('P (hPa)')
        legend('APG','GLORYS','location','northwest')
        datetick('x',6)
        title({[shortname ' [' num2str(Tc(j)) ' day filter]']; 'Gauge 1'})
        set(gca,'fontsize',14)
        box on; grid on
        subplot(212); hold on
        plot(t2,p2f,'linewidth',1)
        plot(t2,sfc2*pmf,'linewidth',1)
        plot(t2,p2_l,'color',[0 114 189]/255,'linewidth',1)
        plot(t2,sfc2*pg_l,'r','linewidth',1)
        plot(t2,p2_l-sfc2*pg_l,'k','linewidth',1)
        text(t2(end-500),1,num2str(p2f_rate),'color',[0 114 189]/255)
        text(t2(end-500),0,num2str(sfc2*gf_rate),'color','r')
        text(t2(end-500),-1,num2str(p2f_rate-sfc2*gf_rate),'color','k')
        ylabel('P (hPa)')
        legend('APG','GLORYS','location','northwest')
        datetick('x',6)
        title('Gauge 2')
        set(gca,'fontsize',14)
        box on; grid on

        fh=gcf;
        fh.PaperUnits='inches';
        fh.PaperPosition=[0 0 8.5 11];
        print(['../figures/trend_comps/glorys/' num2str(Tc(j)) '_day_filter/scaled_glorys_' POBS_dir(i).name(1:end-4)],'-dpng','-r300')

    end

    %% ECCO2
    % find model index that corresponds to current station
    i_eco=find(strcmp(shortname,ecco2.nearest_sta));

    % trim ECCO data to match APG
    te=ecco2.t{i_eco};
    pe=ecco2.bpa_int{i_eco};
    pe=pe(te>=(dataf.tf(1)) & te<=dataf.tf(end));
    te=te(te>=(dataf.tf(1)) & te<=dataf.tf(end));

    % linear fits to the time series
    pm=polyfit(te-te(1),pe,1);
    e_rate(i)=pm(1)*365; % store ECCO2 rate
    pe_l=polyval(pm,te-te(1));

    % figure(49); clf; hold on
    % subplot(211); hold on
    % plot(dataf.tf,dataf.p1_dcor-mean(dataf.p1_dcor),'linewidth',1)
    % plot(te,pe,'linewidth',1)
    % plot(dataf.tf,p1_l,'color',[0 114 189]/255,'linewidth',1)
    % plot(te,pe_l,'r','linewidth',1)
    % ylabel('P (hPa)')
    % legend('APG','ECCO2','location','northwest')
    % datetick('x',6)
    % title({shortname; 'Gauge 1'})
    % set(gca,'fontsize',14)
    % box on; grid on
    % subplot(212); hold on
    % plot(dataf.tf,dataf.p2_dcor-mean(dataf.p2_dcor),'linewidth',1)
    % plot(te,pe,'linewidth',1)
    % plot(dataf.tf,p2_l,'color',[0 114 189]/255,'linewidth',1)
    % plot(te,pe_l,'r','linewidth',1)
    % ylabel('P (hPa)')
    % legend('APG','ECCO2','location','northwest')
    % datetick('x',6)
    % title('Gauge 2')
    % set(gca,'fontsize',14)
    % box on; grid on
    % 
    % fh=gcf;
    % fh.PaperUnits='inches';
    % fh.PaperPosition=[0 0 8.5 11];
    % print(['../figures/trend_comps/ecco2/ecco2_' POBS_dir(i).name(1:end-4)],'-dpng','-r300')

    % interpolate ECCO onto APG time basis
    pei=interp1(te,pe,t1,'linear','extrap');

    % filter parameters
    Tc=[30]; % cutoff periods in days

    for j=1:length(Tc)

        % low-pass filter data and model
        [b,a]=butter(4,(2/(Tc(j)*24)),'low');
        p1f=filtfilt(b,a,p1);
        p2f=filtfilt(b,a,p2);
        pmf=filtfilt(b,a,pei);

        % check1=std(detrend(p1f)); check1=round(check1,1);
        % check2=std(detrend(pmf)); check2=round(check2,1);
        % check3=std(detrend(p1f-pmf)); check3=round(check3,1);

        sfc1=inv(detrend(pmf)'*detrend(pmf))*detrend(pmf)'*detrend(p1f);
        % sfc1=1;
        sfc2=inv(detrend(pmf)'*detrend(pmf))*detrend(pmf)'*detrend(p2f);
        % sfc2=1;

        % linear fits to the filtered time series
        pm=polyfit(t1-t1(1),pmf,1);
        pe_l=polyval(pm,t1-t1(1));
        ef_rate=pm(1)*365; % [hPa/yr]
        pm=polyfit(t1-t1(1),p1f,1);
        p1_l=polyval(pm,t1-t1(1));
        p1f_rate=pm(1)*365; % [hPa/yr]
        pm=polyfit(t2-t2(1),p2f,1);
        p2_l=polyval(pm,t2-t2(1));
        p2f_rate=pm(1)*365; % [hPa/yr]

        figure(48); clf; hold on
        subplot(211); hold on
        plot(t1,p1f,'linewidth',1)
        plot(t1,sfc1*pmf,'linewidth',1)
        plot(t1,p1_l,'color',[0 114 189]/255,'linewidth',1)
        plot(t1,sfc1*pe_l,'r','linewidth',1)
        plot(t1,p1_l-sfc1*pe_l,'k','linewidth',1)
        text(t1(end-500),1,num2str(p1f_rate),'color',[0 114 189]/255)
        text(t1(end-500),0,num2str(sfc1*ef_rate),'color','r')
        text(t1(end-500),-1,num2str(p1f_rate-sfc1*ef_rate),'color','k')
        ylabel('P (hPa)')
        legend('APG','ECCO2','location','northwest')
        datetick('x',6)
        title({[shortname ' [' num2str(Tc(j)) ' day filter]']; 'Gauge 1'})
        set(gca,'fontsize',14)
        box on; grid on
        subplot(212); hold on
        plot(t2,p2f,'linewidth',1)
        plot(t2,sfc2*pmf,'linewidth',1)
        plot(t2,p2_l,'color',[0 114 189]/255,'linewidth',1)
        plot(t2,sfc2*pe_l,'r','linewidth',1)
        plot(t2,p2_l-sfc2*pe_l,'k','linewidth',1)
        text(t2(end-500),1,num2str(p2f_rate),'color',[0 114 189]/255)
        text(t2(end-500),0,num2str(sfc2*ef_rate),'color','r')
        text(t2(end-500),-1,num2str(p2f_rate-sfc2*ef_rate),'color','k')
        ylabel('P (hPa)')
        legend('APG','ECCO2','location','northwest')
        datetick('x',6)
        title('Gauge 2')
        set(gca,'fontsize',14)
        box on; grid on

        fh=gcf;
        fh.PaperUnits='inches';
        fh.PaperPosition=[0 0 8.5 11];
        print(['../figures/trend_comps/ecco2/' num2str(Tc(j)) '_day_filter/scaled_ecco_' POBS_dir(i).name(1:end-4)],'-dpng','-r300')

    end
end
h_rate(1:2)=[];
s_rate(1:2)=[];
g_rate(1:2)=[];
a_rate(1:2)=[];
e_rate(1:2)=[];

%% MAPS OF OCEAN-INDUCED CHANGE

addpath('../../hikurangi/code/m_map')
load('../pressure_data/geometry');

% reconcile different station order
stalist = {POBS_dir.name};
stalist(1:2) =[];
for i=1:length(stalist)
    staID = str2double(stalist{(i)}(5:end-4));
    if isnan(staID)
        staID=7;
    end
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

%---HYCOM
h1=m_quiver(stalon(ii),stalat(ii),zeros(10,1),-h_rate/20,'off','color','r','linewidth',2);
title('HYCOM inferred deformation')

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../figures/trend_comps/hycom/mapview','-dpng','-r300')

%---SATELLITE
delete(h1)
h1=m_quiver(stalon(ii),stalat(ii),zeros(10,1),-s_rate/20,'off','color','r','linewidth',2);
title('Altimetry inferred deformation')

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../figures/trend_comps/satellite/mapview','-dpng','-r300')

%---GLORYS
delete(h1)
h1=m_quiver(stalon(ii),stalat(ii),zeros(10,1),-g_rate/20,'off','color','r','linewidth',2);
title('GLORYS inferred deformation')

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../figures/trend_comps/glorys/mapview','-dpng','-r300')

% this one actually compares well, so let's take a difference
delete(h1)
h1=m_quiver(stalon(ii),stalat(ii),zeros(10,1),-(a_rate-g_rate)/20,'off','color','r','linewidth',2);
title('A0A - GLORYS inferred deformation')

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../figures/trend_comps/glorys/mapview_dif','-dpng','-r300')

%---ECCO2
delete(h1)
h1=m_quiver(stalon(ii),stalat(ii),zeros(10,1),-e_rate/20,'off','color','r','linewidth',2);
title('ECCO2 inferred deformation')

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../figures/trend_comps/ecco2/mapview_alt','-dpng','-r300')