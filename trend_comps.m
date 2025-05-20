% trend_comps.m
%
% Generates stacked plots of pressure data before and after drift
% correction (both gauges), as well as a boxplot (scatter?) that
% demonstrates the change in trend as a result.
%

clear; close all
load('../pressure_data_Y2/geometry');

POBS_dir=dir('../stitched_data_Y2/drift_corrected');
POBS_list={POBS_dir.name}';
file_check=cellfun(@(v)v(1),POBS_list);
i_list=find(eq(file_check,'P')); % subselects only files begining with 'P'

% sort everything by depth
[d,id]=sort(stadepth,'descend');
sname=staname(id);
fname=filename(id);
sla=stalat(id);
slo=stalon(id);
il=i_list(id);

figure(4); clf; hold on
figure(5); clf; hold on
figure(6); clf; hold on
for i=1:length(il)
    k=il(i);
    load(['../stitched_data_Y2/drift_corrected/' POBS_dir(k).name],'dataf');

    ta1=dataf.tf; pa1=dataf.p1f;
    ta1=ta1(~isnan(pa1)); pa1=pa1(~isnan(pa1)); pa1=pa1-mean(pa1);
    ma1=polyfit(ta1,pa1,1);
    fa1=polyval(ma1,ta1);

    ta2=dataf.tf; pa2=dataf.p2f;
    ta2=ta2(~isnan(pa2)); pa2=pa2(~isnan(pa2)); pa2=pa2-mean(pa2);
    ma2=polyfit(ta2,pa2,1);
    fa2=polyval(ma2,ta2);

    tb1=dataf.tf; pb1=dataf.p1_dcor;
    tb1=tb1(~isnan(pb1)); pb1=pb1(~isnan(pb1)); pb1=pb1-mean(pb1);
    mb1=polyfit(tb1,pb1,1);
    fb1=polyval(mb1,tb1);

    tb2=dataf.tf; pb2=dataf.p2_dcor;
    tb2=tb2(~isnan(pb2)); pb2=pb2(~isnan(pb2)); pb2=pb2-mean(pb2);
    mb2=polyfit(tb2,pb2,1);
    fb2=polyval(mb2,tb2);

    figure(4) % Gauge 1, before and after drift correction
    plot(ta1,pa1+(i-1)*10,'r','linewidth',1)
    plot(tb1,pb1+(i-1)*10,'color',[0 114 189]/255,'linewidth',1)
    plot(ta1,fa1+(i-1)*10,'r','linewidth',1)
    plot(tb1,fb1+(i-1)*10,'color',[0 114 189]/255,'linewidth',1)
    if strcmp(sname(i),'POBS-09')
        text(ta1(1)-60,mean([pa1(end-99:end);pb1(end-99:end)])+(i-1)*10,sname(i),'fontsize',12)
    else
        text(ta1(end)+10,mean([pa1(end-99:end);pb1(end-99:end)])+(i-1)*10,sname(i),'fontsize',12)
    end

    figure(5) % Gauge 2, before and after drift correction
    plot(ta2,pa2+(i-1)*10,'r','linewidth',1)
    plot(tb2,pb2+(i-1)*10,'color',[0 114 189]/255,'linewidth',1)
    plot(ta2,fa2+(i-1)*10,'r','linewidth',1)
    plot(tb2,fb2+(i-1)*10,'color',[0 114 189]/255,'linewidth',1)
    if strcmp(sname(i),'POBS-09')
        text(ta2(1)-60,mean([pa2(end-99:end);pb2(end-99:end)])+(i-1)*10,sname(i),'k','fontsize',12)
    else
        text(ta2(end)+10,mean([pa2(end-99:end);pb2(end-99:end)])+(i-1)*10,sname(i),'k','fontsize',12)
    end

    figure(6)
    plot(ma1(1)*365,i-0.1,'ro','markerfacecolor','r','markersize',10)
    plot(mb1(1)*365,i-0.1,'o','color',[0 114 189]/255,'markerfacecolor',[0 114 189]/255,'markersize',10)
    plot(ma2(1)*365,i+0.1,'rs','markerfacecolor','r','markersize',10)
    plot(mb2(1)*365,i+0.1,'s','color',[0 114 189]/255,'markerfacecolor',[0 114 189]/255,'markersize',10)

    a1r(i)=ma1(1)*365;
    b1r(i)=mb1(1)*365;
    a2r(i)=ma2(1)*365;
    b2r(i)=mb2(1)*365;

end

figure(4)
xlim([datenum(2023,11,01) datenum(2025,01,01)])
datetick('x',6,'keeplimits')
ylabel('P (hPa)')
legend('Pre-correction','Post-correction','location','northeast')
title('Gauge 1')
set(gca,'fontsize',14)
ylim([-10 130])
box on; grid on
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print('../figures_Y2/stacks_comps/G1_corrected_pstack','-dpng','-r300')
print('../figures_Y2/stacks_comps/G1_corrected_pstack','-depsc','-vector')

figure(5)
xlim([datenum(2023,11,01) datenum(2025,01,01)])
datetick('x',6,'keeplimits')
ylabel('P (hPa)')
legend('Pre-correction','Post-correction','location','northeast')
title('Gauge 2')
set(gca,'fontsize',14)
ylim([-10 130])
box on; grid on
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print('../figures_Y2/stacks_comps/G2_corrected_pstack','-dpng','-r300')
print('../figures_Y2/stacks_comps/G2_corrected_pstack','-depsc','-vector')

figure(7); clf; hold on
stem((1:12)-0.1,a1r,'ro','markersize',10,'markerfacecolor','r')
stem((1:12)-0.1,b1r','o','markersize',10,'color',[0 114 189]/255,'markerfacecolor',[0 114 189]/255)
stem((1:12)+0.1,a2r,'rs','markersize',10,'markerfacecolor','r')
stem((1:12)+0.1,b2r','s','markersize',10,'color',[0 114 189]/255,'markerfacecolor',[0 114 189]/255)
% view(90,90)
ylabel('Slope (hPa/yr)')
legend('G1 pre','G1 post','G2 pre','G2 post','location','southeast')
set(gca,'fontsize',14)
set(gca,'xtick',1:12)
namelist={POBS_dir(il).name};
for k=1:length(namelist)
    namelist{k}=namelist{k}(1:end-4);
end
set(gca,'xticklabels',namelist)
xlim([0 13])
box on; grid on
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../figures_Y2/stacks_comps/trend_summary','-dpng','-r300')
print('../figures_Y2/stacks_comps/trend_summary','-depsc','-vector')

%% COMPARE AGAINST ECCO2

% load('../../hikurangi/ocean_data/ECCO2/ecco2_pressure_update.mat')
% 
% a1_rate=nan(12,1);
% a2_rate=nan(12,1);
% e1_rate=nan(12,1);
% 
% figure(82); clf; hold on
% figure(83); clf; hold on
% n=0;
% for i=iii
%     k=i_list(i);
%     load(['../stitched_data/drift_corrected/' POBS_dir(k).name],'dataf');
% 
%     % identifier for station
%     shortname=POBS_dir(k).name(1:end-4);
%     if length(shortname)==5
%         shortname=[shortname(1:4) '-0' shortname(5)];
%     elseif length(shortname)==6
%         shortname=[shortname(1:4) '-' shortname(5:6)];
%     else
%         shortname=[shortname(1:4) '-0' shortname(5)];
%     end
% 
%     % gauge 1
%     tb1=dataf.tf; pb1=dataf.p1_dcor;
%     tb1=tb1(~isnan(pb1)); pb1=pb1(~isnan(pb1)); pb1=pb1-mean(pb1);
%     mb1=polyfit(tb1,pb1,1);
%     fb1=polyval(mb1,tb1);
%     a1_rate(i)=mb1(1)*365; % store A-0-A rate
% 
%     % gauge 2
%     tb2=dataf.tf; pb2=dataf.p2_dcor;
%     tb2=tb2(~isnan(pb2)); pb2=pb2(~isnan(pb2)); pb2=pb2-mean(pb2);
%     mb2=polyfit(tb2,pb2,1);
%     fb2=polyval(mb2,tb2);
%     a2_rate(i)=mb2(1)*365; % store A-0-A rate
% 
%     % find model index that corresponds to current station
%     i_eco=find(strcmp(shortname,ecco2.nearest_sta));
% 
%     % trim ECCO data to match APG
%     te1=ecco2.t{i_eco}; pe1=ecco2.bpa_int{i_eco};
%     pe1=pe1(te1>=tb1(1) & te1<=tb1(end));
%     te1=te1(te1>=tb1(1) & te1<=tb1(end));
%     pe1=interp1(te1,pe1,tb1,'linear','extrap');
% 
%     te2=ecco2.t{i_eco}; pe2=ecco2.bpa_int{i_eco};
%     pe2=pe2(te2>=tb2(1) & te2<=tb2(end));
%     te2=te2(te2>=tb2(1) & te2<=tb2(end));
%     pe2=interp1(te2,pe2,tb2,'linear','extrap');
% 
%     % linear fits to the time series
%     pm=polyfit(tb1-tb1(1),pe1,1);
%     e1_rate(i)=pm(1)*365; % store ECCO2 rate
%     pe1_l=polyval(pm,tb1-tb1(1));
%     pm=polyfit(tb2-tb2(1),pe2,1);
%     e2_rate(i)=pm(1)*365; % store ECCO2 rate
%     pe2_l=polyval(pm,tb2-tb2(1));
% 
%     figure(82);
%     plot(tb1,pb1+10*n,'color',[0 114 189]/255,'linewidth',1)
%     plot(tb1,pe1+10*n,'color',[0 100 0]/255,'linewidth',1)
%     % plot(tb1,pb1-pe1+10*n,'k','linewidth',1)
%     plot(tb1,fb1+10*n,'color',[0 114 189]/255,'linewidth',1)
%     plot(tb1,pe1_l+10*n,'color',[0 100 0]/255,'linewidth',1)
%     % plot(tb1,fb1-pe1_l+10*n,'k','linewidth',1)
%     text(tb1(end)+15,mean([pb1(end-99:end);pe1(end-99:end)])+n*10,sname(n+1),'k','fontsize',12)
% 
%     figure(83);
%     plot(tb2,pb2+10*n,'color',[0 114 189]/255,'linewidth',1)
%     plot(tb2,pe2+10*n,'color',[0 100 0]/255,'linewidth',1)
%     % plot(tb2,pb2-pe2+10*n,'k','linewidth',1)
%     plot(tb2,fb2+10*n,'color',[0 114 189]/255,'linewidth',1)
%     plot(tb2,pe2_l+10*n,'color',[0 100 0]/255,'linewidth',1)
%     % plot(tb2,fb2-pe2_l+10*n,'k','linewidth',1)
%     text(tb2(end)+15,mean([pb2(end-99:end);pe2(end-99:end)])+n*10,sname(n+1),'k','fontsize',12)
% 
%     n=n+1;
% end
% 
% figure(82)
% title('Gauge 1')
% ylabel('P (hPa)')
% legend('APG','ECCO2','location','northeast')
% datetick('x',6)
% set(gca,'fontsize',14)
% ylim([-20 110])
% box on; grid on
% fh=gcf;
% fh.PaperUnits='inches';
% fh.PaperPosition=[0 0 8.5 11];
% print('../figures/manuscript/ecco_comparison_G1','-dpng','-r300')
% 
% figure(83)
% title('Gauge 2')
% ylabel('P (hPa)')
% legend('APG','ECCO2','location','northeast')
% datetick('x',6)
% set(gca,'fontsize',14)
% ylim([-20 110])
% box on; grid on
% fh=gcf;
% fh.PaperUnits='inches';
% fh.PaperPosition=[0 0 8.5 11];
% print('../figures/manuscript/ecco_comparison_G2','-dpng','-r300')

%----- map view

% addpath('../../hikurangi/code/m_map')
% 
% % setup base map
% figure(54); clf; hold on
% lat1=-40; lat2=-38.5;
% lon1=177.5; lon2=179.5;
% m_proj('mercator','longitudes',[lon1 lon2],'latitudes',[lat1 lat2]);
% 
% X=linspace(174.5,180,1320);
% Y=linspace(-37.5,-42,1080);
% Z=imread('../../hikurangi/gshhg-bin-2.3.6/NZ_bathymetry.tiff');
% Z(Z>=0)=NaN;
% [~,hc]=m_contourf(X,Y,Z,-5000:100:0,'linecolor','none');
% 
% m_gshhs_i('patch',[.7 .7 .7]);
% 
% m_contour(X,Y,Z,[-150 -150],'k','linewidth',1) % shelf at 150 m?
% m_text(177.75,-39.55,'150 m','fontsize',12,'rotation',35)
% m_contour(X,Y,Z,[-2750 -2750],'k','linewidth',1) % 2750 m does a decent job visually
% m_text(178.4,-39.9,'2750 m','fontsize',12,'rotation',80)
% 
% colormap(gca,cmocean('-deep'))
% cb1 = colorbar(gca,'eastoutside');
% ylabel(cb1,'Depth (m)')
% set(cb1,'fontsize',14)
% m_grid('xlabeldir','end','fontsize',14);
% 
% % station markers
% m_plot(stalon,stalat,'^k','markersize',10)
% m_text(stalon-0.21,stalat,staname,'fontsize',12)
% 
% % scale arrow
% m_quiver(177.6,-39.95,0,5/20,'off','color','k','linewidth',2)
% m_text(177.65,-39.8,'5 cm/yr','fontsize',12)
% 
% % ECCO2
% h1=m_quiver(stalon(ii),stalat(ii),zeros(10,1),-e_rate/20,'off','color','r','linewidth',2);
% title('ECCO2 inferred deformation')
% 
% fh=gcf;
% fh.PaperUnits='inches';
% fh.PaperPosition=[0 0 11 8.5];
% print('../figures/trend_comps/ecco2/mapview','-dpng','-r300')
% 
% % this one actually compares well, so let's take a difference
% delete(h1)
% h1=m_quiver(stalon(ii),stalat(ii),zeros(10,1),-(a_rate-e_rate)/20,'off','color','r','linewidth',2);
% title('A0A - ECCO2 inferred deformation')
% 
% fh=gcf;
% fh.PaperUnits='inches';
% fh.PaperPosition=[0 0 11 8.5];
% print('../figures/trend_comps/ecco2/mapview_dif','-dpng','-r300')