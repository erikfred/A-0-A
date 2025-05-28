% manuscript_figs_1.m
%
% Generates stacked plots of pressure data before and after drift
% correction (both gauges), as well as a boxplot (scatter?) that
% demonstrates the change in trend as a result.
%

clear; close all

f3=true;
f4=false;

load('../pressure_data/geometry');

POBS_dir=dir('../stitched_data/drift_corrected');
f_list={POBS_dir.name}';
file_check=cellfun(@(v)v(1),f_list);
i_list=find(eq(file_check,'P')); % subselects only files begining with 'P'
POBS_list=f_list(i_list);

% first, sort geometry info by depth
[sd,id]=sort(stadepth,'descend');
sla=stalat(id);
slo=stalon(id);
sname=staname(id);
fname=filename(id);

% map from one sorting to the other
for i=1:length(fname)
    if length(fname{i})>6
        fname{i}='POBS7';
    end
    ii(i)=find(strcmp([fname{i} '.mat'],POBS_list));
end
P_list=POBS_list(ii); % should be depth-sorted

if f3
    figure(4); clf; hold on
    figure(5); clf; hold on
    figure(6); clf; hold on
    for i=1:length(P_list)
        load(['../stitched_data/drift_corrected/' P_list{i}],'dataf');
        if strcmp(P_list{i},'POBS7.mat')
            load('../pressure_data/POBS7_co-located-w-QA15.mat','calInfoAll1')
        else
            load(['../pressure_data/' P_list{i}],'calInfoAll1')
        end

        % exclude interval spaning first 2 calibrations
        if strcmp(P_list{i},'POBS2.mat') % POBS2 need an extra week
            ikeep=dataf.tf>calInfoAll1.t0p(3);
        else
            ikeep=dataf.tf>calInfoAll1.t0p(2);
        end

        ta1=dataf.tf(ikeep); pa1=dataf.p1f(ikeep);
        ta1=ta1(~isnan(pa1)); pa1=pa1(~isnan(pa1)); pa1=pa1-mean(pa1);
        ma1=polyfit(ta1,pa1,1);
        fa1=polyval(ma1,ta1);

        ta2=dataf.tf(ikeep); pa2=dataf.p2f(ikeep);
        ta2=ta2(~isnan(pa2)); pa2=pa2(~isnan(pa2)); pa2=pa2-mean(pa2);
        ma2=polyfit(ta2,pa2,1);
        fa2=polyval(ma2,ta2);

        tb1=dataf.tf(ikeep); pb1=dataf.p1_dcor(ikeep);
        tb1=tb1(~isnan(pb1)); pb1=pb1(~isnan(pb1)); pb1=pb1-mean(pb1);
        mb1=polyfit(tb1,pb1,1);
        fb1=polyval(mb1,tb1);

        tb2=dataf.tf(ikeep); pb2=dataf.p2_dcor(ikeep);
        tb2=tb2(~isnan(pb2)); pb2=pb2(~isnan(pb2)); pb2=pb2-mean(pb2);
        mb2=polyfit(tb2,pb2,1);
        fb2=polyval(mb2,tb2);

        figure(4) % Gauge 1
        plot(ta1,pa1+(i-1)*10,'r','linewidth',1)
        plot(tb1,pb1+(i-1)*10,'color',[0 114 189]/255,'linewidth',1)
        plot(ta1,fa1+(i-1)*10,'r','linewidth',1)
        plot(tb1,fb1+(i-1)*10,'color',[0 114 189]/255,'linewidth',1)
        if strcmp(sname(i),'POBS-09')
            text(ta1(1)-60,mean([pa1(end-99:end);pb1(end-99:end)])+(i-1)*10,sname(i),'fontsize',12)
        else
            text(ta1(end)+10,mean([pa1(end-99:end);pb1(end-99:end)])+(i-1)*10,sname(i),'fontsize',12)
        end

        figure(5) % Gauge 2
        plot(ta2,pa2+(i-1)*10,'r','linewidth',1)
        plot(tb2,pb2+(i-1)*10,'color',[0 114 189]/255,'linewidth',1)
        plot(ta2,fa2+(i-1)*10,'r','linewidth',1)
        plot(tb2,fb2+(i-1)*10,'color',[0 114 189]/255,'linewidth',1)
        if strcmp(sname(i),'POBS-09')
            text(ta2(1)-60,mean([pa2(end-99:end);pb2(end-99:end)])+(i-1)*10,sname(i),'k','fontsize',12)
        else
            text(ta2(end)+10,mean([pa2(end-99:end);pb2(end-99:end)])+(i-1)*10,sname(i),'k','fontsize',12)
        end

        figure(6) % Best Gauge (just 2 for all except POBS4)
        if strcmp(sname(i),'POBS-04')
            plot(ta1,pa1+(i-1)*10,'r','linewidth',0.5)
            plot(tb1,pb1+(i-1)*10,'color',[0 114 189]/255,'linewidth',0.5)
            plot(ta1,fa1+(i-1)*10,'r','linewidth',1)
            plot(tb1,fb1+(i-1)*10,'color',[0 114 189]/255,'linewidth',1)
        else
            plot(ta2,pa2+(i-1)*10,'r','linewidth',0.5)
            plot(tb2,pb2+(i-1)*10,'color',[0 114 189]/255,'linewidth',0.5)
            plot(ta2,fa2+(i-1)*10,'r','linewidth',1)
            plot(tb2,fb2+(i-1)*10,'color',[0 114 189]/255,'linewidth',1)
        end
        if strcmp(sname(i),'POBS-09')
            text(ta2(1)-60,mean([pa2(end-99:end);pb2(end-99:end)])+(i-1)*10,sname(i),'k','fontsize',10)
        else
            text(ta2(end)+10,mean([pa2(end-99:end);pb2(end-99:end)])+(i-1)*10,sname(i),'k','fontsize',10)
        end

        a1r(i)=ma1(1)*365;
        b1r(i)=mb1(1)*365;
        a2r(i)=ma2(1)*365;
        b2r(i)=mb2(1)*365;
    end

    figure(4)
    xlim([datenum(2022,10,01) datenum(2024,04,15)])
    datetick('x',3,'keeplimits')
    ylabel('P (hPa)')
    legend('Pre-correction','Post-correction','location','northeast')
    title('Gauge 1')
    set(gca,'fontsize',14)
    ylim([-10 130])
    box on; grid on
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print('../figures/manuscript/G1_corrected_pstack','-dpng','-r300')
    print('../figures/manuscript/G1_corrected_pstack','-depsc','-vector')

    figure(5)
    xlim([datenum(2022,10,01) datenum(2024,04,15)])
    datetick('x',3,'keeplimits')
    ylabel('P (hPa)')
    legend('Pre-correction','Post-correction','location','northeast')
    title('Gauge 2')
    set(gca,'fontsize',14)
    ylim([-10 130])
    box on; grid on
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print('../figures/manuscript/G2_corrected_pstack','-dpng','-r300')
    print('../figures/manuscript/G2_corrected_pstack','-depsc','-vector')

    figure(6)
    xlim([datenum(2022,10,01) datenum(2024,04,15)])
    datetick('x',3,'keeplimits')
    ylabel('P (hPa)')
    legend('Pre-correction','Post-correction','location','northeast')
    set(gca,'fontsize',12)
    ylim([-10 130])
    box on; grid on
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 5.5 8.5];
    print('../figures/manuscript/corrected_pstack','-dpng','-r300')
    print('../figures/manuscript/corrected_pstack','-depsc','-vector')

    figure(7); clf; hold on
    stem((1:12)-0.1,a1r,'ro','markersize',7,'markerfacecolor','r')
    stem((1:12)-0.1,b1r','o','markersize',7,'color',[0 114 189]/255,'markerfacecolor',[0 114 189]/255)
    stem((1:12)+0.1,a2r,'rs','markersize',7,'markerfacecolor','r')
    stem((1:12)+0.1,b2r','s','markersize',7,'color',[0 114 189]/255,'markerfacecolor',[0 114 189]/255)
    view(90,-90)
    ylabel('Slope (hPa/yr)')
    legend('G1 pre','G1 post','G2 pre','G2 post','location','northwest')
    set(gca,'fontsize',12)
    set(gca,'xtick',1:12)
    set(gca,'xticklabels',sname)
    xlim([0 13])
    box on; grid on
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 5.5 8.5];
    print('../figures/manuscript/trend_summary','-dpng','-r300')
    print('../figures/manuscript/trend_summary','-depsc','-vector')
end

%% COMPARE AGAINST ECCO2, GLORYS
if f4
    load('../../hikurangi/ocean_data/ECCO2/ecco2_pressure_update2.mat')
    load('../../hikurangi/ocean_data/GLORYS/glorys_pressure_update2.mat')

    % low-pass filter parameters
    [b,a]=butter(4,(2/(30*24)),'low'); % assumes hourly sample rate

    a_rate=nan(12,1);
    e_rate=nan(12,1);
    g_rate=nan(12,1);

    figure(82); clf; hold on
    figure(83); clf; hold on
    figure(84); clf; hold on
    for i=1:length(P_list)
        load(['../stitched_data/drift_corrected/' P_list{i}],'dataf');

        % preferred gauge (just 2 for all except POBS4)
        if strcmp(sname(i),'POBS-04')
            ta=dataf.tf; pa=dataf.p1_dcor;
        else
            ta=dataf.tf; pa=dataf.p2_dcor;
        end

        % apply filter defined at start of file
        ta=ta(~isnan(pa)); pa=pa(~isnan(pa)); pa=pa-mean(pa);
        pa=filtfilt(b,a,pa);

        % linear fit to observation
        ta=ta(~isnan(pa)); pa=pa(~isnan(pa)); pa=pa-mean(pa);
        ma=polyfit(ta-ta(1),pa,1);
        fa=polyval(ma,ta-ta(1));
        a_rate(i)=ma(1)*365; % store A-0-A rate

        % find model index that corresponds to current station
        i_eco=find(strcmp(sname{i},ecco2.nearest_sta));
        i_glo=find(strcmp(sname{i},p_glo.nearest_sta));

        % trim model data to match APG
        te=ecco2.t{i_eco}; pe=ecco2.bpa_int{i_eco};
        pe=pe(te>=ta(1) & te<=ta(end));
        te=te(te>=ta(1) & te<=ta(end));
        pe=interp1(te,pe,ta,'linear','extrap');
        if strcmp(sname(i),'POBS-09') % ECCO wasn't available for full interval
            % apply filter defined at start of file
            pe=filtfilt(b,a,pe(ta<=te(end)));
            pe(ta>te(end))=NaN;
        else
            % apply filter defined at start of file
            pe=filtfilt(b,a,pe);
        end

        tg=p_glo.t{i_glo}; pg=p_glo.bpa_int{i_glo};
        pg=pg(tg>=ta(1) & tg<=ta(end));
        tg=tg(tg>=ta(1) & tg<=ta(end));
        pg=interp1(tg,pg,ta,'linear','extrap');
        % apply filter defined at start of file
        pg=filtfilt(b,a,pg);

        % linear fits to models
        if strcmp(sname(i),'POBS-09') % ECCO wasn't available for full interval
            ikeep=~isnan(pe);
            me=polyfit(ta(ikeep)-ta(1),pe(ikeep),1);
            fe=polyval(me,ta-ta(1));
            fe(~ikeep)=NaN;

            %% THIS PART MAKES INCIDENTALLY MESSES UP APG COMP AGAINST GLORYS
            % re-do A-0-A calc for apples-to-apples comparison
            ma=polyfit(ta(ikeep)-ta(1),pa(ikeep),1);
            fa=polyval(ma,ta(ikeep)-ta(1));
            a_rate(i)=ma(1)*365; % store A-0-A rate
            fa(~ikeep)=NaN;
        else
            me=polyfit(ta-ta(1),pe,1);
            fe=polyval(me,ta-ta(1));
        end
        e_rate(i)=me(1)*365; % store ECCO2 rate
        mg=polyfit(ta-ta(1),pg,1);
        g_rate(i)=mg(1)*365; % store GLORYS rate
        fg=polyval(mg,ta-ta(1));

        figure(82); % plot ECCO comparisons
        plot(ta,pa+10*(i-1),'color',[0 114 189]/255,'linewidth',0.5)
        plot(ta,pe+10*(i-1),'r','linewidth',0.5)
        plot(ta,pa-pe+10*(i-1),'k','linewidth',1)
        plot(ta,fa+10*(i-1),'color',[0 114 189]/255,'linewidth',1)
        plot(ta,fe+10*(i-1),'r','linewidth',1)
        plot(ta,fa-fe+10*(i-1),'k','linewidth',2)
        if strcmp(sname(i),'POBS-09')
            text(ta(1)-60,mean(pa(1:100)-pe(1:100))+(i-1)*10,{sname{i};[num2str(round(a_rate(i),1)) ' hPa/yr'];...
                [num2str(round(a_rate(i)-e_rate(i),1)) ' hPa/yr']},'k','fontsize',12)
        else
            text(ta(end)+10,mean(pa(end-99:end)-pe(end-99:end))+(i-1)*10,{sname{i};[num2str(round(a_rate(i),1)) ' hPa/yr'];...
                [num2str(round(a_rate(i)-e_rate(i),1)) ' hPa/yr']},'k','fontsize',12)
        end

        figure(83); % plot GLORYS comparisons
        plot(ta,pa+10*(i-1),'color',[0 114 189]/255,'linewidth',0.5)
        plot(ta,pg+10*(i-1),'r','linewidth',0.5)
        plot(ta,pa-pg+10*(i-1),'k','linewidth',1)
        plot(ta,fa+10*(i-1),'color',[0 114 189]/255,'linewidth',1)
        plot(ta,fg+10*(i-1),'r','linewidth',1)
        plot(ta,fa-fg+10*(i-1),'k','linewidth',2)
        if strcmp(sname(i),'POBS-09')
            text(ta(1)-60,mean(pa(1:100)-pg(1:100))+(i-1)*10,{sname{i};[num2str(round(a_rate(i),1)) ' hPa/yr'];...
                [num2str(round(a_rate(i)-g_rate(i),1)) ' hPa/yr']},'k','fontsize',12)
        else
            text(ta(end)+10,mean(pa(end-99:end)-pg(end-99:end))+(i-1)*10,{sname{i};[num2str(round(a_rate(i),1)) ' hPa/yr'];...
                [num2str(round(a_rate(i)-g_rate(i),1)) ' hPa/yr']},'k','fontsize',12)
        end

        figure(84); % plot mixed comparisons
        if i>4
            plot(ta,pa+10*(i-1),'color',[0 114 189]/255,'linewidth',0.5)
            plot(ta,pg+10*(i-1),'r','linewidth',0.5)
            plot(ta,pa-pg+10*(i-1),'k','linewidth',1)
            plot(ta,fa+10*(i-1),'color',[0 114 189]/255,'linewidth',1)
            plot(ta,fg+10*(i-1),'r','linewidth',1)
            plot(ta,fa-fg+10*(i-1),'k','linewidth',2)
            if strcmp(sname(i),'POBS-09')
                text(ta(1)-60,mean(pa(1:100)-pg(1:100))+(i-1)*10,{sname{i};[num2str(round(a_rate(i),1)) ' hPa/yr'];...
                [num2str(round(a_rate(i)-g_rate(i),1)) ' hPa/yr']},'k','fontsize',12)
            else
                text(ta(end)+10,mean(pa(end-99:end)-pg(end-99:end))+(i-1)*10,{sname{i};[num2str(round(a_rate(i),1)) ' hPa/yr'];...
                [num2str(round(a_rate(i)-g_rate(i),1)) ' hPa/yr']},'k','fontsize',12)
            end
        else
            plot(ta,pa+10*(i-1),'color',[0 114 189]/255,'linewidth',0.5)
            plot(ta,pe+10*(i-1),'r','linewidth',0.5)
            plot(ta,pa-pe+10*(i-1),'k','linewidth',1)
            plot(ta,fa+10*(i-1),'color',[0 114 189]/255,'linewidth',1)
            plot(ta,fe+10*(i-1),'r','linewidth',1)
            plot(ta,fa-fe+10*(i-1),'k','linewidth',2)
            if strcmp(sname(i),'POBS-09')
                text(ta(1)-60,mean(pa(1:100)-pe(1:100))+(i-1)*10,{sname{i};[num2str(round(a_rate(i),1)) ' hPa/yr'];...
                [num2str(round(a_rate(i)-e_rate(i),1)) ' hPa/yr']},'k','fontsize',12)
            else
                text(ta(end)+10,mean(pa(end-99:end)-pe(end-99:end))+(i-1)*10,{sname{i};[num2str(round(a_rate(i),1)) ' hPa/yr'];...
                [num2str(round(a_rate(i)-e_rate(i),1)) ' hPa/yr']},'k','fontsize',12)
            end
        end
        
    end

    figure(82)
    ylabel('P (hPa)')
    legend('APG','ECCO2','Difference','location','northeast')
    xlim([datenum(2022,10,01) datenum(2024,04,15)])
    datetick('x',3,'keeplimits')
    set(gca,'fontsize',12)
    ylim([-20 120])
    box on; grid on
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 5.5 8.5];
    % print('../figures/manuscript/supplement/figureSY/ecco_comparison','-dpng','-r300')
    % print('../figures/manuscript/supplement/figureSY/ecco_comparison','-depsc','-vector')

    figure(83)
    ylabel('P (hPa)')
    legend('APG','GLORYS','Difference','location','northeast')
    xlim([datenum(2022,10,01) datenum(2024,04,15)])
    datetick('x',3,'keeplimits')
    set(gca,'fontsize',12)
    ylim([-20 120])
    box on; grid on
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 5.5 8.5];
    % print('../figures/manuscript/supplement/figureSY/glorys_comparison','-dpng','-r300')
    % print('../figures/manuscript/supplement/figureSY/glorys_comparison','-depsc','-vector')

    figure(84)
    ylabel('P (hPa)')
    legend('APG','Model','Difference','location','northeast')
    xlim([datenum(2022,10,01) datenum(2024,04,15)])
    datetick('x',3,'keeplimits')
    set(gca,'fontsize',12)
    ylim([-20 120])
    box on; grid on
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    % print('../figures/manuscript/mixed_comparison','-dpng','-r300')
    % print('../figures/manuscript/mixed_comparison','-depsc','-vector')

    %----- map view

    addpath('../../hikurangi/code/m_map')

    % setup base map
    figure(54); clf; hold on
    lat1=-40; lat2=-38.5;
    lon1=177.5; lon2=179.5;
    m_proj('mercator','longitudes',[lon1 lon2],'latitudes',[lat1 lat2]);

    X=linspace(174.5,180,1320);
    Y=linspace(-37.5,-42,1080);
    Z=imread('../../hikurangi/gshhg-bin-2.3.6/NZ_bathymetry.tiff');
    Z(Z>=0)=NaN;
    % [~,hc]=m_contourf(X,Y,Z,-5000:100:0,'linecolor','none');

    m_gshhs_i('patch',[.7 .7 .7]);

    m_contour(X,Y,Z,[-150 -150],'k','linewidth',1) % shelf at 150 m?
    m_text(177.75,-39.55,'150 m','fontsize',12,'rotation',35)
    m_contour(X,Y,Z,[-2750 -2750],'k','linewidth',1) % 2750 m does a decent job visually
    m_text(178.4,-39.9,'2750 m','fontsize',12,'rotation',80)

    % colormap(gca,cmocean('-deep'))
    % cb1 = colorbar(gca,'eastoutside');
    % ylabel(cb1,'Depth (m)')
    % set(cb1,'fontsize',12)
    m_grid('xlabeldir','end','fontsize',12);

    % station markers
    m_plot(slo,sla,'^k','markersize',10)
    m_text(slo-0.21,sla,sname,'fontsize',12)

    % scale arrow
    m_quiver(177.6,-39.95,0,1/20,'off','color','k','linewidth',2)
    m_text(177.65,-39.925,'1 cm/yr','fontsize',12)

    % A-0-A observations
    h1=m_quiver(slo-0.01,sla,zeros(12,1),-a_rate/20,'off','color','c','linewidth',2);

    % model-corrected observations
    c_rate=a_rate-g_rate; c_rate(1:4)=a_rate(1:4)-e_rate(1:4);
    h2=m_quiver(slo+0.01,sla,zeros(12,1),-c_rate/20,'off','color','k','linewidth',2);

    legend([h1 h2],'Drift Corrected','Drift+Model Corrected','location','southeast','fontsize',12)

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 5.5];
    print('../figures/manuscript/mapview_def','-dpng','-r300')
    print('../figures/manuscript/mapview_def','-depsc','-vector')
end

%% Analog to Figure 3b that shows model corrections

if fz
    figure(8); clf; hold on
    glo=[5.9; 6.7; 4.4; 4.8; -1.4; 0; 1.9; -0.4; 3.7; 0.3; 1.8; -1.7]';
    hyc=[4.6; 3.5; 2.7; 4.2; -5.3; 1.4; 1.7; 1.2; 4.9; 1.9; 2.8; -0.5];
    stem((1:12),b2r','s','markersize',10,'color',[0 114 189]/255,'markerfacecolor',[0 114 189]/255)
    stem((1:12)-0.1,glo,'kp','markersize',15,'markerfacecolor','y')
    stem((1:12)+0.1,hyc,'kd','markersize',10,'markerfacecolor','y')
    % view(90,90)
    ylabel('Slope (hPa/yr)')
    legend('drift corrected','GLORYS-corrected','ECCO-corrected','location','southeast')
    set(gca,'fontsize',14)
    set(gca,'xtick',1:12)
    namelist={POBS_dir(i_list).name};
    for k=1:length(namelist)
        namelist{k}=namelist{k}(1:end-4);
    end
    set(gca,'xticklabels',namelist)
    xlim([0 13])
    ylim([-6 10])
    box on; grid on
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print('../figures/manuscript/modtrend_summary','-dpng','-r300')
    print('../figures/manuscript/modtrend_summary','-depsc','-vector')
end