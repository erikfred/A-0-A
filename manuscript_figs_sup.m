% manuscript_figs_sup.m
%
% one place to build all my supplemental figures
%

clear; close all

% set which to run
s1=false;
s2=false;
s5=false;
s8=true;

%% FIGURE S1,S3,S4 -- gauge justification

if s1
    %----- POBS1
    load('../stitched_data/drift_corrected/POBS1.mat')
    alt=load('../stitched_data/drift_corrected/POBS3.mat','dataf');
    [tdif,ia,ib]=intersect(round(dataf.tf,6),round(alt.dataf.tf,6));
    pdif1=(dataf.p1_dcor(ia))-(alt.dataf.p1_dcor(ib)); pdif1=pdif1-pdif1(1);
    pdif2=(dataf.p2_dcor(ia))-(alt.dataf.p1_dcor(ib)); pdif2=pdif2-pdif2(1);

    m1=polyfit((1:length(tdif))',pdif1,1);
    fit1=polyval(m1,(1:length(tdif))');
    m2=polyfit((1:length(tdif))',pdif2,1);
    fit2=polyval(m2,(1:length(tdif))');

    figure(56); clf;
    subplot(311); hold on
    plot(calInfoAll1.t0p,calInfoAll1.pCal-calInfoAll1.pCal(1),'o-','markersize',12,'linewidth',1)
    plot(calInfoAll2.t0p,calInfoAll2.pCal-calInfoAll2.pCal(1),'s-','markersize',12,'linewidth',1)
    ylabel('P (hPa)')
    yyaxis right
    plot(calInfoAll1.t0p,(calInfoAll1.pCal-calInfoAll1.pCal(1))-(calInfoAll2.pCal-calInfoAll2.pCal(1)),'xk-','markersize',10,'linewidth',2)
    set(gca,'YColor','k')
    legend('Gauge 1','Gauge 2','Difference','location','east')
    datetick('x',3)
    set(gca,'fontsize',12)
    ylabel('\DeltaP (hPa)')
    box on; grid on

    subplot(312); hold on
    plot(dataf.tf,dataf.p1_dcor-dataf.p1_dcor(1),'linewidth',1)
    plot(dataf.tf,dataf.p2_dcor-dataf.p2_dcor(1),'linewidth',1)
    ylabel('P (hPa)')
    yyaxis right
    plot(dataf.tf,(dataf.p1_dcor-dataf.p1_dcor(1))-(dataf.p2_dcor-dataf.p2_dcor(1)),'k','linewidth',2)
    set(gca,'YColor','k')
    legend('Gauge 1','Gauge 2','Difference','location','northeast')
    datetick('x',3)
    set(gca,'fontsize',12)
    ylabel('\DeltaP (hPa)')
    box on; grid on

    subplot(313); hold on
    plot(tdif,pdif1,'linewidth',1)
    plot(tdif,pdif2,'linewidth',1)
    plot(tdif,fit1,'color',[0 114 189]/255, 'linewidth',1)
    plot(tdif,fit2,'r','linewidth',1)
    ylabel('\DeltaP (hPa)')
    legend('Gauge 1 Difference','Gauge 2 Difference','location','northeast')
    datetick('x',3)
    set(gca,'fontsize',12)
    box on; grid on

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print('../figures/manuscript/supplement/figureSX/POBS1_combined','-dpng','-r300')

    %----- POBS4
    load('../stitched_data/drift_corrected/POBS4.mat')
    alt=load('../stitched_data/drift_corrected/POBS2.mat','dataf');
    [tdif,ia,ib]=intersect(round(dataf.tf,6),round(alt.dataf.tf,6));
    pdif1=(dataf.p1_dcor(ia))-(alt.dataf.p1_dcor(ib)); pdif1=pdif1-pdif1(1);
    pdif2=(dataf.p2_dcor(ia))-(alt.dataf.p1_dcor(ib)); pdif2=pdif2-pdif2(1);

    m1=polyfit((1:length(tdif))',pdif1,1);
    fit1=polyval(m1,(1:length(tdif))');
    m2=polyfit((1:length(tdif))',pdif2,1);
    fit2=polyval(m2,(1:length(tdif))');

    figure(56); clf;
    subplot(311); hold on
    plot(calInfoAll1.t0p,calInfoAll1.pCal-calInfoAll1.pCal(1),'o-','markersize',12,'linewidth',1)
    plot(calInfoAll2.t0p,calInfoAll2.pCal-calInfoAll2.pCal(1),'s-','markersize',12,'linewidth',1)
    ylabel('P (hPa)')
    yyaxis right
    plot(calInfoAll1.t0p,(calInfoAll1.pCal-calInfoAll1.pCal(1))-(calInfoAll2.pCal-calInfoAll2.pCal(1)),'xk-','markersize',10,'linewidth',2)
    set(gca,'YColor','k')
    legend('Gauge 1','Gauge 2','Difference','location','east')
    datetick('x',3)
    set(gca,'fontsize',12)
    ylabel('\DeltaP (hPa)')
    box on; grid on

    subplot(312); hold on
    plot(dataf.tf,dataf.p1_dcor-dataf.p1_dcor(1),'linewidth',1)
    plot(dataf.tf,dataf.p2_dcor-dataf.p2_dcor(1),'linewidth',1)
    ylabel('P (hPa)')
    yyaxis right
    plot(dataf.tf,(dataf.p1_dcor-dataf.p1_dcor(1))-(dataf.p2_dcor-dataf.p2_dcor(1)),'k','linewidth',2)
    set(gca,'YColor','k')
    legend('Gauge 1','Gauge 2','Difference','location','northeast')
    datetick('x',3)
    set(gca,'fontsize',12)
    ylabel('\DeltaP (hPa)')
    box on; grid on

    subplot(313); hold on
    plot(tdif,pdif1,'linewidth',1)
    plot(tdif,pdif2,'linewidth',1)
    plot(tdif,fit1,'color',[0 114 189]/255, 'linewidth',1)
    plot(tdif,fit2,'r','linewidth',1)
    ylabel('\DeltaP (hPa)')
    legend('Gauge 1 Difference','Gauge 2 Difference','location','northeast')
    datetick('x',3)
    set(gca,'fontsize',12)
    box on; grid on

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print('../figures/manuscript/supplement/figureSX/POBS4_combined','-dpng','-r300')

    %----- POBS16
    load('../stitched_data/drift_corrected/POBS16.mat')
    alt=load('../stitched_data/drift_corrected/POBS5.mat','dataf');
    [tdif,ia,ib]=intersect(round(dataf.tf,6),round(alt.dataf.tf,6));
    pdif1=(dataf.p1_dcor(ia))-(alt.dataf.p1_dcor(ib)); pdif1=pdif1-pdif1(1);
    pdif2=(dataf.p2_dcor(ia))-(alt.dataf.p1_dcor(ib)); pdif2=pdif2-pdif2(1);

    m1=polyfit((1:length(tdif))',pdif1,1);
    fit1=polyval(m1,(1:length(tdif))');
    m2=polyfit((1:length(tdif))',pdif2,1);
    fit2=polyval(m2,(1:length(tdif))');

    figure(56); clf;
    subplot(311); hold on
    plot(calInfoAll1.t0p,calInfoAll1.pCal-calInfoAll1.pCal(1),'o-','markersize',12,'linewidth',1)
    plot(calInfoAll2.t0p,calInfoAll2.pCal-calInfoAll2.pCal(1),'s-','markersize',12,'linewidth',1)
    ylabel('P (hPa)')
    yyaxis right
    plot(calInfoAll1.t0p,(calInfoAll1.pCal-calInfoAll1.pCal(1))-(calInfoAll2.pCal-calInfoAll2.pCal(1)),'xk-','markersize',10,'linewidth',2)
    set(gca,'YColor','k')
    legend('Gauge 1','Gauge 2','Difference','location','east')
    datetick('x',3)
    set(gca,'fontsize',12)
    ylabel('\DeltaP (hPa)')
    box on; grid on

    subplot(312); hold on
    plot(dataf.tf,dataf.p1_dcor-dataf.p1_dcor(1),'linewidth',1)
    plot(dataf.tf,dataf.p2_dcor-dataf.p2_dcor(1),'linewidth',1)
    ylabel('P (hPa)')
    yyaxis right
    plot(dataf.tf,(dataf.p1_dcor-dataf.p1_dcor(1))-(dataf.p2_dcor-dataf.p2_dcor(1)),'k','linewidth',2)
    set(gca,'YColor','k')
    legend('Gauge 1','Gauge 2','Difference','location','northeast')
    datetick('x',3)
    set(gca,'fontsize',12)
    ylabel('\DeltaP (hPa)')
    box on; grid on

    subplot(313); hold on
    plot(tdif,pdif1,'linewidth',1)
    plot(tdif,pdif2,'linewidth',1)
    plot(tdif,fit1,'color',[0 114 189]/255, 'linewidth',1)
    plot(tdif,fit2,'r','linewidth',1)
    ylabel('\DeltaP (hPa)')
    legend('Gauge 1 Difference','Gauge 2 Difference','location','northeast')
    datetick('x',3)
    set(gca,'fontsize',12)
    box on; grid on

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print('../figures/manuscript/supplement/figureSX/POBS16_combined','-dpng','-r300')
end

%% FIGURE S2 -- show all pressure records and calibrations

if s2
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

    figure(4); clf;
    for i=1:length(P_list)
        load(['../stitched_data/drift_corrected/' P_list{i}]);
        if strcmp(sname{i},'POBS-07')
            load('../pressure_data/POBS7_co-located-w-QA15.mat','barInfoAll');
        else
            load(['../pressure_data/' P_list{i}],'barInfoAll');
        end

        % pressure time series
        tplot=dataf.tf;
        plot1a=dataf.p1f; plot1a=plot1a-mean(plot1a,'omitmissing');
        plot1b=dataf.p1_dcor; plot1b=plot1b-mean(plot1b,'omitmissing');
        plot2a=dataf.p2f; plot2a=plot2a-mean(plot2a,'omitmissing');
        plot2b=dataf.p2_dcor; plot2b=plot2b-mean(plot2b,'omitmissing');

        subplot(121); hold on
        plot(tplot,plot1a+10*(i-1),'r','linewidth',1)
        plot(tplot,plot1b+10*(i-1),'color',[0 114 189]/255, 'linewidth',1)
        if strcmp(sname{i},'POBS-09')
            xpo=tplot(1)-60; ypo=mean([plot1a(end-99:end);plot1b(end-99:end)],'omitmissing')+(i-1)*10;
        elseif any(isnan(plot1a))
            inan=find(isnan(plot1a),1)-1;
            xpo=tplot(inan)+10; ypo=mean([plot1a(inan-99:inan);plot1b(inan-99:inan)],'omitmissing')+(i-1)*10;
        else
            xpo=tplot(end)+10; ypo=mean([plot1a(end-99:end);plot1b(end-99:end)],'omitmissing')+(i-1)*10;
        end
        text(xpo,ypo,sname{i},'fontsize',12)

        subplot(122); hold on
        plot(tplot,plot2a+10*(i-1),'r','linewidth',1)
        plot(tplot,plot2b+10*(i-1),'color',[0 114 189]/255, 'linewidth',1)
        if strcmp(sname{i},'POBS-09')
            xpo=tplot(1)-60; ypo=mean([plot2a(end-99:end);plot2b(end-99:end)],'omitmissing')+(i-1)*10;
        elseif any(isnan(plot2a))
            inan=find(isnan(plot2a),1)-1;
            xpo=tplot(inan)+10; ypo=mean([plot2a(inan-99:inan);plot2b(inan-99:inan)],'omitmissing')+(i-1)*10;
        else
            xpo=tplot(end)+10; ypo=mean([plot2a(end-99:end);plot2b(end-99:end)],'omitmissing')+(i-1)*10;
        end
        text(xpo,ypo,sname{i},'fontsize',12)

        % calibrations
        ical=calInfoAll1.t0p>=tplot(1) & calInfoAll1.t0p<=tplot(end);
        tcal=calInfoAll1.t0p(ical);
        
        if strcmp(sname{i},'POBS-02')
            cal1=calInfoAll1.pCal(ical)-barInfoAll.T(ical)*3.71; cal1=cal1-mean(cal1,'omitmissing');
            cal2=calInfoAll2.pCal(ical)-barInfoAll.T(ical)*3.71; cal2=cal2-mean(cal2,'omitmissing');
        elseif any(isnan(plot1a))
            inan=find(isnan(plot1a),1); tnan=tplot(inan);
            cal1=calInfoAll1.pCal(ical)-barInfoAll.pCal(ical); cal1=cal1-mean(cal1,'omitmissing'); cal1(tcal>=tnan)=NaN;
            cal2=calInfoAll2.pCal(ical)-barInfoAll.pCal(ical); cal2=cal2-mean(cal2,'omitmissing');
        elseif any(isnan(plot2a))
            inan=find(isnan(plot2a),1); tnan=tplot(inan);
            cal1=calInfoAll1.pCal(ical)-barInfoAll.pCal(ical); cal1=cal1-mean(cal1,'omitmissing');
            cal2=calInfoAll2.pCal(ical)-barInfoAll.pCal(ical); cal2=cal2-mean(cal2,'omitmissing'); cal2(tcal>=tnan)=NaN;
        else
            cal1=calInfoAll1.pCal(ical)-barInfoAll.pCal(ical); cal1=cal1-mean(cal1,'omitmissing');
            cal2=calInfoAll2.pCal(ical)-barInfoAll.pCal(ical); cal2=cal2-mean(cal2,'omitmissing');
        end

        subplot(121); hold on
        plot(tcal,cal1+10*(i-1),'ko:','linewidth',1)

        subplot(122); hold on
        plot(tcal,cal2+10*(i-1),'ko:','linewidth',1)
    end

    subplot(121)
    xlim([datenum(2022,10,01) datenum(2024,04,15)])
    datetick('x',3,'keeplimits')
    ylabel('P (hPa)')
    legend('Pre-correction','Post-correction','Calibrations','location','northeast')
    set(gca,'fontsize',12)
    ylim([-10 130])
    box on; grid on

    subplot(122)
    xlim([datenum(2022,10,01) datenum(2024,04,15)])
    datetick('x',3,'keeplimits')
    ylabel('P (hPa)')
    legend('Pre-correction','Post-correction','Calibrations','location','northeast')
    set(gca,'fontsize',12)
    ylim([-10 130])
    box on; grid on

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print('../figures/manuscript/supplement/figureS2','-dpng','-r300')
end

%% FIGURE S5 -- compare barometer to ideal gas law

if s5
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

    figure(5); clf;
    n=0;
    for i=1:length(P_list)
        if strcmp(sname{i},'POBS-07')
            load('../pressure_data/POBS7_co-located-w-QA15.mat','calInfoAll1','calInfoAll2','barInfoAll');
        else
            load(['../pressure_data/' P_list{i}],'calInfoAll1','calInfoAll2','barInfoAll');
        end

        if strcmp(sname{i},'POBS-15')
            tplot=barInfoAll.t0p(1:end-10);
            bplot=barInfoAll.pCal(1:end-10); bplot=bplot-mean(bplot,'omitmissing');
            Tplot=barInfoAll.T(1:end-10)*3.71; Tplot=Tplot-mean(Tplot,'omitmissing');
        else
            tplot=barInfoAll.t0p;
            bplot=barInfoAll.pCal; bplot=bplot-mean(bplot,'omitmissing');
            Tplot=barInfoAll.T*3.71; Tplot=Tplot-mean(Tplot,'omitmissing');
        end

        if any(i==[1 3 7 10 11 12])
            ip=find(i==[1 3 7 10 11 12]);
            imap=[2 4 6 8 10 12]; imap=imap(ip);
            subplot(6,2,imap); hold on
            plot(tplot,bplot,'ko-','linewidth',1)
            plot(tplot,Tplot,'m^-','linewidth',1)

            if strcmp(sname{i},'POBS-09')
                xpo=tplot(1)-80; ypo=mean(Tplot(end-4:end)+3,'omitmissing');
            else
                xpo=tplot(end)+20; ypo=mean(Tplot(end-4:end),'omitmissing');
            end
        else
            n=n+1;
            subplot(121); hold on
            plot(tplot,bplot+5*(n-1),'ko-','linewidth',1)
            plot(tplot,Tplot+5*(n-1),'m^-','linewidth',1)

            xpo=tplot(end)+20; ypo=mean(Tplot(end-4:end),'omitmissing')+5*(n-1);
        end        
        text(xpo,ypo,sname{i},'fontsize',12)
    end

    subplot(121)
    xlim([datenum(2022,10,01) datenum(2024,04,15)])
    datetick('x',3,'keeplimits')
    ylabel('P (hPa)')
    legend('Barometer','Ideal Gas Law','location','northeast')
    set(gca,'fontsize',12)
    box on; grid on

    subplot(622)
    xlim([datenum(2022,10,01) datenum(2024,04,15)])
    datetick('x',3,'keeplimits')
    ylabel('P (hPa)')
    % legend('Barometer','Ideal Gas Law','location','northeast')
    set(gca,'fontsize',12)
    box on; grid on

    subplot(624)
    xlim([datenum(2022,10,01) datenum(2024,04,15)])
    datetick('x',3,'keeplimits')
    ylabel('P (hPa)')
    % legend('Barometer','Ideal Gas Law','location','northeast')
    set(gca,'fontsize',12)
    box on; grid on

    subplot(626)
    xlim([datenum(2022,10,01) datenum(2024,04,15)])
    datetick('x',3,'keeplimits')
    ylabel('P (hPa)')
    % legend('Barometer','Ideal Gas Law','location','northeast')
    set(gca,'fontsize',12)
    box on; grid on

    subplot(628)
    xlim([datenum(2022,10,01) datenum(2024,04,15)])
    datetick('x',3,'keeplimits')
    ylabel('P (hPa)')
    % legend('Barometer','Ideal Gas Law','location','northeast')
    set(gca,'fontsize',12)
    box on; grid on

    subplot(6,2,10)
    xlim([datenum(2022,10,01) datenum(2024,04,15)])
    datetick('x',3,'keeplimits')
    ylabel('P (hPa)')
    % legend('Barometer','Ideal Gas Law','location','northeast')
    set(gca,'fontsize',12)
    box on; grid on

    subplot(6,2,12)
    xlim([datenum(2022,10,01) datenum(2024,04,15)])
    datetick('x',3,'keeplimits')
    ylabel('P (hPa)')
    % legend('Barometer','Ideal Gas Law','location','northeast')
    set(gca,'fontsize',12)
    box on; grid on

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print('../figures/manuscript/supplement/figureS5','-dpng','-r300')
end

%% FIGURE S7 -- analog to Figure 4 that shows effect of further trimming

if s8
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

    load('../../hikurangi/ocean_data/ECCO2/ecco2_pressure_update2.mat')
    load('../../hikurangi/ocean_data/GLORYS/glorys_pressure_update2.mat')

    % low-pass filter parameters
    [b,a]=butter(4,(2/(30*24)),'low'); % assumes hourly sample rate

    a_rate=nan(12,1);
    e_rate=nan(12,1);
    g_rate=nan(12,1);

    figure(4); clf; hold on
    for i=1:length(P_list)
        load(['../stitched_data/drift_corrected/' P_list{i}],'dataf');

        % preferred gauge (just 2 for all except POBS4)
        if strcmp(sname(i),'POBS-04')
            ta=dataf.tf; pa=dataf.p1_dcor;
        else
            ta=dataf.tf; pa=dataf.p2_dcor;
        end

        % remove all data prior to December 1 2022
        ta=ta(~isnan(pa)); pa=pa(~isnan(pa)); pa=pa-mean(pa);
        itc=ta<datenum(2022,12,01);
        ta(itc)=[]; pa(itc)=[];

        % apply filter defined at start of file
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

        figure(4); % plot mixed comparisons
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

    figure(4)
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
    print('../figures/manuscript/supplement/truncated_comparison','-dpng','-r300')
    print('../figures/manuscript/supplement/truncated_comparison','-depsc','-vector')
end