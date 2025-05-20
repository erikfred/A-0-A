% manuscript_figs_misc.m
%
% figuring some things out before finalizing in other script
%

clear; close all

calcomp=false;
calproc=true;
diffall=false;

if calcomp
    %% CALIBRATION COMPARISONS FOR REDUNDANT GAUGES

    %----- POBS12
    load('../pressure_data/POBS12.mat','calInfoAll1','calInfoAll2','barInfoAll')

    figure(22); clf; hold on
    plot(calInfoAll1.t0p,(calInfoAll1.pCal-calInfoAll1.pCal(1)),'o-','markersize',7,'linewidth',0.5)
    plot(calInfoAll2.t0p,(calInfoAll2.pCal-calInfoAll2.pCal(1)),'s-','markersize',7,'linewidth',0.5)
    plot(barInfoAll.t0p,(barInfoAll.pCal-barInfoAll.pCal(1)),'xk-','markersize',8,'linewidth',0.5)
    ylabel('P (hPa)')
    legend('Gauge 1','Gauge 2','Barometer','location','southwest')
    datetick('x',3)
    set(gca,'fontsize',12)
    box on; grid on

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 4.25 2.75];
    print('../figures/manuscript/cal_demo/G1_G2_barometer','-dpng','-r300')
    print('../figures/manuscript/cal_demo/G1_G2_barometer','-depsc','-vector')

    figure(23); clf; hold on
    plot(calInfoAll1.t0p,(calInfoAll1.pCal-calInfoAll1.pCal(1))-(barInfoAll.pCal-barInfoAll.pCal(1)),'o-','markersize',7,'linewidth',0.5)
    plot(calInfoAll2.t0p,(calInfoAll2.pCal-calInfoAll2.pCal(1))-(barInfoAll.pCal-barInfoAll.pCal(1)),'s-','markersize',7,'linewidth',0.5)
    plot(calInfoAll1.t0p,(calInfoAll1.pCal-calInfoAll1.pCal(1))-(calInfoAll2.pCal-calInfoAll2.pCal(1)),'^k-','markersize',8,'linewidth',0.5)
    ylabel('P (hPa)')
    legend('Gauge 1','Gauge 2','Difference','location','southwest')
    datetick('x',3)
    set(gca,'fontsize',12)
    box on; grid on

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 4.25 2.75];
    print('../figures/manuscript/cal_demo/G1-G2_diff','-dpng','-r300')
    print('../figures/manuscript/cal_demo/G1-G2_diff','-depsc','-vector')
end

if calproc
    %% CALIBRATION PROCESS

    %----- POBS15
    load('../pressure_data/POBS15.mat')

    j=11;
    i0p=calInfoAll1.i0p(j);
    i1p=calInfoAll1.i1p(j);

    figure(5); clf
    subplot(211); hold on
    plot((data.ta(i0p-136:i1p+100) - data.ta(i0p-136))*86400,data.a1(i0p-136:i1p+100),'linewidth',0.5)
    plot((data.ta(i0p-136:i1p+100) - data.ta(i0p-136))*86400,data.a2(i0p-136:i1p+100),'linewidth',0.5)
    ylim([149550 149850])
    xlim([0 505])
    set(gca,'fontsize',12)
    set(gca,'xticklabels',[])
    box on; grid on
    subplot(212); hold on
    plot((data.ta(i0p-136:i1p+100) - data.ta(i0p-136))*86400,data.a2(i0p-136:i1p+100),'linewidth',0.5)
    plot((data.ta(i0p-136:i1p+100) - data.ta(i0p-136))*86400,data.a2(i0p-136:i1p+100),'linewidth',0.5)
    ylim([525 825])
    xlim([0 505])
    set(gca,'fontsize',12)
    ylabel('P (hPa)')
    xlabel('Time (s)')
    box on; grid on

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 4.25 2.75];
    print('../figures/calibrations/demo/cal_a','-dpng','-r300')
    print('../figures/calibrations/demo/cal_a','-depsc','-vector')

    figure(6); clf; hold on
    fill([327 386 386 327],[833.5 833.5 834 834],[0.95 0.95 0.95],'edgecolor','none')
    plot((data.ta(i0p+75:i1p-1) - data.ta(i0p-136))*86400,data.a1(i0p+75:i1p-1),'linewidth',0.5,'color',[0 114 189]/255)
    yline(calInfoAll1.pCal(j),'--','color',[0 114 189]/255)
    ylim([833.5 834])
    set(gca,'fontsize',12)
    ylabel('P1 (hPa)')
    yyaxis right
    plot((data.ta(i0p+75:i1p-1) - data.ta(i0p-136))*86400,data.a2(i0p+75:i1p-1)-0.05,'r','linewidth',0.5)
    yline(calInfoAll2.pCal(j)-0.05,'--r')
    ylim([799 799.5])
    set(gca,'YColor','r')
    ylabel('P2 (hPa)')
    xlabel('Time (s)')
    xlim([205 410])
    box on; grid on

    figure(6); clf; hold on
    fill([327 386 386 327],[833.5 833.5 834 834],[0.95 0.95 0.95],'edgecolor','none')
    h1=plot((data.ta(i0p+75:i1p-1) - data.ta(i0p-136))*86400,data.a1(i0p+75:i1p-1),'linewidth',0.5,'color',[0 114 189]/255);
    yline(calInfoAll1.pCal(j),'--','color',[0 114 189]/255)
    h2=plot((data.ta(i0p+75:i1p-1) - data.ta(i0p-136))*86400,data.a2(i0p+75:i1p-1)+34.45,'r','linewidth',0.5);
    yline(calInfoAll2.pCal(j)+34.45,'--r')
    legend([h1 h2],'P1','P2+34.45 hPa','location','southeast')
    ylim([833.5 834])
    set(gca,'fontsize',12)
    ylabel('P (hPa)')
    xlabel('Time (s)')
    xlim([205 410])
    box on; grid on

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 4.25 2.75];
    print('../figures/calibrations/demo/cal_b','-dpng','-r300')
    print('../figures/calibrations/demo/cal_b','-depsc','-vector')
end

%% COMPARE AGAINST ECCO2, GLORYS
if diffall

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

    figure(82); clf; hold on
    figure(83); clf; hold on
    figure(84); clf; hold on
    for i=1:length(P_list)
        load(['../stitched_data/drift_corrected/' P_list{i}],'dataf');

        % preferred gauge (just 2 for all except POBS4)
        if strcmp(sname(i),'POBS-04')
            ta{i}=dataf.tf; pa{i}=dataf.p1_dcor;
        else
            ta{i}=dataf.tf; pa{i}=dataf.p2_dcor;
        end

        % apply filter defined at start of file
        ta{i}=ta{i}(~isnan(pa{i})); pa{i}=pa{i}(~isnan(pa{i})); pa{i}=pa{i}-mean(pa{i});
        pa{i}=filtfilt(b,a,pa{i});

        % find model index that corresponds to current station
        i_eco=find(strcmp(sname{i},ecco2.nearest_sta));
        i_glo=find(strcmp(sname{i},p_glo.nearest_sta));

        % trim model data to match APG
        te{i}=ecco2.t{i_eco}; pe{i}=ecco2.bpa_int{i_eco};
        pe{i}=pe{i}(te{i}>=ta{i}(1) & te{i}<=ta{i}(end));
        te{i}=te{i}(te{i}>=ta{i}(1) & te{i}<=ta{i}(end));
        pe{i}=interp1(te{i},pe{i},ta{i},'linear','extrap');
        if strcmp(sname(i),'POBS-09') % ECCO wasn't available for full interval
            % apply filter defined at start of file
            pe{i}=filtfilt(b,a,pe{i}(ta{i}<=te{i}(end)));
            pe{i}(ta{i}>te{i}(end))=NaN;
        else
            % apply filter defined at start of file
            pe{i}=filtfilt(b,a,pe{i});
        end

        tg{i}=p_glo.t{i_glo}; pg{i}=p_glo.bpa_int{i_glo};
        pg{i}=pg{i}(tg{i}>=ta{i}(1) & tg{i}<=ta{i}(end));
        tg{i}=tg{i}(tg{i}>=ta{i}(1) & tg{i}<=ta{i}(end));
        pg{i}=interp1(tg{i},pg{i},ta{i},'linear','extrap');
        % apply filter defined at start of file
        pg{i}=filtfilt(b,a,pg{i});
    end

    % plots
    for i=1:length(pa)
        for j=1:length(pa)
            if j==i
                continue
            end

            % differences between APGs
            [tdif,ia,ib]=intersect(round(ta{i},6),round(ta{j},6));
            p1=pa{i}(ia);
            p2=pa{j}(ib);
            pdif=p1-p2;

            figure(82); clf; hold on
            plot(tdif,p1,'linewidth',1)
            plot(tdif,p2,'linewidth',1)
            plot(tdif,pdif,'k','linewidth',2)
            title([sname{i} '-' sname{j} ' (APG only)'])
            ylabel('P (hPa)')
            legend(sname{i},sname{j})
            ylim([-10 10])
            datetick('x',3)
            set(gca,'fontsize',12)
            box on; grid on

            fh=gcf;
            fh.PaperUnits='inches';
            fh.PaperPosition=[0 0 11 8.5];
            print(['../figures/differences/test/' sname{i} '-' sname{j} '_apg'],'-dpng','-r300')

            % differences between ECCO predictions
            [tdif,ia,ib]=intersect(round(ta{i},5),round(ta{j},5));
            p1=pe{i}(ia);
            p2=pe{j}(ib);
            pdif=p1-p2;

            figure(83); clf; hold on
            plot(tdif,p1,'linewidth',1)
            plot(tdif,p2,'linewidth',1)
            plot(tdif,pdif,'k','linewidth',2)
            title([sname{i} '-' sname{j} ' (ECCO only)'])
            ylabel('P (hPa)')
            legend(sname{i},sname{j})
            ylim([-10 10])
            datetick('x',3)
            set(gca,'fontsize',12)
            box on; grid on

            fh=gcf;
            fh.PaperUnits='inches';
            fh.PaperPosition=[0 0 11 8.5];
            print(['../figures/differences/test/' sname{i} '-' sname{j} '_ecco'],'-dpng','-r300')

            % differences between GLORYS predictions
            [tdif,ia,ib]=intersect(round(ta{i},6),round(ta{j},6));
            p1=pg{i}(ia);
            p2=pg{j}(ib);
            pdif=p1-p2;

            figure(84); clf; hold on
            plot(tdif,p1,'linewidth',1)
            plot(tdif,p2,'linewidth',1)
            plot(tdif,pdif,'k','linewidth',2)
            title([sname{i} '-' sname{j} ' (GLORYS only)'])
            ylabel('P (hPa)')
            legend(sname{i},sname{j})
            ylim([-10 10])
            datetick('x',3)
            set(gca,'fontsize',12)
            box on; grid on

            fh=gcf;
            fh.PaperUnits='inches';
            fh.PaperPosition=[0 0 11 8.5];
            print(['../figures/differences/test/' sname{i} '-' sname{j} '_glorys'],'-dpng','-r300')

            % differences between model-corrected APGs
            if i>4 % GLORYS preferred
                pcor1=pa{i}-pg{i};
            else
                pcor1=pa{i}-pe{i};
            end
            if j>4 % GLORYS preferred
                pcor2=pa{j}-pg{j};
            else
                pcor2=pa{j}-pe{j};
            end

            [tdif,ia,ib]=intersect(round(ta{i},6),round(ta{j},6));
            p1=pcor1(ia);
            p2=pcor2(ib);
            pdif=p1-p2;

            figure(85); clf; hold on
            plot(tdif,p1,'linewidth',1)
            plot(tdif,p2,'linewidth',1)
            plot(tdif,pdif,'k','linewidth',2)
            title([sname{i} '-' sname{j} ' (model corrected)'])
            ylabel('P (hPa)')
            legend(sname{i},sname{j})
            ylim([-10 10])
            datetick('x',3)
            set(gca,'fontsize',12)
            box on; grid on

            fh=gcf;
            fh.PaperUnits='inches';
            fh.PaperPosition=[0 0 11 8.5];
            print(['../figures/differences/test/' sname{i} '-' sname{j} '_modcor'],'-dpng','-r300')
        end
    end

end