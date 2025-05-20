% inspect_stitch.m
%
% Compares stitched A-0-A against tidal model and against depth-matched
% reference to help evaluate smoothness of data segments in vicinity of
% calibration interval.
%
% NOTE: DIFFERENCING SECTIONS HAVE BEEN COMMENTED OUT UNTIL KU DATA COME IN
% AND WE CAN FIND MEANINGFUL MATCHES FOR ALL DEPTHS
%

clear; close all

G2=load('../../hikurangi/processed_data/GONDOR_II_1Hz.mat');
load('../pressure_data_Y2/geometry');

path(path,'../../EKF_spotl/matlab')

POBS_dir=dir('../pressure_data_Y2');
flist={POBS_dir.name}';
file_check=cellfun(@(v)v(1),flist);
i_list=find(eq(file_check,'P')); % subselects only files begining with 'P'
POBS_list=flist(i_list);

for i=1:length(POBS_list)

    % if ~strcmp(POBS_dir(i).name,'POBS09.mat') && ~strcmp(POBS_dir(i).name,'POBS10.mat')
    %     continue
    % end

    load([POBS_dir(1).folder '/' POBS_list{i}])
    t=data.ta(data.lNormState1);
    p1=data.a1(data.lNormState1);
    p2=data.a2(data.lNormState1);
    T1=data.Ta1(data.lNormState1);
    T2=data.Ta2(data.lNormState1);

    % identifier for station
    shortname=staname{i};
    
    %% redundant gauge differences

    % expand calibrations onto 1 Hz sampling
    cal=load(['../stitched_data_Y2/drift_corrected/' POBS_list{i}]);
    c1=interp1(cal.dataf.tf,cal.dataf.p1f-cal.dataf.p1_dcor,t);
    c2=interp1(cal.dataf.tf,cal.dataf.p2f-cal.dataf.p2_dcor,t);

    % plot pressure over full time interval
    figure(200); clf; hold on
    plot(t,(p1-mean(p1))-(p2-mean(p2)),'linewidth',1)
    plot(t,(c1-nanmean(c1))-(c2-nanmean(c2)),'linewidth',1)
    plot(t,(p1-mean(p1))-(p2-mean(p2))-((c1-nanmean(c1))-(c2-nanmean(c2))),'k','linewidth',2)
    ylabel('P (hPa)')
    legend('P_e_x_t','Calibrations','P_c_o_r','location','northwest')
    datetick('x')
    set(gca,'fontsize',14)
    box on; grid on

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    keyboard % adjust ylim, etc.
    if ~exist(['../figures_Y2/stitching/stitch_assessment/redundant_gauge/' POBS_list{i}(1:end-4)],'dir')
        mkdir(['../figures_Y2/stitching/stitch_assessment/redundant_gauge/' POBS_list{i}(1:end-4)])
    end
    print(['../figures_Y2/stitching/stitch_assessment/redundant_gauge/' POBS_list{i}(1:end-4) '/' POBS_list{i}(1:end-4)],'-dpng','-r300')

    % plot temperature over full time interval
    figure(201); clf; hold on
    plot(t,(T1-mean(T1)),'linewidth',1)
    plot(t,(T2-mean(T2)),'linewidth',1)
    ylabel('T (C)')
    yyaxis right
    plot(t,(T1-mean(T1))-(T2-mean(T2)),'k','linewidth',1)
    ylabel('\DeltaT (C)')
    legend('Gauge 1','Gauge 2',['Gauge 1' char(8212) 'Gauge 2'],'location','northwest')
    datetick('x',6)
    title([shortname ' (Temperature)'])
    set(gca,'fontsize',14)
    box on; grid on

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print(['../figures_Y2/stitching/stitch_assessment/redundant_gauge/' POBS_list{i}(1:end-4) '/tempdif'],'-dpng','-r300')

    % zoom in on each calibration
    % for j=1:length(calInfoAll1.t0p)
    %     xlim([-2 2]+calInfoAll1.t0p(j)); ylim auto
    %     datetick('x','keeplimits')
    % 
    %     fh=gcf;
    %     fh.PaperUnits='inches';
    %     fh.PaperPosition=[0 0 11 8.5];
    %     print(['../figures_Y2/stitching/stitch_assessment/redundant_gauge/' POBS_list{i}(1:end-4) '/cal_' num2str(j)],'-djpeg','-r100')
    % end

    %% tidal filter
    % % calculate tides
    % [h1,th]=oceantide(stalat(i),stalon(i),-t(1),-t(end),60,'osu','datenum');
    % h1s=interp1(th,h1,t,'linear','extrap');
    % h1s=(h1s-mean(h1s))*100;
    % 
    % figure(22); clf; hold on
    % subplot(211); hold on
    % plot(t,(p1-mean(p1))-h1s,'linewidth',1)
    % ylabel('P (hPa)')
    % legend('de-tided APG','location','northwest')
    % datetick('x',6)
    % title({shortname; 'Gauge 1'})
    % set(gca,'fontsize',14)
    % box on; grid on
    % subplot(212); hold on
    % plot(t,(p2-mean(p2))-h1s,'linewidth',1)
    % ylabel('P (hPa)')
    % legend('de-tided APG','location','northwest')
    % datetick('x',6)
    % title('Gauge 2')
    % set(gca,'fontsize',14)
    % box on; grid on
    % 
    % % zoom in on each calibration
    % for j=1:length(calInfoAll1.t0p)
    %     subplot(211)
    %     xlim([-0.4 0.4]+calInfoAll1.t0p(j)); ylim auto
    %     datetick('x','keeplimits')
    %     subplot(212)
    %     xlim([-0.4 0.4]+calInfoAll1.t0p(j)); ylim auto
    %     datetick('x','keeplimits')
    % 
    %     fh=gcf;
    %     fh.PaperUnits='inches';
    %     fh.PaperPosition=[0 0 8.5 11];
    %     if ~exist(['../figures_Y2/stitching/stitch_assessment/tidal_model/' POBS_list{i}(1:end-4)],'dir')
    %         mkdir(['../figures_Y2/stitching/stitch_assessment/tidal_model/' POBS_list{i}(1:end-4)])
    %     end
    %     print(['../figures_Y2/stitching/stitch_assessment/tidal_model/' POBS_list{i}(1:end-4) '/cal_' num2str(j)],'-djpeg','-r100')
    % end
    
    %% depth-matched differences
    % % find GONDOR station at closest depth
    % [dmin,i_g]=min(abs(stadepth(i)-G2.stadepth));
    % 
    % % modifications
    % % if i==8 % depth too different to be useful
    % %     load([POBS_dir(k+2).folder '/' POBS_dir(k+2).name])
    % %     t3=data.ta(data.lNormState1);
    % %     p3=data.a1(data.lNormState1);
    % % elseif i==9 % depth too different to be useful
    % %     load([POBS_dir(k+1).folder '/' POBS_dir(k+1).name])
    % %     t3=data.ta(data.lNormState1);
    % %     p3=data.a1(data.lNormState1);
    % % elseif i==10 % depth too different to be useful
    % %     load([POBS_dir(k-2).folder '/' POBS_dir(k-2).name])
    % %     t3=data.ta(data.lNormState1);
    % %     p3=data.a1(data.lNormState1);
    % % elseif i_g==7 % this time series (KU22-PB) is very short
    % %     i_g=8;
    % %     t3=G2.ts{i_g}';
    % %     p3=G2.ps{i_g}';
    % % else
    % %     t3=G2.ts{i_g}';
    % %     p3=G2.ps{i_g}';
    % % end
    % 
    % t3=G2.ts{i_g}';
    % p3=G2.ps{i_g}';
    % 
    % % find temporal overlap
    % [~,ia,ib]=intersect(round(t*86400),round(t3*86400));
    % tt=t(ia);
    % p1=p1(ia);
    % p2=p2(ia);
    % p3=p3(ib);
    % 
    % figure(21); clf; hold on
    % subplot(211); hold on
    % plot(tt,(p1-mean(p1))-(p3-mean(p3)),'linewidth',1)
    % ylabel('P (hPa)')
    % legend(['APG' char(8212) 'depth match'],'location','northwest')
    % datetick('x',6)
    % title({shortname; 'Gauge 1'})
    % set(gca,'fontsize',14)
    % box on; grid on
    % subplot(212); hold on
    % plot(tt,(p2-mean(p2))-(p3-mean(p3)),'linewidth',1)
    % ylabel('P (hPa)')
    % legend(['APG' char(8212) 'depth match'],'location','northwest')
    % datetick('x',6)
    % title('Gauge 2')
    % set(gca,'fontsize',14)
    % box on; grid on
    % 
    % for j=1:length(calInfoAll1.t0p)
    %     if calInfoAll1.t0p(j)<tt(1) || calInfoAll1.t0p(j)>tt(end)
    %         continue
    %     end
    %     subplot(211)
    %     xlim([-0.4 0.4]+calInfoAll1.t0p(j)); ylim auto
    %     datetick('x','keeplimits')
    %     subplot(212)
    %     xlim([-0.4 0.4]+calInfoAll1.t0p(j)); ylim auto
    %     datetick('x','keeplimits')
    % 
    %     fh=gcf;
    %     fh.PaperUnits='inches';
    %     fh.PaperPosition=[0 0 8.5 11];
    %     if ~exist(['../figures_Y2/stitching/stitch_assessment/depth_matching/' POBS_list{i}(1:end-4)],'dir')
    %         mkdir(['../figures_Y2/stitching/stitch_assessment/depth_matching/' POBS_list{i}(1:end-4)])
    %     end
    %     print(['../figures_Y2/stitching/stitch_assessment/depth_matching/' POBS_list{i}(1:end-4) '/cal_' num2str(j)],'-djpeg','-r100')
    % end

    %% difference detided data
    % % remove tide
    % h1s=h1s(ia);
    % p1_cor=(p1-mean(p1))-h1s;
    % p2_cor=(p2-mean(p2))-h1s;
    % 
    % % modifications should match prior section
    % % if i==8
    % %     [h3,th]=oceantide(stalat(i+2),stalon(i+2),-tt(1),-tt(end),60,'osu','datenum');
    % % elseif i==9
    % %     [h3,th]=oceantide(stalat(i+1),stalon(i+1),-tt(1),-tt(end),60,'osu','datenum');
    % % elseif i==10
    % %     [h3,th]=oceantide(stalat(i-2),stalon(i-2),-tt(1),-tt(end),60,'osu','datenum');
    % % else
    % %     [h3,th]=oceantide(G2.stalat(i_g),G2.stalon(i_g),-tt(1),-tt(end),60,'osu','datenum');
    % % end
    % 
    % [h3,th]=oceantide(G2.stalat(i_g),G2.stalon(i_g),-tt(1),-tt(end),60,'osu','datenum');
    % 
    % h3s=interp1(th,h3,tt,'linear','extrap');
    % h3s=(h3s-mean(h3s))*100;
    % p3_cor=(p3-mean(p3))-h3s;
    % 
    % figure(24); clf; hold on
    % subplot(211); hold on
    % 
    % plot(tt,p1_cor-p3_cor,'linewidth',1)
    % ylabel('P (hPa)')
    % legend(['APG' char(8212) 'depth match (detided)'],'location','northwest')
    % datetick('x',6)
    % title({shortname; 'Gauge 1'})
    % set(gca,'fontsize',14)
    % box on; grid on
    % subplot(212); hold on
    % plot(tt,p2_cor-p3_cor,'linewidth',1)
    % ylabel('P (hPa)')
    % legend(['APG' char(8212) 'depth match (detided)'],'location','northwest')
    % datetick('x',6)
    % title('Gauge 2')
    % set(gca,'fontsize',14)
    % box on; grid on
    % 
    % for j=1:length(calInfoAll1.t0p)
    %     if calInfoAll1.t0p(j)<tt(1) || calInfoAll1.t0p(j)>tt(end)
    %         continue
    %     end
    %     subplot(211)
    %     xlim([-2 2]+calInfoAll1.t0p(j)); ylim auto
    %     datetick('x','keeplimits')
    %     subplot(212)
    %     xlim([-2 2]+calInfoAll1.t0p(j)); ylim auto
    %     datetick('x','keeplimits')
    % 
    %     fh=gcf;
    %     fh.PaperUnits='inches';
    %     fh.PaperPosition=[0 0 8.5 11];
    %     if ~exist(['../figures_Y2/stitching/stitch_assessment/tide&match/' POBS_list{i}(1:end-4)],'dir')
    %         mkdir(['../figures_Y2/stitching/stitch_assessment/tide&match/' POBS_list{i}(1:end-4)])
    %     end
    %     print(['../figures_Y2/stitching/stitch_assessment/tide&match/' POBS_list{i}(1:end-4) '/cal_' num2str(j)],'-djpeg','-r100')
    % end

end